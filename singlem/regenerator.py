import logging
import os

import extern
import dendropy
from graftm.graftm_package import GraftMPackage

from .graftm_result import GraftMResult
from .singlem_package import SingleMPackageVersion3, SingleMPackageVersion2, SingleMPackage
from .sequence_classes import SeqReader
from .dereplicator import Dereplicator
from .sequence_extractor import SequenceExtractor


class Regenerator:
    def regenerate(self, **kwargs):
        input_singlem_package = kwargs.pop('input_singlem_package')
        output_singlem_package = kwargs.pop('output_singlem_package')
        working_directory = kwargs.pop('working_directory')
        euk_sequences = kwargs.pop('euk_sequences')
        euk_taxonomy = kwargs.pop('euk_taxonomy')
        intermediate_archaea_graftm_package = kwargs.pop('intermediate_archaea_graftm_package')
        intermediate_bacteria_graftm_package = kwargs.pop('intermediate_bacteria_graftm_package')
        input_taxonomy = kwargs.pop('input_taxonomy')
        min_aligned_percent = kwargs.pop('min_aligned_percent')

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        original_pkg = SingleMPackage.acquire(input_singlem_package)
        original_hmm_path = original_pkg.hmm_path()
        basename = original_pkg.graftm_package_basename()

        # Run GraftM on the euk sequences with the bacterial set
        euk_graftm_output = os.path.join(working_directory,
                                         "%s-euk_graftm" % basename)
        cmd = "graftM graft --graftm_package '%s' --search_and_align_only --forward '%s' --output %s --force" % (
            original_pkg.graftm_package_path(),
            euk_sequences,
            euk_graftm_output)
        extern.run(cmd)

        # Extract hit sequences from that set
        euk_result = GraftMResult(euk_graftm_output, False)
        hit_paths = euk_result.unaligned_sequence_paths(require_hits=True)
        if len(hit_paths) != 1: raise Exception(
                "Unexpected number of hits against euk in graftm")
        euk_hits_path = next(iter(hit_paths.values())) #i.e. first

        # Concatenate euk, archaea and bacterial sequences
        archaeal_intermediate_pkg = GraftMPackage.acquire(
            intermediate_archaea_graftm_package)
        bacterial_intermediate_pkg = GraftMPackage.acquire(
            intermediate_bacteria_graftm_package)
        num_euk_hits = 0
        final_sequences_path = os.path.join(working_directory,
                                            "%s_final_sequences.faa" % basename)

        with open(final_sequences_path, 'w') as final_seqs_fp:
            with open(euk_hits_path) as euk_seqs_fp:
                for name, seq, _ in SeqReader().readfq(euk_seqs_fp):
                    if name.find('_split_') == -1:
                        num_euk_hits += 1
                        final_seqs_fp.write(">%s\n%s\n" % (name, seq))
            logging.info("Found %i eukaryotic sequences to include in the package" % \
                         num_euk_hits)

            for gpkg in [archaeal_intermediate_pkg, bacterial_intermediate_pkg]:
                num_total = 0
                num_written = 0
                with open(gpkg.unaligned_sequence_database_path()) as seqs:
                    for name, seq, _ in SeqReader().readfq(seqs):
                        num_total += 1
                        # if name in species_dereplicated_ids:
                        final_seqs_fp.write(">%s\n%s\n" % (name, seq))
                        num_written += 1
                logging.info(
                    "Of %i sequences in gpkg %s, %i species-dereplicated were included in the final package." %(
                        num_total, gpkg, num_written))

        # Concatenate euk and input taxonomy
        final_taxonomy_file = os.path.join(working_directory,
                                            "%s_final_taxonomy.csv" % basename)
        extern.run("cat %s %s > %s" % (
            euk_taxonomy, input_taxonomy, final_taxonomy_file))

        # Run graftm create to get the final package
        final_gpkg = os.path.join(working_directory,
                                  "%s_final.gpkg" % basename)
        cmd = "graftM create --force --min_aligned_percent %s --sequences %s --taxonomy %s --search_hmm_files %s %s --hmm %s --output %s" % (
            min_aligned_percent,
            final_sequences_path,
            final_taxonomy_file,
            ' '.join(archaeal_intermediate_pkg.search_hmm_paths()),
            ' '.join(bacterial_intermediate_pkg.search_hmm_paths()),
            original_hmm_path,
            final_gpkg)
        try:
            extern.run(cmd)
        except Exception:
            logging.info("Automatically retrying graftM after taxit create failure")
            
            rerooted_tree = "graftm_create_tree." + basename.split(".")[0] + ".tree"
            
            cmd = "graftM create --force --min_aligned_percent %s --sequences %s --taxonomy %s --search_hmm_files %s %s --hmm %s --output %s --rerooted_tree %s" % (
            min_aligned_percent,
            final_sequences_path,
            final_taxonomy_file,
            ' '.join(archaeal_intermediate_pkg.search_hmm_paths()),
            ' '.join(bacterial_intermediate_pkg.search_hmm_paths()),
            original_hmm_path,
            final_gpkg,
            rerooted_tree)
            extern.run(cmd)

        ##############################################################################
        # Remove sequences from the diamond DB that are not in the tree i.e.
        # those that are exact duplicates, so that the diamond_example hits are
        # always in the tree.
        # Read the list of IDs in the tree with dendropy
        final_gpkg_object = GraftMPackage.acquire(final_gpkg)
        unaligned_seqs = final_gpkg_object.unaligned_sequence_database_path()
        tree = dendropy.Tree.get(path=final_gpkg_object.reference_package_tree_path(),
                                 schema='newick')
        leaf_names = [l.taxon.label.replace(' ','_') for l in tree.leaf_node_iter()]
        logging.debug("Read in final tree with %i leaves" % len(leaf_names))

        # Extract out of the sequences file in the graftm package
        final_seqs = SequenceExtractor().extract_and_read(
            leaf_names, unaligned_seqs)
        if len(final_seqs) != len(leaf_names):
            raise Exception("Do not appear to have extracted the expected number of sequences from the unaligned fastat file")

        # Write the reads into sequences file in place
        with open(unaligned_seqs, 'w') as f:
            for s in final_seqs:
                f.write(">%s\n" % s.name)
                f.write(s.seq)
                f.write("\n")

        # Regenerate the diamond DB
        final_gpkg_object.create_diamond_db()

        ##############################################################################
        # Run singlem create to put the final package together
        if original_pkg.version == 2:
            SingleMPackageVersion2.compile(
                output_singlem_package,
                final_gpkg,
                original_pkg.singlem_position(),
                original_pkg.window_size())
        elif original_pkg.version == 3:
            SingleMPackageVersion3.compile(
                output_singlem_package,
                final_gpkg,
                original_pkg.singlem_position(),
                original_pkg.window_size(),
                original_pkg.target_domains(),
                original_pkg.gene_description())
        logging.info("SingleM package generated.")


