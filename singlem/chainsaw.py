import logging
import os
import shutil
from singlem.sequence_classes import SeqReader
from .singlem_package import SingleMPackage
from graftm.graftm_package import GraftMPackage

class Chainsaw:
    @staticmethod
    def chainsaw(**kwargs):
        '''Renew an OTU table'''
        input_singlem_package_path = kwargs.pop('input_singlem_package_path')
        output_singlem_package_path = kwargs.pop('output_singlem_package_path')
        sequence_prefix = kwargs.pop('sequence_prefix')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        # Open the input package
        input_spkg = SingleMPackage.acquire(input_singlem_package_path)

        # Ensure v2
        if input_spkg.version != 2:
            raise Exception("Only works with v2 packages for the moment.")
        
        # mkdir output package folder
        os.mkdir(output_singlem_package_path)
        shutil.copyfile(
            os.path.join(input_singlem_package_path,'CONTENTS.json'),
            os.path.join(output_singlem_package_path,'CONTENTS.json'))
        # mkdir graftm folder
        graftm_path = os.path.join(output_singlem_package_path, os.path.basename(input_spkg.graftm_package_path()))
        os.mkdir(graftm_path)
        # mkdir refpkg folder
        refpkg_path = os.path.join(graftm_path, os.path.basename(input_spkg.graftm_package().reference_package_path()))
        os.mkdir(refpkg_path)
        shutil.copyfile(
            os.path.join(input_spkg.graftm_package().reference_package_path(), 'CONTENTS.json'),
            os.path.join(refpkg_path,'CONTENTS.json'))
        # Copy the taxonomy file and the contents file
        taxonomy = input_spkg.graftm_package().taxtastic_taxonomy_path()
        shutil.copyfile(taxonomy, os.path.join(refpkg_path, os.path.basename(taxonomy)))
        for f in [
            os.path.join(input_spkg.graftm_package_path(),'CONTENTS.json'),
            input_spkg.graftm_package().alignment_hmm_path()]:

            shutil.copyfile(f, os.path.join(graftm_path, os.path.basename(f)))

        for f in input_spkg.graftm_package().search_hmm_paths():
            shutil.copyfile(f, os.path.join(graftm_path, os.path.basename(f)))

        # Change the sequences adding the prefix
        unaligned_basename = os.path.basename(input_spkg.graftm_package().unaligned_sequence_database_path())
        with open(os.path.join(graftm_path,unaligned_basename),'w') as out:
            for (name, seq, _) in SeqReader().readfq(open(input_spkg.graftm_package().unaligned_sequence_database_path())):
                out.write(">{}{}\n{}\n".format(sequence_prefix,name,seq))

        # Change the seqinfo adding the prefix
        seqinfo_basename = os.path.basename(input_spkg.graftm_package().taxtastic_seqinfo_path())
        with open(os.path.join(refpkg_path,seqinfo_basename),'w') as out:
            first = True
            for line in open(input_spkg.graftm_package().taxtastic_seqinfo_path()):
                if first:
                    out.write(line)
                    first = False
                else:
                    out.write("{}{}".format(sequence_prefix,line))

        logging.info("Remaking DIAMOND db ..")
        new_graftm = GraftMPackage.acquire(graftm_path)
        dmnd = new_graftm.create_diamond_db()

        logging.info("Chainsaw stopping.")