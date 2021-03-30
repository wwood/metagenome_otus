#!/usr/bin/env python3

#=======================================================================
# Authors: Ben Woodcroft, Tim Lamberton.
#
# Unit tests.
#
# Copyright
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import unittest
import subprocess
import os.path
import tempfile
import tempdir
import extern
import sys
import json
import re

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.pipe import SearchPipe
from singlem.sequence_classes import SeqReader

class Tests(unittest.TestCase):
    headers = str.split('gene sample sequence num_hits coverage taxonomy')
    headers_with_extras = headers + str.split('read_names nucleotides_aligned taxonomy_by_known? read_unaligned_sequences')
    maxDiff = None
    two_packages = '%s %s' % (
        os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
        os.path.join(path_to_data, '4.12.22seqs.spkg'))

    def assertEqualOtuTable(self, expected_array, observed_string):
        observed_array = list([line.split("\t") for line in observed_string.split("\n")])
        if expected_array[-1] != ['']:
            expected_array.append([''])

        # make sure headers are OK
        self.assertEqual(expected_array[0], observed_array[0])

        # sort the rest of the table and compare that
        self.assertEqual(sorted(expected_array[1:]), sorted(observed_array[1:]))
    
    def test_fast_protein_package(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_fast_protein_package_prefilter(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --diamond_prefilter --otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))


    def test_fast_protein_package_diamond_package_assignment(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --diamond-package-assignment --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_fast_protein_package_diamond_package_assignment_paired_both_hit(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write(inseqs)
                n2.flush()

                cmd = "%s pipe --include-off-target-hits --forward %s --reverse %s --otu_table /dev/stdout --diamond-package-assignment --singlem_packages %s" % (
                    path_to_script, n.name, n2.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_fast_protein_package_diamond_package_assignment_paired_one_hits(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write('''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA''')
                n2.flush()

                cmd = "%s pipe --include-off-target-hits --forward %s --reverse %s --otu_table /dev/stdout --diamond-package-assignment --singlem_packages %s" % (
                    path_to_script, n.name, n2.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))
                cmd = "%s pipe --include-off-target-hits --reverse %s --forward %s --otu_table /dev/stdout --diamond-package-assignment --singlem_packages %s" % (
                    path_to_script, n.name, n2.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n2.name).replace('.fa',''),''))

    def test_fast_protein_package_prefilter_with_diamond_assignment(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --diamond_prefilter --otu_table /dev/stdout --singlem_packages %s --assignment-method diamond" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_fast_protein_package_prefilter_with_diamond_assignment_2_samples(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='sample1.fa') as n:
            n.write(inseqs)
            n.flush()

            with tempfile.NamedTemporaryFile(mode='w',suffix='sample2.fa') as n2:
                n2.write(inseqs)
                n2.flush()

                cmd = "%s pipe --include-off-target-hits --sequences %s %s --diamond_prefilter --otu_table /dev/stdout --singlem_packages %s --assignment-method diamond" % (
                    path_to_script, n.name, n2.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),'').replace(os.path.basename(n2.name).replace('.fa',''),''))

    def test_fast_protein_package_prefilter_with_diamond_assignment_paired(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --forward %s --reverse %s --diamond_prefilter --otu_table /dev/stdout --singlem_packages %s --assignment-method diamond" % (
                path_to_script, n.name, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))


    def test_fast_protein_package_prefilter_with_diamond_example_assignment(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	2524614704',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --forward %s --reverse %s --diamond_prefilter --otu_table /dev/stdout --singlem_packages %s --assignment-method diamond_example" % (
                path_to_script, n.name, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))
    
    def test_minimal(self):
        expected = [
            self.headers,
            ['S1.5.ribosomal_protein_L11_rplK','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s --debug pipe --include-off-target-hits --sequences %s/1_pipe/minimal.fa --otu_table /dev/stdout --threads 4" % (path_to_script,
                                                                                                    path_to_data)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))
        
    def test_minimal_prefilter(self):
        expected = [
            self.headers,
            ['S1.5.ribosomal_protein_L11_rplK','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s --debug pipe --include-off-target-hits --sequences %s/1_pipe/minimal.fa --diamond_prefilter --otu_table /dev/stdout --threads 4" % (path_to_script,
                                                                                                    path_to_data)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

    def test_insert(self):
        expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','2','4.95','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s --quiet pipe --include-off-target-hits --sequences %s/1_pipe/insert.fna --otu_table /dev/stdout --threads 4" % (path_to_script,
                                                                                                    path_to_data)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))
    
    def test_insert_prefilter(self):
        expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','2','4.95','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s --quiet pipe --include-off-target-hits --sequences %s/1_pipe/insert.fna --diamond_prefilter --otu_table /dev/stdout --threads 4" % (path_to_script,
                                                                                                    path_to_data)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

    def test_print_insert(self):
        expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','1','2.44','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['S1.5.ribosomal_protein_L11_rplK','insert','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTtttCAAGCAGGTGTG','1','2.51','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s --debug pipe --include-off-target-hits --sequences %s/1_pipe/insert.fna --otu_table /dev/stdout --threads 4 --include_inserts" % (path_to_script,
                                                                                                    path_to_data)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

    def test_known_tax_table(self):
        expected = [
            self.headers,
            ['4.12.22seqs','small',
             'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
             '4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
            ['4.11.22seqs','small',
             'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA',
             '2','4.88','Root; d__Bacteria; p__Firmicutes']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s --quiet pipe --include-off-target-hits --sequences %s/1_pipe/small.fa --otu_table /dev/stdout --threads 4 --singlem_packages %s" % (
            path_to_script,
            path_to_data,
            self.two_packages)
        self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

        expected = [self.headers,
                    ['4.12.22seqs','small',
                     'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                     '4','9.76','some1'],
                    ['4.11.22seqs','small',
                     'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA',
                     '2','4.88','Root; d__Bacteria; p__Firmicutes']]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        with tempfile.NamedTemporaryFile(mode='w',prefix='singlem_test_known') as t:
            t.write('\n'.join(["\t".join(x) for x in expected[:2]]))
            t.flush()

            cmd = "%s --quiet pipe --include-off-target-hits --sequences %s/1_pipe/small.fa --otu_table /dev/stdout --threads 4 --known_otu_tables %s --singlem_packages %s"\
                 % (path_to_script,
                    path_to_data,
                    t.name,
                    self.two_packages)
            self.assertEqual(exp, sorted(extern.run(cmd).split("\n")))

    def test_diamond_assign_taxonomy(self):
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fasta') as f:
            query = "\n".join(['>HWI-ST1243:156:D1K83ACXX:7:1109:18214:9910 1:N:0:TCCTGAGCCTAAGCCT',
                'GTTAAATTACAAATTCCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTGAACATCATGGGATTCTGTAAAGAGT',''])
            f.write(query)
            f.flush()

            cmd = "%s --debug pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --assignment_method diamond --threads 4" % (path_to_script,
                                                            f.name)

            expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK',os.path.basename(f.name)[:-6],'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','1','2.44','Root; d__Bacteria; p__Firmicutes; c__Bacilli_A; o__Thermoactinomycetales; f__Thermoactinomycetaceae']]
            expected = ["\t".join(x) for x in expected]+['']
            observed = extern.run(cmd).split("\n")
            r = re.compile('; c__.*') # Do not test beyond phylum level because updated diamond version change slightly.
            self.assertEqual([r.sub('',e) for e in expected], [r.sub('',e) for e in observed])

    def test_diamond_example_assign_taxonomy(self):
        expected = [self.headers,['S1.5.ribosomal_protein_L11_rplK','minimal','CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG','4','9.76','2513237297']
                    ]
        exp = sorted(["\t".join(x) for x in expected]+[''])

        cmd = "%s --debug pipe --include-off-target-hits --sequences %s/1_pipe/minimal.fa --otu_table /dev/stdout --threads 4 --assignment_method diamond_example" % (path_to_script,
                                                                                                    path_to_data)
        observed = sorted(extern.run(cmd).split("\n"))
        r = re.compile('\t.*?$') # Do not test the exact genome number because updated diamond version change this slightly.
        self.assertEqual([r.sub('',e) for e in exp], [r.sub('',e) for e in observed])

    def test_one_read_two_orfs_two_diamond_hits(self):
        # what a pain the real world is
        seq = '''>HWI-ST1240:128:C1DG3ACXX:7:2204:6599:65352 1:N:0:GTAGAGGATAGATCGC
ACCCACAGCTCGGGGTTGCCCTTGCCCGACCCCATGCGTGTCTCGGCGGGCTTCTGGTGACGGGCTTGTCCGGGAAGACGCGGATCCAGACCTTGCCTCCGCGCTTGACGTGCCGGGTCATCGCGATACGGGCCGCCTCGATCTGACGTGC
'''
        expected = [
            self.headers,
            ['S1.7.ribosomal_protein_L16_L10E_rplP		CGCGTCTTCCCGGACAAGCCCGTCACCAGAAGCCCGCCGAGACACGCATGGGGTCGGGCA	1	1.64	GCA_000949295.1']]
        exp = sorted(["\t".join(x) for x in expected]+[''])
        with tempfile.NamedTemporaryFile(mode='w',prefix='singlem_test',suffix='.fa') as t:
            t.write(seq)
            t.flush()
            cmd = "%s --quiet pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --threads 4 --assignment_method diamond_example" % (path_to_script,
                                                                                                    t.name)
            self.assertEqual(exp,
                             sorted(extern.run(cmd).
                                    replace(
                                        os.path.basename(t.name).replace('.fa',''),
                                        '').
                                    split("\n")))

    def test_jplace_output(self):
        expected_jpace = {'fields': ['classification',
                                     'distal_length',
                                     'edge_num',
                                     'like_weight_ratio',
                                     'likelihood',
                                     'pendant_length'],
                          'metadata': 'the_metadata',
                          'placements':
                          [{
                              'nm': [['CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                                      2]],
                              'p': [
                                  [
                                      "o__Bacillales",
                                      0.0874346630859,
                                      13,
                                      0.333350512423,
                                      -608.20180926,
                                      6.11351501465e-06
                                  ],
                                  [
                                      "o__Bacillales",
                                      0.0643521435547,
                                      14,
                                      0.333326884837,
                                      -608.201880142,
                                      6.11351501465e-06
                                  ],
                                  [
                                      "p__Firmicutes",
                                      5.97534179688e-06,
                                      15,
                                      0.33332260274,
                                      -608.201892989,
                                      6.11351501465e-06
                                  ]
                              ]}],
                          'tree': 'tree_thanks',
                          'version': 3}

        with tempdir.TempDir() as d:
            cmd = "%s pipe --include-off-target-hits --sequences %s --otu_table /dev/null --output_jplace %s"\
                  " --singlem_packages %s" % (
                      path_to_script,
                      os.path.join(path_to_data,'1_pipe','jplace_test.fna'),
                      os.path.join(d, "my_jplace"),
                      os.path.join(path_to_data,'4.12.22seqs.spkg'))
            extern.run(cmd)
            jplace_path = os.path.join(d, 'my_jplace_jplace_test_4.12.22seqs.jplace')
            with open(jplace_path) as f:
                j = json.load(f)
            j['tree'] = 'tree_thanks'
            j['metadata'] = 'the_metadata'
            self.assertEqual(expected_jpace, j)

            # Make sure the guppy sing does not croak
            extern.run("guppy sing -o /dev/null '%s'" % jplace_path)

    def test_nucleotide_package(self):
        expected = [
            "\t".join(self.headers),
            '61_otus.v3		GGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGACTGACGCTGATGTGCGAAAGCG	1	2.56	Root; k__Bacteria; p__Proteobacteria',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483 1:N:0:AAGAGGCAAAGGAGTA
GATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGGATTAGATACCCTGGTAGT
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))
            
    def test_nucleotide_package_prefilter(self):
        """ correct behaviour is to fail, as DIAMOND does not have a blastn capability """
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483 1:N:0:AAGAGGCAAAGGAGTA
GATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGGATTAGATACCCTGGTAGT
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            try:
                cmd = "%s pipe --include-off-target-hits --sequences %s --diamond_prefilter --otu_table /dev/stdout --singlem_packages %s" % (
                    path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'))
                extern.run(cmd)
                self.fail() # this is meant to fail
            except:
                pass
                

    def test_revcom_nucleotide_package(self):
        expected = [
            "\t".join(self.headers),
            '61_otus.v3		GGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGACTGACGCTGATGTGCGAAAGCG	1	2.56	Root; k__Bacteria; p__Proteobacteria',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483_revcom
ACTACCAGGGTATCTAATCCTGTTTGATCCCCACGCTTTCGCACATCAGCGTCAGTTACAGACCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCATATC
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))

    def test_fwd_and_revcom_nucleotide_package(self):
        expected = [
            "\t".join(self.headers),
            '61_otus.v3		GGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGACTGACGCTGATGTGCGAAAGCG	2	5.13	Root; k__Bacteria; p__Proteobacteria',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483 1:N:0:AAGAGGCAAAGGAGTA
GATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGGATTAGATACCCTGGTAGT
>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483_revcom
ACTACCAGGGTATCTAATCCTGTTTGATCCCCACGCTTTCGCACATCAGCGTCAGTTACAGACCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCATATC
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))

    def test_two_nucleotide_packages(self):
        expected = [
            "\t".join(self.headers),
            '61_otus.v3		GGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGACTGACGCTGATGTGCGAAAGCG	2	5.13	Root; k__Bacteria; p__Proteobacteria',
            '61_otus.second.v3		TTAGGTAGTTGCTGGGGTAACGTCCCAACAAGCCGATAATCGGTACGGGTTGTGAGAGCA	1	1.66	Root; k__Archaea; p__Euryarchaeota',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483 1:N:0:AAGAGGCAAAGGAGTA
GATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGGATTAGATACCCTGGTAGT
>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483_revcom
ACTACCAGGGTATCTAATCCTGTTTGATCCCCACGCTTTCGCACATCAGCGTCAGTTACAGACCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCATATC
>NS500333:10:H0V2GAGXX:2:13211:8623:16289 1:N:0:GATCAG
ATTAGGTAGTTGCTGGGGTAACGTCCCAACAAGCCGATAATCGGTACGGGTTGTGAGAGCAAGAGCCCGGAGATGGATTCTGAGACACGAATCCAGGTCCTACGGGGCGCAGCAGGCGCGAAAACTTTACACTGCGCGAAAGCGCGATA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --singlem_packages %s %s" % (
                path_to_script,
                n.name,
                os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'),
                os.path.join(path_to_data,'second_packge.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_no_taxonomy(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --singlem_packages %s --no_assign_taxonomy" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))

    def test_known_sequence_taxonomy(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	2	4.88	mytax; yeh',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
>another
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            with tempfile.NamedTemporaryFile(mode='w',) as taxf:
                taxf.write("HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482\tmytax; yeh\n")
                taxf.write("another\tmytax; yeh; 2\n")
                taxf.flush()

                cmd = "%s pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --singlem_packages %s "\
                      "--no_assign_taxonomy --known_sequence_taxonomy %s"% (
                          path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'),
                          taxf.name)
                self.assertEqual(expected,
                                 extern.run(cmd).replace(
                                     os.path.basename(n.name).replace('.fa',''),
                                     '').split("\n"))

    def test_sample_name_strange_characters(self):
        expected = [self.headers,
                    ['4.12.22seqs','contigs.fasta.metabat-bins-_-t20_--superspecific.8',
                     'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                     '4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.22seqs','contigs.fasta.metabat-bins-_-t20_--superspecific.8',
                     'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA',
                     '2','4.88','Root; d__Bacteria; p__Firmicutes'],
                    ['4.12.22seqs','contigs.fasta.metabat-bins-_-t20_--superspecific.9',
                     'CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG',
                     '4','9.76','Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales'],
                    ['4.11.22seqs','contigs.fasta.metabat-bins-_-t20_--superspecific.9',
                     'TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA',
                     '2','4.88','Root; d__Bacteria; p__Firmicutes']]
        cmd = "%s --quiet pipe --include-off-target-hits --sequences "\
        "%s/1_pipe/contigs.fasta.metabat-bins-_-t20_--superspecific.9.fa "\
        "%s/1_pipe/contigs.fasta.metabat-bins-_-t20_--superspecific.8.fa "\
        "--otu_table /dev/stdout --threads 4 --singlem_packages %s" %(
            path_to_script,
            path_to_data,
            path_to_data,
            self.two_packages)
        self.assertEqualOtuTable(expected, extern.run(cmd))

    def test_archive_otu_groopm_compatibility(self):
        """This tests for API stability, where the API is used by GroopM 2.0"""
        expected = [('contig_1', '4.11.22seqs', 'Root; d__Bacteria; p__Firmicutes')]

        inseqs = '''>contig_1 abcd
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --archive_otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'))

            j = json.loads(extern.run(cmd))
            fields = j['fields']
            data = j['otus']
            self.assertEqual(expected,
                             [(name, row[fields.index('gene')], row[fields.index('taxonomy')]) for row in data for name in row[fields.index('read_names')]]
                            )

    def test_protein_package_non60_length(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTT	1	1.41	Root; d__Bacteria; p__Firmicutes',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.length30.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_nucleotide_package_non60_length(self):
        expected = [
            "\t".join(self.headers),
            '61_otus.v3		GGAGGAACAC	1	1.10	Root; k__Bacteria; p__Proteobacteria',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1105:6981:63483 1:N:0:AAGAGGCAAAGGAGTA
GATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAGGATTAGATACCCTGGTAGT
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.length10.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))

    def test_nucleotide_matches_forward_and_reverse(self):
        '''Not the best test since the read does not actually match the alignment
        region, but at least no error should be thrown.'''
        expected = [
            "\t".join(self.headers),
            '']
        inseqs = '''>seq18975201
GTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACGCCATGGAGTCGAGTTGCAGACTCCAATCCGAACTGGGGCCGGTTTTTATGGATTGGCTTCCCCTCGCGGGTTCGCGACCCTTTGTACCGGCCATTGTAACACGTGTGTAGC
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --include-off-target-hits --sequences %s --otu_table /dev/stdout --singlem_packages %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'61_otus.v3.gpkg.spkg'))
            self.assertEqual(expected,
                             extern.run(cmd).replace(
                                 os.path.basename(n.name).replace('.fa',''),
                                 '').split("\n"))

    def test_paired_reads_hello_world(self):
        # Reads should be merged
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        inseqs_reverse = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT
''' # reverse complement of the forward, so should collapse.
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write(inseqs_reverse)
                n2.flush()

                cmd = "{} pipe --include-off-target-hits --sequences {} --otu_table /dev/stdout --singlem_packages {} --reverse {}".format(
                    path_to_script,
                    n.name,
                    os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'),
                    n2.name)
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_paired_reads_one_read_each(self):
        # Reads should be merged
        expected = [
            "\t".join(self.headers_with_extras),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	2	4.88	Root; d__Bacteria; p__Firmicutes	HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 seq2	60 60	False	ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
>seq2
AAAAAAAAAAAAAAAAA
'''
        inseqs_reverse = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
AAAAAAAAAAAAAAAAA
>seq2
TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT
''' # reverse complement of the forward, so should collapse.
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write(inseqs_reverse)
                n2.flush()

                cmd = "{} pipe --include-off-target-hits --sequences {} --otu_table /dev/stdout --singlem_packages {} --reverse {} --output_extras".format(
                    path_to_script,
                    n.name,
                    os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'),
                    n2.name)
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_paired_reads_two_reads_each(self):
        # Reads should be merged
        expected = [
            "\t".join(self.headers_with_extras),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	3	7.32	Root; d__Bacteria; p__Firmicutes	HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 seq2 seq3	60 60 60	False	ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
>seq2
AAAAAAAAAAAAAAAAA
>seq3
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        inseqs_reverse = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
AAAAAAAAAAAAAAAAA
>seq2
TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT
>seq3
TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT
''' # reverse complement of the forward, so should collapse.
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write(inseqs_reverse)
                n2.flush()

                cmd = "{} pipe --include-off-target-hits --sequences {} --otu_table /dev/stdout --singlem_packages {} --reverse {} --output_extras".format(
                    path_to_script,
                    n.name,
                    os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'),
                    n2.name)
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))


    def test_paired_reads_one_read_each_diamond(self):
        # Reads should be merged
        expected = [
            "\t".join(self.headers_with_extras),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	2	4.88	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179	HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 seq2	60 60	False	ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
>seq2
AAAAAAAAAAAAAAAAA
'''
        inseqs_reverse = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
AAAAAAAAAAAAAAAAA
>seq2
TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT
''' # reverse complement of the forward, so should collapse.
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write(inseqs_reverse)
                n2.flush()

                cmd = "{} pipe --include-off-target-hits --sequences {} --otu_table /dev/stdout --singlem_packages {} --reverse {} --output_extras --assignment_method diamond".format(
                    path_to_script,
                    n.name,
                    os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'),
                    n2.name)
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))


    def test_paired_reads_one_read_each_diamond_example(self):
        # Reads should be merged
        expected = [
            "\t".join(self.headers_with_extras),
            '4.11.22seqs		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	2	4.88	2524614704	HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 seq2	60 60	False	ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
>seq2
AAAAAAAAAAAAAAAAA
'''
        inseqs_reverse = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
AAAAAAAAAAAAAAAAA
>seq2
TTCAGCTGCACGACGTACCATAGTGTTTTTGTATACTTTATACTCAACACCAGCTTCACGTAATTGTGAACGTAAGTCAGTAACTTCAGCTACTGTTAAT
''' # reverse complement of the forward, so should collapse.
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()
            with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n2:
                n2.write(inseqs_reverse)
                n2.flush()

                cmd = "{} pipe --include-off-target-hits --sequences {} --otu_table /dev/stdout --singlem_packages {} --reverse {} --output_extras --assignment_method diamond_example".format(
                    path_to_script,
                    n.name,
                    os.path.join(path_to_data,'4.11.22seqs.gpkg.spkg'),
                    n2.name)
                self.assertEqualOtuTable(
                    list([line.split("\t") for line in expected]),
                    extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_two_orfs_in_same_read(self):
        # Read finds 2 ORFs, pretty rare for reads (but happens for genomes
        # more frequently)
        expected = [
            "\t".join(self.headers_with_extras),
            'S1.12.ribosomal_protein_S12_S23		CGTGGTGTCTGCACCCGGGTGTACACCACCACCCGAAGA---------AGCCGAACTCGG	1	1.50	Root; d__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Propionibacteriales	A00178:38:H5NYYDSXX:2:1552:32524:1517	51	False	CGGGATGTAGGCAGTGACCTCCACGCCTGAGGAGAGCCGGACGCGTGCGACCTTGCGCAACGCCGAGTTCGGCTTCTTCGGGTGGTGGTGTACACCCGGGTGCAGACACCACGGCGCTGGGGCGAACCCTTGAGCGCAGGGGTGTTGGTCT',
            '']
        inseqs = '''>A00178:38:H5NYYDSXX:2:1552:32524:1517 1:N:0:CAACGGA+ATCCGTT
CGGGATGTAGGCAGTGACCTCCACGCCTGAGGAGAGCCGGACGCGTGCGACCTTGCGCAACGCCGAGTTCGGCTTCTTCGGGTGGTGGTGTACACCCGGGTGCAGACACCACGGCGCTGGGGCGAACCCTTGAGCGCAGGGGTGTTGGTCT
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "{} pipe --include-off-target-hits --sequences {} --otu_table /dev/stdout --singlem_packages {} --output_extras".format(
                path_to_script,
                n.name,
                os.path.join(path_to_data, 'S1.12.ribosomal_protein_S12_S23.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))

    def test_split_genes(self):
        # 2 ORFs found causing GraftM to do its "split" thing
        expected = [
            "\t".join(self.headers),
            "S1.2.ribosomal_protein_L3_rplC	aa_orf_split_bug	GTTGACGTGGCGGCCATCACAAAGGGCAAGGGATGGCAGGGCGTCCTGAAGCGGTGGAAC	1	1.05	Root; d__Archaea; p__Crenarchaeota; c__Nitrososphaeria; o__Nitrososphaerales; f__Nitrosopumilaceae; g__Nitrosopumilus",
            '']
        cmd = "{} pipe --include-off-target-hits --sequences {} --otu_table /dev/stdout --singlem_packages {}".format(
            path_to_script,
            os.path.join(path_to_data, 'aa_orf_split_bug.fna'),
            os.path.join(path_to_data, 'S1.2.ribosomal_protein_L3_rplC.gpkg.spkg'))
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_no_results(self):
        # ORFs but no hits at all 
        expected = [
            "\t".join(self.headers),
            '']
        cmd = "{} pipe --include-off-target-hits --sequences {} --otu_table /dev/stdout --singlem_packages {}".format(
            path_to_script,
            os.path.join(path_to_data, 'random.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'))
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_no_results_diamond_prefilter(self):
        # ORFs but no hits at all 
        expected = [
            "\t".join(self.headers),
            '']
        cmd = "{} pipe --include-off-target-hits --sequences {} --otu_table /dev/stdout --singlem_packages {} --diamond-prefilter".format(
            path_to_script,
            os.path.join(path_to_data, 'random.fna'),
            os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'))
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_genome(self):
        expected = [
            "\t".join(self.headers),
            'S1.13.ribosomal_protein_S15P_S13e	uap2	AAGGATTTGAGTGCAAAAAGAGGACTCGATTTTACAGAGGCAAAGATAAGAAAACTTGGA	1	1.15	Root; d__Archaea; p__Halobacterota; c__Methanomicrobia; o__Methanomicrobiales; f__Methanocullaceae; g__Methanoculleus',
            '']

        cmd = '{} pipe --include-off-target-hits --genome-fasta-files {}//uap2.fna --singlem-packages {}/S1.13.chainsaw.gpkg.spkg --otu_table /dev/stdout --assignment-method diamond --diamond-prefilter --assignment-method diamond'.format(
            path_to_script,
            path_to_data,
            path_to_data,
        )
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_prefilter_piped_in_input(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs	63	TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	2	4.88	Root; d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Lachnospiraceae_bacterium_NK4A179]; s__Lachnospiraceae_bacterium_NK4A179',
            '']

        cmd = '{} pipe --include-off-target-hits --forward <(cat {}/1_pipe/small.fa) --singlem-packages {}/4.11.22seqs.gpkg.spkg --otu_table /dev/stdout --assignment-method diamond --diamond-prefilter --assignment-method diamond'.format(
            path_to_script,
            path_to_data,
            path_to_data,
        )
        self.assertEqualOtuTable(
            list([line.split("\t") for line in expected]),
            extern.run(cmd))

    def test_keep_target_hits(self):
        expected = [
            "\t".join(self.headers),
            '4.11.22seqs.v3_bacteria		TTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTA	1	2.44	Root; d__Bacteria; p__Firmicutes',
            '']
        inseqs = '''>HWI-ST1243:156:D1K83ACXX:7:1106:18671:79482 1:N:0:TAAGGCGACTAAGCCT
ATTAACAGTAGCTGAAGTTACTGACTTACGTTCACAATTACGTGAAGCTGGTGTTGAGTATAAAGTATACAAAAACACTATGGTACGTCGTGCAGCTGAA
'''
        with tempfile.NamedTemporaryFile(mode='w',suffix='.fa') as n:
            n.write(inseqs)
            n.flush()

            cmd = "%s pipe --sequences %s --otu_table /dev/stdout --singlem_packages %s %s" % (
                path_to_script, n.name, os.path.join(path_to_data,'4.11.22seqs.v3_bacteria_targetted.gpkg.spkg'), os.path.join(path_to_data,'4.11.22seqs.v3_archaea_targetted.gpkg.spkg'))
            self.assertEqualOtuTable(
                list([line.split("\t") for line in expected]),
                extern.run(cmd).replace(os.path.basename(n.name).replace('.fa',''),''))




if __name__ == "__main__":
    unittest.main()
