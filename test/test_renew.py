
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
import os.path
import extern
import sys

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','singlem')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
from singlem.renew import Renew

class Tests(unittest.TestCase):
    headers = str.split('gene sample sequence num_hits coverage taxonomy')
    maxDiff = None
    two_packages = '%s %s' % (
        os.path.join(path_to_data, '4.11.22seqs.gpkg.spkg'),
        os.path.join(path_to_data, '4.12.22seqs.spkg'))

    def test_hello_word_cmdline(self):
        cmd = "{} renew --include_off_target_hits --input_archive_otu_table {}/small.otu_table.json --singlem_package {}/4.12.22seqs.spkg --sequences {}/1_pipe/small.fa --assignment_method diamond --otu_table /dev/stdout".format(
            path_to_script,
            path_to_data,
            path_to_data,
            path_to_data)
        output = extern.run(cmd)
        expected = [
            self.headers,
            ["4.12.22seqs","small","CCTGCAGGTAAAGCGAATCCAGCACCACCAGTTGGTCCAGCATTAGGTCAAGCAGGTGTG","4","9.76","Root; d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Gracilibacillus; s__Gracilibacillus_lacisalsi"]
        ]
        self.assertEqualOtuTable(expected, output)

    def assertEqualOtuTable(self, expected_array, observed_string):
        observed_array = list([line.split("\t") for line in observed_string.split("\n")])
        if expected_array[-1] != ['']:
            expected_array.append([''])

        # make sure headers are OK
        self.assertEqual(expected_array[0], observed_array[0])

        # sort the rest of the table and compare that
        self.assertEqual(sorted(expected_array[1:]), sorted(observed_array[1:]))


if __name__ == "__main__":
    unittest.main()
