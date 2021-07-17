#!/usr/bin/env python3

#=======================================================================
# Authors: Rossen Zhao
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
import tempdir
import sys
import io

from singlem.condense import Condenser

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data/condense')

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path

class Tests(unittest.TestCase):
    maxDiff = None
    
    def test_condense_small(self): #TODO add singlem packages, example input and output comparators
        with tempdir.in_tempdir():
            Condenser.condense(
                input_otu_table = os.path.join(
                    path_to_data, 'small_condense_input.csv'),
                singlem_packages = [
                    os.path.join(path_to_data, 'S2.1.ribosomal_protein_L2_rplB.gpkg.spkg'),
                    os.path.join(path_to_data, 'S2.10.ribosomal_protein_S7.gpkg.spkg'),
                    os.path.join(path_to_data, 'S2.11.ribosomal_protein_S10_rpsJ.gpkg.spkg')],
                trim_percent = 5,
                output_otu_table = 'small_condense_output.csv',
                krona = '/dev/null')
            self.assertListEqual(
                list(io.open('small_condense_output.csv')), 
                list(io.open(os.path.join(path_to_data, 'small_condense_output.csv'))))

if __name__ == "__main__":
    unittest.main()