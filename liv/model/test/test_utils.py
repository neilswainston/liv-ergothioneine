'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
import os
import random
import unittest

from cobra.io import read_sbml_model
from synbiochem.utils.chem_utils import get_molecular_mass

from liv.model.utils import get_mw


class Test(unittest.TestCase):
    '''Test class for ice_utils.'''

    def test_get_mw(self):
        '''Tests get_mw method.'''
        curr_dir = os.path.dirname(os.path.realpath(__file__))

        model = read_sbml_model(
            os.path.join(curr_dir, '../../../data/models/yeastGEM.xml'))

        tests = 0

        while tests < 10:
            met = random.choice(model.metabolites)
            # met = model.metabolites.get_by_id('s_3071__91__lp__93__')
            print(met.id, met.name)

            # Get existing formula and mw:
            formula = met.formula

            if formula:
                existing_mw = get_molecular_mass(formula, r_mass=2**16)

                # Unset formula:
                met.formula = None

                calc_mw = get_mw(model, met.id)

                if calc_mw:
                    self.assertAlmostEqual(existing_mw, calc_mw, 0)
                    tests += 1

                # Reset formula:
                met.formula = formula


if __name__ == "__main__":
    unittest.main()
