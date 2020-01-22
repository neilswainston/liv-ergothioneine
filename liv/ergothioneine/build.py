'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
import os
import sys
from cobra.io import read_sbml_model, write_sbml_model

from liv.model.utils import add_creator, add_met, add_reaction, makedirs


def build(model):
    '''Build model.'''

    # Add ergothioneine pathway:
    _add_pathway(model)

    # Add transport / media:
    _add_transport(model)

    return model


def _add_pathway(model):
    '''Add ergothioneine pathway.'''
    add_met(model, 'hercynine_c', 'hercynine', 'C9H15N3O2', 'c')
    add_met(model, 'herncystsulfox_c', 'hercynylcysteine sulfoxide',
            'C12H20N4O5S', 'c')
    add_met(model, 'sulfhercyn_c', '2-sulfenohercynine', 'C9H15N3O3S', 'c')
    add_met(model, 'ergothioneine_c', 'ergothioneine', 'C9H15N3O2S', 'c')

    # L-histidine + 3 S-adenosyl-L-methionine -->
    # hercynine + S-adenosyl-L-homocysteine + 3 H+
    reac_str = '''
        s_1006[c] + 3 s_1416[c] -->
        hercynine_c + 3 s_1413[c] + 3 s_0794[c]
        '''
    add_reaction(model, 'Egt1a', 'Egt1a', reac_str, 'Egt1')

    # hercynine + L-cysteine + oxygen --> hercynylcysteine sulfoxide + H20
    reac_str = '''
        hercynine_c + s_0981[c] + s_1275[c] -->
        herncystsulfox_c + s_0803[c]
        '''
    add_reaction(model, 'Egt1b', 'Egt1b', reac_str, 'Egt1')

    # hercynylcysteine sulfoxide + H20 -->
    # 2-sulfenohercynine + ammonium pyruvate
    reac_str = '''
        herncystsulfox_c + s_0803[c] -->
        sulfhercyn_c + s_0419[c] + s_4184[c]
        '''
    add_reaction(model, 'Egt2', 'Egt2', reac_str, 'Egt2')

    # 2-sulfenohercynine + electron donor (glutathione) ->
    # ergothioneine + electron_acceptor (glutathione disulfide) + H20
    reac_str = '''
        sulfhercyn_c + 2 s_0750[c] -->
        ergothioneine_c + s_0754[c] + s_0803[c]
        '''
    add_reaction(model, 'ergothioneine_synthesis', 'ergothioneine synthesis',
                 reac_str)

    # ergothioneine ->
    reac_str = '''
        ergothioneine_c -->
        '''
    add_reaction(model, 'ergothioneine_sink', 'ergothioneine sink',
                 reac_str, check=False)


def _add_transport(model):
    '''Add transport bounds.'''
    duration = 84  # hr

    # Glucose transport:
    # 175g / 180.2 g/mol / 84h * 1000 = mmol/h
    glucose_mass = 180.2
    glucose_uptake_min = (175 - 3.5) / glucose_mass / duration * 1000
    glucose_uptake_max = (175 + 3.5) / glucose_mass / duration * 1000
    model.reactions.get_by_id('r_1714').lower_bound = -glucose_uptake_max
    model.reactions.get_by_id('r_1714').upper_bound = -glucose_uptake_min

    # ATP maintenance:
    atp_main_fact = 0.7

    model.reactions.get_by_id('r_4046').lower_bound = \
        atp_main_fact * glucose_uptake_min
    model.reactions.get_by_id('r_4046').upper_bound = \
        atp_main_fact * glucose_uptake_max

    # Ergothioneine:
    ergoth_mass = 229.3
    ergoth_prod_min = (0.598 - 0.018) / ergoth_mass / duration * 1000
    ergoth_prod_max = (0.598 + 0.018) / ergoth_mass / duration * 1000
    model.reactions.get_by_id('ergothioneine_sink').lower_bound = \
        ergoth_prod_min
    model.reactions.get_by_id('ergothioneine_sink').upper_bound = \
        ergoth_prod_max

    # L-histidine transport:
    # model.reactions.get_by_id('r_1893').lower_bound = -1

    # L-cysteine transport:
    # model.reactions.get_by_id('r_1883').lower_bound = -1


def main(args):
    '''main method.'''
    model = read_sbml_model(args[0])

    # Create updated model:
    build(model)
    add_creator(model, *args[3:7])

    # Write updated model:
    out_filename = os.path.join(args[1], '%s.xml' % args[2])
    makedirs(out_filename)
    write_sbml_model(model, out_filename)


if __name__ == '__main__':
    main(sys.argv[1:])
