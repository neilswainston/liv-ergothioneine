'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
import os.path
import sys

from cobra import Metabolite, Reaction
from cobra.io import read_sbml_model, write_sbml_model


def update_model(model):
    '''Get model.'''
    add_met(model, 'hercynine_c', 'hercynine', 'C9H15N3O2', 'c')
    add_met(model, 'herncystsulfox_c', 'hercynylcysteine sulfoxide',
            'C12H20N4O5S', 'c')
    add_met(model, 'sulfhercyn_c', '2-sulfenohercynine', 'C9H15N3O3S', 'c')
    add_met(model, 'ergothioneine_c', 'ergothioneine', 'C9H16N3O2S', 'c')

    # L-histidine + 3 S-adenosyl-L-methionine --> hercynine +
    # S-adenosyl-L-homocysteine
    reac_str = '''
        s_1006__91__c__93__ + 3 s_1416__91__c__93__ -->
        hercynine_c + s_1413__91__c__93__
        '''
    add_reaction(model, 'Egt1a', 'Egt1a', reac_str, 'Egt1')

    # hercynine + L-cysteine + oxygen --> hercynylcysteine sulfoxide + H20
    reac_str = '''
        hercynine_c + s_0981__91__c__93__ + s_1275__91__c__93__ -->
        herncystsulfox_c + s_0803__91__c__93__
        '''
    add_reaction(model, 'Egt1b', 'Egt1b', reac_str, 'Egt1')

    # hercynylcysteine sulfoxide + H20 -->
    # 2-sulfenohercynine + ammonium pyruvate
    reac_str = '''
        herncystsulfox_c + s_0803__91__c__93__ -->
        sulfhercyn_c + s_0419__91__c__93__ + s_4184__91__c__93__
        '''
    add_reaction(model, 'Egt2', 'Egt2', reac_str, 'Egt2')

    # 2-sulfenohercynine + electron donor (NADH) ->
    # ergothioneine + electron_acceptor (NAD) + H20
    reac_str = '''
        sulfhercyn_c + s_1203__91__c__93__ -->
        ergothioneine_c + s_1198__91__c__93__ + s_0803__91__c__93__
        '''
    add_reaction(model, 'ergothioneine_synthesis', 'ergothioneine synthesis',
                 reac_str)

    return model


def add_creator(model, family_name, given_name, organisation, email):
    '''Add creator.'''
    model._sbml['creators'] = [{'familyName': family_name,
                                'givenName': given_name,
                                'organisation': organisation,
                                'email': email}]


def add_met(model, met_id, name, formula, compartment):
    '''Add metabolite.'''
    met = Metabolite(met_id, formula=formula,
                     name=name, compartment=compartment)

    return model.add_metabolites([met])


def add_reaction(model, reaction_id, name, reac_str,
                 gene_reaction_rule=None, subsystem=None):
    '''Add reaction.'''
    reaction = Reaction(reaction_id)
    model.add_reaction(reaction)
    reaction.name = name
    reaction.build_reaction_from_string(reac_str)

    if gene_reaction_rule:
        reaction.gene_reaction_rule = gene_reaction_rule

    reaction.subsystem = subsystem
    return reaction


def main(args):
    '''main method.'''
    model_in = read_sbml_model(args[0])
    model_out = update_model(model_in)
    add_creator(model_out, *args[2:6])

    if not os.path.exists(os.path.dirname(args[1])):
        os.makedirs(os.path.dirname(args[1]))

    write_sbml_model(model_out, args[1])


if __name__ == '__main__':
    main(sys.argv[1:])
