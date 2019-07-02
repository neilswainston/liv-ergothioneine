'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
from cobra import Metabolite, Reaction


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
                 gene_reaction_rule=None, subsystem=None, check=True):
    '''Add reaction.'''
    reaction = Reaction(reaction_id)
    model.add_reaction(reaction)
    reaction.name = name
    reaction.build_reaction_from_string(reac_str)

    if gene_reaction_rule:
        reaction.gene_reaction_rule = gene_reaction_rule

    reaction.subsystem = subsystem

    balance_error = reaction.check_mass_balance()

    if check and balance_error:
        raise ValueError(balance_error)

    return reaction
