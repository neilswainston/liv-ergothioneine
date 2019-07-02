'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
from functools import partial
import os.path

from cobra import Metabolite, Reaction

import pandas as pd


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


def save(model, solution, filename):
    '''Save solution.'''
    makedirs(filename)
    df = solution.to_frame()

    # Remove zero fluxes and sort:
    df = df[df['fluxes'].abs() > 1e-12]
    df.sort_values('fluxes', ascending=False, inplace=True)

    # Produce user-friendly output:
    df.reset_index(level=0, inplace=True)
    response = df['index'].apply(partial(_get_react_details, model))
    df[['reaction_name', 'reaction_def']] = response

    df.to_csv(filename)


def makedirs(filename):
    '''Make directories.'''
    dir_name = os.path.dirname(filename)

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)


def _get_react_details(model, row):
    '''Get reaction details.'''
    react = model.reactions.get_by_id(row)
    return pd.Series([react.name,
                      react.build_reaction_string(use_metabolite_names=True)])
