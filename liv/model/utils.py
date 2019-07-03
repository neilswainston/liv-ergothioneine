'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=protected-access
# pylint: disable=too-many-arguments
from functools import partial
import os.path

from cobra import Metabolite, Reaction

import pandas as pd


def to_df(model):
    '''Convert model to DataFrame.'''
    data = [[react.id, react.name, react.lower_bound, react.upper_bound,
             react.build_reaction_string(use_metabolite_names=True)]
            for react in model.reactions]

    return pd.DataFrame(data, columns=['id', 'name',
                                       'lower_bound', 'upper_bound',
                                       'def']).set_index('id')


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


def get_flux_df(model, solution):
    '''Get flux DataFrame.'''
    df = solution.to_frame()

    # Remove zero fluxes:
    df = df[df['fluxes'].abs() > 1e-12]

    # Produce user-friendly output:
    df.reset_index(level=0, inplace=True)
    response = df.apply(partial(_get_react_details, model), axis=1)
    df[['reaction_name', 'reaction_def', 'fluxes']] = response
    df.set_index('index', inplace=True)

    # Sort:
    df.sort_values('fluxes', ascending=False, inplace=True)

    return df


def makedirs(filename):
    '''Make directories.'''
    dir_name = os.path.dirname(filename)

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)


def _get_react_details(model, row):
    '''Get reaction details.'''
    react = model.reactions.get_by_id(row['index'])

    if row['fluxes'] < 0:
        react = _reverse_react(react)

    return pd.Series([react.name,
                      react.build_reaction_string(use_metabolite_names=True),
                      abs(row['fluxes'])])


def _reverse_react(react):
    '''Reverse reaction.'''
    rev_react = Reaction()
    rev_react.name = react.name + ' (rev)'
    rev_react.lower_bound = -react.upper_bound
    rev_react.upper_bound = -react.lower_bound

    coeffs = {met: -react.get_coefficient(met.id)
              for met in react.reactants}
    coeffs.update({met: -react.get_coefficient(met.id)
                   for met in react.products})

    rev_react.add_metabolites(coeffs)

    return rev_react
