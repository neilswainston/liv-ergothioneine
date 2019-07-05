'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=wrong-import-order
import os.path
import sys

from cobra.io import read_sbml_model, write_sbml_model

from liv.ergothioneine.build import build
from liv.model.plot import plot
from liv.model.utils import add_creator, get_flux_df, get_mw, makedirs, to_df
import numpy as np
import pandas as pd


def simulate(model, out_dir):
    '''Simulate model.'''
    # Biomass:
    biomass_max = _simulate(model).objective_value
    biomass_react = model.reactions.get_by_id('r_2111')

    # Biomass, ergothioneine, glucose, L-histidine, L-cysteine:
    react_ids = ['r_2111', 'ergothioneine_sink', 'r_1714', 'r_1201', 'r_1192']

    react_names = [model.reactions.get_by_id(react_id).name
                   for react_id in react_ids]

    react_fluxes = []

    for prop in np.arange(0.0, 1.1, 0.1):
        biomass_flux = biomass_max * prop
        biomass_react.lower_bound = biomass_flux
        biomass_react.upper_bound = biomass_flux

        model.objective = 'ergothioneine_sink'
        solution = _simulate(model)

        react_fluxes.append([prop] + solution.fluxes[react_ids].tolist())

        # Format and save flux solution:
        filename = os.path.join(out_dir, '%.1f_biomass_fluxes.csv' % prop)
        flux_df = get_flux_df(model, solution)
        makedirs(filename)
        flux_df.to_csv(filename)

    # Save reaction fluxes:
    index_name = 'Max biomass proportion'
    react_flux_df = pd.DataFrame(react_fluxes,
                                 columns=[index_name] + react_names)

    # react_flux_df = react_flux_df.abs()
    react_flux_df.where(react_flux_df > 1e-12, 0, inplace=True)
    react_flux_df.set_index(index_name, inplace=True)
    react_flux_df.name = 'react_fluxes'

    return react_flux_df


def _simulate(model):
    '''Simulate model.'''
    solution = model.optimize()
    # print(model.summary(names=True))
    return solution


def main(args):
    '''main method.'''
    model = read_sbml_model(args[0])

    # Get biomass MW:
    print(get_mw(model, 's_0450__91__c__93__'))

    # Create updated model:
    build(model)
    add_creator(model, *args[3:7])

    # Write updated model:
    makedirs(args[1])
    write_sbml_model(model, os.path.join(args[1], '%s.xml' % args[2]))
    to_df(model).to_csv(os.path.join(args[1], '%s.csv' % args[2]))

    # Simulate updated model:
    react_flux_df = simulate(model, args[1])

    # Save and plot:
    react_flux_df.to_csv(os.path.join(args[1], '%s.csv' % react_flux_df.name))

    plot(react_flux_df, 'flux / mmol h-1',
         os.path.join(args[1], '%s.png' % react_flux_df.name))


if __name__ == '__main__':
    main(sys.argv[1:])
