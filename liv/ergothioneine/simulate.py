'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
import os.path
import sys

from cobra.io import read_sbml_model, write_sbml_model

from liv.ergothioneine.build import build
from liv.ergothioneine.utils import add_creator
import numpy as np


def simulate(model):
    '''Simulate model.'''
    # Biomass:
    biomass_max = _simulate(model).objective_value
    biomass_react = model.reactions.get_by_id('r_2111')

    reacts = ['r_2111', 'ergothioneine_sink']

    for prop in np.arange(1.0, -0.1, -0.1):
        biomass_flux = biomass_max * prop
        biomass_react.lower_bound = biomass_flux
        biomass_react.upper_bound = biomass_flux

        model.objective = 'ergothioneine_sink'
        solution = _simulate(model)
        print('%.1f' % prop, ['%.3f' % val
                              for val in solution.fluxes[reacts].values])


def _simulate(model):
    '''Simulate model.'''
    solution = model.optimize()
    # print(model.summary(names=True))
    return solution


def main(args):
    '''main method.'''
    model_in = read_sbml_model(args[0])
    model_out = build(model_in)
    add_creator(model_out, *args[2:6])

    if not os.path.exists(os.path.dirname(args[1])):
        os.makedirs(os.path.dirname(args[1]))

    write_sbml_model(model_out, args[1])

    simulate(model_out)


if __name__ == '__main__':
    main(sys.argv[1:])
