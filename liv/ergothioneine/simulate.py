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

    reacts = ['r_2111', 'ergothioneine_sink']

    for prop in np.arange(0, 1.1, 0.1):
        model.objective = {model.reactions.get_by_id(react): coeff
                           for react, coeff in zip(reacts, [prop, 1 - prop])}
        solution = _simulate(model)
        print('%.1f' % prop, solution.fluxes[reacts].values)


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
