'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
import matplotlib.pyplot as plt


def plot(df, ylabel, filename):
    '''Plot.'''
    for column in df:
        plt.plot(df.index.values, df[column])

    plt.legend(df.columns)
    plt.xlabel(df.index.name)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(filename)
