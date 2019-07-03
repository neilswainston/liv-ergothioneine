'''
(c) University of Liverpool 2019

All rights reserved.

@author: neilswainston
'''
from liv.model.utils import add_met, add_reaction


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
        s_1006__91__c__93__ + 3 s_1416__91__c__93__ -->
        hercynine_c + 3 s_1413__91__c__93__ + 3 s_0794__91__c__93__
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

    # 2-sulfenohercynine + electron donor (glutathione) ->
    # ergothioneine + electron_acceptor (glutathione disulfide) + H20
    reac_str = '''
        sulfhercyn_c + 2 s_0750__91__c__93__ -->
        ergothioneine_c + s_0754__91__c__93__ + s_0803__91__c__93__
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

    # ATP maintenance:
    # r_4046

    # L-histidine transport:
    model.reactions.get_by_id('r_1893').lower_bound = -1

    # L-cysteine transport:
    model.reactions.get_by_id('r_1883').lower_bound = -1
