# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 13:56:22 2025

@author: biss3
"""

from data_loader import load_pure_comp_properties, load_liquid_mixture_model

class SystemParameters:
    def __init__(self, components, model_liquid="ideal"):
        if isinstance(components, str):
            # Convert to string for iterable
            components = [components]
        self.components = components  # list of Component objects
        parameters_pure_comp = load_pure_comp_properties(components)
        parameters_mixture_liquid = load_liquid_mixture_model(components, model_liquid)
        self.parameters = parameters_pure_comp | parameters_mixture_liquid
        self.parameters["number_of_components"] = len(components)
        self.column = {
            'number_of_stages':         20,  #[-]
            'stage_feed':               [11],  #[-]
            'constant_liquid':          1/0.2**1.5,  #[1/s]
            'constant_vapour':          1/(700/101325)**0.5,  #[1/s]
            'weir_height':              0.02,  #[m]
            'weir_crossarea':           (0.2)**2*3.14/4,  #[m^2]
            'stage_height':             0.4  #[m]
            }