# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 19:15:35 2025

@author: biss3
"""

import numpy as np

from sep2p import mixture_prop
from sep2p import activity_coefficient_models
from sep2p import system_parameters


def test_molecular_mass_mixture_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = system_parameters.SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = mixture_prop.molecular_mass_mixture(system.parameters, x)
    expected = (78.11*0.25 + 92.14*0.75)/1000
    assert abs(calculated - expected) < 1e-1
    
def test_density_liquid_mixture_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = system_parameters.SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = mixture_prop.density_liquid_mixture(system.parameters, T, x)
    expected = 816.0
    assert abs(calculated - expected) < 1e-1

def test_density_gas_mixture_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = system_parameters.SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = mixture_prop.density_gas_mixture(system.parameters, T, P, x)
    expected = 3.086
    assert abs(calculated - expected) < 1e-2

def test_enthalpy_gas_mixture_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = system_parameters.SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = mixture_prop.enthalpy_gas_mixture(system.parameters, T, P, x)
    expected = 61140.76
    assert abs(calculated - expected) < 1e-0

def test_enthalpy_liquid_mixture_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = system_parameters.SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = mixture_prop.enthalpy_liquid_mixture(system.parameters, T, P, x)
    expected = 26781.176
    assert abs(calculated - expected) < 1e-0

def test_vapour_liquid_equilibrium_constant_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = system_parameters.SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = mixture_prop.vapour_liquid_equilibrium_constant(system.parameters, T, P, x)[0]
    expected = np.array([0.90300, 0.34317])
    assert abs(np.sum(calculated - expected)) < 1e-2
    
def test_activity_coefficient_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = system_parameters.SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = mixture_prop.activity_coefficient(system.parameters, T, P, x)
    expected = np.array([1, 1])
    assert abs(np.sum(calculated - expected)) < 1e-0
    
def test_const_pressure_heat_capacity_liquid_mixture_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = system_parameters.SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = mixture_prop.const_pressure_heat_capacity_liquid_mixture(system.parameters, T, P, x)
    expected = 175.119
    assert abs(calculated - expected) < 1e-1


def test_activity_coefficient_unifac_paper():
    T = 307
    x1 = 0.047
    P = 101325
    components = ["Acetone", "Pentane"]
    system = system_parameters.SystemParameters(components, model_liquid="unifac1p")
    x = np.array([x1, 1 - x1])
    
    calculated = activity_coefficient_models.unifac(system.parameters, T, x)
    expected = np.array([4.9, 1.0])
    
    assert (abs(calculated - expected) < 1e-1).all()