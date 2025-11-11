# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 19:15:35 2025

@author: biss3
"""

import sys
import os
import pytest
import numpy as np

from .mixture_prop import molecular_mass_mixture, density_liquid_mixture, density_gas_mixture, enthalpy_gas_mixture, enthalpy_liquid_mixture, vapour_liquid_equilibrium_constant, activity_coefficient, const_pressure_heat_capacity_liquid_mixture
from .activity_coefficient_models import unifac
from .system_parameters import SystemParameters


def test_molecular_mass_mixture_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = molecular_mass_mixture(system.parameters, x)
    expected = (78.11*0.25 + 92.14*0.75)/1000
    assert abs(calculated - expected) < 1e-1
    
def test_density_liquid_mixture_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = density_liquid_mixture(system.parameters, T, x)
    expected = 816.0
    assert abs(calculated - expected) < 1e-1

def test_density_gas_mixture_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = density_gas_mixture(system.parameters, T, P, x)
    expected = 3.086
    assert abs(calculated - expected) < 1e-2

def test_enthalpy_gas_mixture_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = enthalpy_gas_mixture(system.parameters, T, P, x)
    expected = 61140.76
    assert abs(calculated - expected) < 1e-0

def test_enthalpy_liquid_mixture_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = enthalpy_liquid_mixture(system.parameters, T, P, x)
    expected = 26781.176
    assert abs(calculated - expected) < 1e-0

def test_vapour_liquid_equilibrium_constant_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = vapour_liquid_equilibrium_constant(system.parameters, T, P, x)[0]
    expected = np.array([0.90300, 0.34317])
    assert abs(np.sum(calculated - expected)) < 1e-2
    
def test_activity_coefficient_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = activity_coefficient(system.parameters, T, P, x)
    expected = np.array([1, 1])
    assert abs(np.sum(calculated - expected)) < 1e-0
    
def test_const_pressure_heat_capacity_liquid_mixture_component():
    T = 350
    P = 101325
    components = ["Benzene", "Toluene"]
    system = SystemParameters(components)
    x = np.array([0.25, 0.75])
    calculated = const_pressure_heat_capacity_liquid_mixture(system.parameters, T, P, x)
    expected = 175.119
    assert abs(calculated - expected) < 1e-1


def test_activity_coefficient_unifac_paper():
    T = 307
    x1 = 0.047
    P = 101325
    components = ["Acetone", "Pentane"]
    system = SystemParameters(components, model_liquid="unifac1p")
    x = np.array([x1, 1 - x1])
    
    calculated = unifac(system.parameters, T, x)
    expected = np.array([4.9, 1.0])
    
    assert (abs(calculated - expected) < 1e-1).all()