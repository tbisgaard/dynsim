# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 15:29:17 2025

@author: biss3
"""

import sys
import os
import pytest


from .pure_comp_prop import density_liquid, density_gas, surface_tension, saturation_pressure, enthalpy_gas, enthalpy_liquid, heat_of_vaporisation, const_pressure_heat_capacity_gas, const_pressure_heat_capacity_liquid
from .system_parameters import SystemParameters


def test_density_liquid_known_component():
    """
    Test density_liquid with a component that exists in the database.
    """
    # Example compound in CSV
    components = "Water"
    system = SystemParameters(components)
    T = 25 + 273.15

    calculated = density_liquid(system.parameters, T)

    # Example expected value (depends on CSV coefficients)
    expected = 995.0# 5.459/(0.30542**(1+(1-(25+273.15)/647.13)**0.081))*1000*0.018015
    
    assert abs(calculated - expected) < 1e-1


def test_density_gas_known_conditions():
    """
    Test density_vapour at normal conditions
    """
    components = "Water"
    system = SystemParameters(components)
    T = 25 + 273.15
    P = 101325
    
    calculated = density_gas(system.parameters, T, P)
    
    expected = 101325/8.314/(25+273.15)*0.018015
    
    assert abs(calculated - expected) < 1e-2


def test_saturation_pressure_known_component():
    components = "Water"
    system = SystemParameters(components)
    T = 100 + 273.15
    
    calculated = saturation_pressure(system.parameters, T)
    
    expected = 101325
    
    assert abs(calculated - expected) < 1e1


def test_surface_tension_known_component():
    components = "Water"
    system = SystemParameters(components)
    T = 25 + 273.15
    
    calculated = surface_tension(system.parameters, T)
    
    expected = 0.0727
    
    assert abs(calculated - expected) < 1e-2


def test_enthalpy_gas_component():
    T = 350
    P = 101325
    components = "Water"
    system = SystemParameters(components)
    
    calculated = enthalpy_gas(system.parameters, T, P)
    
    expected = -241792.04  # [J/mol]
    assert abs(calculated - expected) < 1e-0
    
def test_heat_of_vaporisation_component():
    T = 350
    components = "Water"
    system = SystemParameters(components)
    
    calculated = heat_of_vaporisation(system.parameters, T)
    
    expected = 41836.777  # [J/mol]
    assert abs(calculated - expected) < 1e-1


def test_enthalpy_liquid():
    T = 350
    P = 101325
    components = "Water"
    system = SystemParameters(components)
    
    calculated = enthalpy_liquid(system.parameters, T, P)
    
    expected = -283628.82  # [J/mol]
    assert abs(calculated - expected) < 1e-0
    
def test_const_pressure_heat_capacity_gas():
    T = 350
    P = 101325
    components = "Water"
    system = SystemParameters(components)
    
    calculated = const_pressure_heat_capacity_gas(system.parameters, T, P)
    
    expected = 33.86  # [J/mol/K]
    assert abs(calculated - expected) < 1e-2
    
def test_const_pressure_heat_capacity_liquid():
    T = 350
    P = 101325
    components = "Water"
    system = SystemParameters(components)
    
    calculated = const_pressure_heat_capacity_liquid(system.parameters, T, P)
    
    expected = 76.76  # [J/mol/K]
    assert abs(calculated - expected) < 1e-2