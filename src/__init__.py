# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 12:53:12 2025

@author: biss3
"""

# Re-export public API from submodules
from .util_general import gas_constant, acceleration_constant
from .pure_comp_prop import density_liquid, density_gas, surface_tension, saturation_pressure, enthalpy_gas, enthalpy_liquid, heat_of_vaporisation, const_pressure_heat_capacity_gas, const_pressure_heat_capacity_liquid
from .mixture_prop import molecular_mass_mixture, density_liquid_mixture, density_gas_mixture, enthalpy_gas_mixture, enthalpy_liquid_mixture, vapour_liquid_equilibrium_constant, activity_coefficient, const_pressure_heat_capacity_liquid_mixture
from .data_loader import load_pure_comp_properties, load_liquid_mixture_model
from .system_parameters import SystemParameters
from .activity_coefficient_models import unifac
from .equilibrium_calculations import solve_flash, binary_phase_diagram_Ty
#from .models import User

# Optional metadata
__version__ = "0.0.1"
__author__ = "Thomas Bisgaard Astrup"

# Define what gets imported with: from my_project import *
__all__ = ["add", "User"]

