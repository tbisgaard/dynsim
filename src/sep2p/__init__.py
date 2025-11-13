# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 12:53:12 2025

@author: biss3
"""

# Re-export public API from submodules
from sep2p import utils
from sep2p import pure_comp_prop
from sep2p import mixture_prop
from sep2p import data_loader
from sep2p import system_parameters
from sep2p import activity_coefficient_models
from sep2p import equilibrium_calculations
from sep2p import column_models
from sep2p import phys_constants
#from .models import User

# Optional metadata
__version__ = "0.0.1"
__author__ = "Thomas Bisgaard Astrup"

# Define what gets imported with: from my_project import *
__all__ = ["add", "User"]

