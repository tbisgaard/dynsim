# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 12:58:07 2025

@author: biss3
"""


import os
import numpy as np
from util_general import gas_constant

def molecular_mass(parameters):
    """
    

    Parameters
    ----------
    components : String or list of strings.
        Components required for property calculation of length (NC).

    Returns
    -------
    molecular_mass_components : Numpy array.
        Moleculat mass in unit [kg/mol] of length (NC).

    """
    MW = parameters["molecular_mass"]#get_component_values(molecular_mass_cache, components)#np.array([molecular_mass_cache[c] for c in components])#pd.to_numeric(db.loc[components]["molecular_mass"])/1000
    
    MW = MW/1000  # [kg/mol]
    
    molecular_mass_components = np.array(MW).astype(float)
    if molecular_mass_components.shape == (1, 1):
        return molecular_mass_components.item()
    return  np.array(molecular_mass_components)



def density_liquid(parameters, temperature):
    """
    

    Parameters
    ----------
    temperature : Numpy array of length (NS).
        Temperature in [K].
    par_physical : TYPE
        Pure component database.

    Returns
    -------
    density_liquid : Numpy array of length (NS).
        Density of liquid according to DIPPR 105, in [kg/m3].

    """
    
    MW = molecular_mass(parameters)
    A = parameters["molar_density_constant_A"]#get_component_values(molar_density_constant_A_cache, components)#pd.to_numeric(db.loc[components, "molar_density_constant_A"])
    B = parameters["molar_density_constant_B"]#get_component_values(molar_density_constant_B_cache, components)#pd.to_numeric(db.loc[components, "molar_density_constant_B"])
    C = parameters["molar_density_constant_C"]#get_component_values(molar_density_constant_C_cache, components)#pd.to_numeric(db.loc[components, "molar_density_constant_C"])
    D = parameters["molar_density_constant_D"]#get_component_values(molar_density_constant_D_cache, components)#pd.to_numeric(db.loc[components, "molar_density_constant_D"])
    
    MW = np.array(MW, ndmin=1)[None, :]
    A = np.array(A, ndmin=1)[None, :]
    B = np.array(B, ndmin=1)[None, :]
    C = np.array(C, ndmin=1)[None, :]
    D = np.array(D, ndmin=1)[None, :]
    
    T = np.array(temperature, ndmin=1)[:, None]
    
    rho = A/(B**(1+(1-T/C)**D))*MW
    
    rho = rho*1000  # [kg/m3]
    
    density_liquid_components = rho
    
    if density_liquid_components.shape == (1, 1):
        return density_liquid_components.item()
    return density_liquid_components


def density_gas(parameters, temperature, pressure):
    """
    
    
    Parameters
    ----------
    temperature : Numpy array of length (NS).
        Temperature in [K].
    par_physical : Dictionary.
        Pure component database.

    Returns
    -------
    density_gas : Numpy array of length (NS).
        Density of gas in [kg/m3].
    """
    acentric_factor = 1
    MW = molecular_mass(parameters)
    
    MW = np.array(MW, ndmin=1)[None, :]
    T = np.array(temperature, ndmin=1)[:, None]
    P = np.array(pressure, ndmin=1)[:, None]
    
    density_gas_components = MW*P/T/gas_constant()/acentric_factor
    
    if density_gas_components.shape == (1, 1):
        return density_gas_components.item()
    return density_gas_components


def surface_tension(parameters, temperature):
    """
    
    
    Parameters
    ----------
    temperature : Numpy array of length (NS).
        Temperature in [K].
    par_physical : Dictionary.
        Pure component database.

    Returns
    -------
    surface_tension : Numpy array of length (NS).
        Surface tension of mixture by Full DIPPR 106 in [N/m]

    """

    A = parameters["surface_tension_constant_A"]#get_component_values(surface_tension_constant_A_cache, components)#pd.to_numeric(db.loc[components]["surface_tension_constant_A"])
    B = parameters["surface_tension_constant_B"]#get_component_values(surface_tension_constant_B_cache, components)#pd.to_numeric(db.loc[components]["surface_tension_constant_B"])
    C = parameters["surface_tension_constant_C"]#get_component_values(surface_tension_constant_C_cache, components)#pd.to_numeric(db.loc[components]["surface_tension_constant_C"])
    D = parameters["surface_tension_constant_D"]#get_component_values(surface_tension_constant_D_cache, components)#pd.to_numeric(db.loc[components]["surface_tension_constant_D"])
    Tc = parameters["temperature_critical"]#get_component_values(temperature_critical_cache, components)#pd.to_numeric(db.loc[components]["temperature_critical"])
    
    A = np.array(A, ndmin=1)[None, :]
    B = np.array(B, ndmin=1)[None, :]
    C = np.array(C, ndmin=1)[None, :]
    D = np.array(D, ndmin=1)[None, :]
    Tc = np.array(Tc, ndmin=1)[None, :]
    
    T = np.array(temperature, ndmin=1)[:, None]
    
    Tr = T/Tc
    
    exponent = B + C*Tr + D*Tr**2
    surface_tension_components = A*(1 - Tr)**exponent
    
    if surface_tension_components.shape == (1, 1):
        return surface_tension_components.item()
    return surface_tension_components


def saturation_pressure(parameters, temperature):
    """
    

    Parameters
    ----------
    temperature : Numpy array of length (NS).
        Temperature in [K].
    par_physical : TYPE
        Pure component database.

    Returns
    -------
    saturation_pressure : Numpy array of length (NS).
        Saturation pressure based on DIPPR 101 in [Pa].

    """
    # Vapor pressure based on DIPPR 101
    A = parameters["saturation_pressure_constant_A"]#get_component_values(saturation_pressure_constant_A_cache, components)#pd.to_numeric(db.loc[components]["saturation_pressure_constant_A"])
    B = parameters["saturation_pressure_constant_B"]#get_component_values(saturation_pressure_constant_B_cache, components)#pd.to_numeric(db.loc[components]["saturation_pressure_constant_B"])
    C = parameters["saturation_pressure_constant_C"]#get_component_values(saturation_pressure_constant_C_cache, components)#pd.to_numeric(db.loc[components]["saturation_pressure_constant_C"])
    D = parameters["saturation_pressure_constant_D"]#get_component_values(saturation_pressure_constant_D_cache, components)#pd.to_numeric(db.loc[components]["saturation_pressure_constant_D"])
    E = parameters["saturation_pressure_constant_E"]#get_component_values(saturation_pressure_constant_E_cache, components)#pd.to_numeric(db.loc[components]["saturation_pressure_constant_E"])
    
    A = np.array(A, ndmin=1)[None, :]
    B = np.array(B, ndmin=1)[None, :]
    C = np.array(C, ndmin=1)[None, :]
    D = np.array(D, ndmin=1)[None, :]
    E = np.array(E, ndmin=1)[None, :]
    
    T = np.array(temperature, ndmin=1)[:, None]


    log_saturation_pressure_components = (A + B/T + C * np.log(T) + D * T**E)
    saturation_pressure_components = np.exp(log_saturation_pressure_components)  
    
    if saturation_pressure_components.shape == (1, 1):
        return saturation_pressure_components.item()
    return saturation_pressure_components



def enthalpy_gas(parameters, temperature, pressure):
    h_f = parameters["enthalpy_of_formation_ideal_gas"]#get_component_values(enthalpy_of_formation_ideal_gas_cache, components)#pd.to_numeric(db.loc[components]["enthalpy_of_formation_ideal_gas"])
    A = parameters["heat_capacity_ideal_gas_constant_A"]#get_component_values(heat_capacity_ideal_gas_constant_A_cache, components)#pd.to_numeric(db.loc[components]["heat_capacity_ideal_gas_constant_A"])
    B = parameters["heat_capacity_ideal_gas_constant_B"]#get_component_values(heat_capacity_ideal_gas_constant_B_cache, components)#pd.to_numeric(db.loc[components]["heat_capacity_ideal_gas_constant_B"])
    C = parameters["heat_capacity_ideal_gas_constant_C"]#get_component_values(heat_capacity_ideal_gas_constant_C_cache, components)#pd.to_numeric(db.loc[components]["heat_capacity_ideal_gas_constant_C"])
    D = parameters["heat_capacity_ideal_gas_constant_D"]#get_component_values(heat_capacity_ideal_gas_constant_D_cache, components)#pd.to_numeric(db.loc[components]["heat_capacity_ideal_gas_constant_D"])
    E = parameters["heat_capacity_ideal_gas_constant_E"]#get_component_values(heat_capacity_ideal_gas_constant_E_cache, components)#pd.to_numeric(db.loc[components]["heat_capacity_ideal_gas_constant_E"])
    
    A = np.array(A, ndmin=1)[None, :]
    B = np.array(B, ndmin=1)[None, :]
    C = np.array(C, ndmin=1)[None, :]
    D = np.array(D, ndmin=1)[None, :]
    E = np.array(E, ndmin=1)[None, :]
    h_f = np.array(h_f, ndmin=1)[None, :]

    T = np.array(temperature, ndmin=1)[:, None]
    Tref = 298.15*np.ones_like(T)
    
    # Integration of Cp
    h_T_ub = T*A + B*C*np.cosh(C/T)/np.sinh(C/T) \
        - D*E*np.sinh(E/T)/np.cosh(E/T)
    h_T_lb = T*A + B*C*np.cosh(C/Tref)/np.sinh(C/Tref) \
        - D*E*np.sinh(E/Tref)/np.cosh(E/Tref)
    h_T = h_T_ub - h_T_lb
    h_T = h_T/1000  # [J/mol]
    h_P = 0
    h_E = 0
    h = h_f + h_T + h_P + h_E
    
    
    enthalpy_gas_components = h
    if enthalpy_gas_components.shape == (1, 1):
        return enthalpy_gas_components.item()
    return enthalpy_gas_components

def heat_of_vaporisation(parameters, temperature):
    A = parameters["heat_of_vaporisation_constant_A"]#get_component_values(heat_of_vaporisation_constant_A_cache, components)#pd.to_numeric(db.loc[components]["heat_of_vaporisation_constant_A"])
    B = parameters["heat_of_vaporisation_constant_B"]#get_component_values(heat_of_vaporisation_constant_B_cache, components)#pd.to_numeric(db.loc[components]["heat_of_vaporisation_constant_B"])
    C = parameters["heat_of_vaporisation_constant_C"]#get_component_values(heat_of_vaporisation_constant_C_cache, components)#pd.to_numeric(db.loc[components]["heat_of_vaporisation_constant_C"])
    D = parameters["heat_of_vaporisation_constant_D"]#get_component_values(heat_of_vaporisation_constant_D_cache, components)#pd.to_numeric(db.loc[components]["heat_of_vaporisation_constant_D"])
    Tc =  parameters["temperature_critical"]#get_component_values(temperature_critical_cache, components)#pd.to_numeric(db.loc[components]["temperature_critical"])
    
    A = np.array(A, ndmin=1)[None, :]
    B = np.array(B, ndmin=1)[None, :]
    C = np.array(C, ndmin=1)[None, :]
    D = np.array(D, ndmin=1)[None, :]
    Tc = np.array(Tc, ndmin=1)[None, :]
    
    T = np.array(temperature, ndmin=1)[:, None]
    
    Tr = T/Tc
    
    exponent = B + C*Tr  +D*Tr**2
    dhvap = A*(1-Tr)**(exponent)*(Tr<1)
    dhvap = dhvap/1000  # [J/mol]

    heat_of_vaporisation_components = dhvap
    if heat_of_vaporisation_components.shape == (1, 1):
        return heat_of_vaporisation_components.item()
    return heat_of_vaporisation_components


def enthalpy_liquid(parameters, temperature, pressure):
    
    h_G = enthalpy_gas(parameters, temperature, pressure)
    dhvap = heat_of_vaporisation(parameters, temperature)
    h_L = h_G - dhvap
    
    enthalpy_liquid_components = h_L
    
    if np.isscalar(enthalpy_liquid_components):
        return enthalpy_liquid_components    
    if enthalpy_liquid_components.shape == (1, 1):
        return enthalpy_liquid_components.item()
    return enthalpy_liquid_components


def const_pressure_heat_capacity_gas(parameters, temperature, pressure):
    A = parameters["heat_capacity_ideal_gas_constant_A"]#get_component_values(heat_capacity_ideal_gas_constant_A_cache, components)#pd.to_numeric(db.loc[components]["heat_capacity_ideal_gas_constant_A"])
    B = parameters["heat_capacity_ideal_gas_constant_B"]#get_component_values(heat_capacity_ideal_gas_constant_B_cache, components)#pd.to_numeric(db.loc[components]["heat_capacity_ideal_gas_constant_B"])
    C = parameters["heat_capacity_ideal_gas_constant_C"]#get_component_values(heat_capacity_ideal_gas_constant_C_cache, components)#pd.to_numeric(db.loc[components]["heat_capacity_ideal_gas_constant_C"])
    D = parameters["heat_capacity_ideal_gas_constant_D"]#get_component_values(heat_capacity_ideal_gas_constant_D_cache, components)#pd.to_numeric(db.loc[components]["heat_capacity_ideal_gas_constant_D"])
    E = parameters["heat_capacity_ideal_gas_constant_E"]#get_component_values(heat_capacity_ideal_gas_constant_E_cache, components)#pd.to_numeric(db.loc[components]["heat_capacity_ideal_gas_constant_E"])

    A = np.array(A, ndmin=1)[None, :]
    B = np.array(B, ndmin=1)[None, :]
    C = np.array(C, ndmin=1)[None, :]
    D = np.array(D, ndmin=1)[None, :]
    E = np.array(E, ndmin=1)[None, :]

    T = np.array(temperature, ndmin=1)[:, None]  
    
    CP = A+B*(C/T/np.sinh(C/T))**2 +D*(E/T/np.cosh(E/T))**2
    
    CP = CP/1000  # [J/mol/K]
    
    const_pressure_heat_capacity_gas_components = CP
    if const_pressure_heat_capacity_gas_components.shape == (1, 1):
        return const_pressure_heat_capacity_gas_components.item()
    return const_pressure_heat_capacity_gas_components

def const_pressure_heat_capacity_liquid(parameters, temperature, pressure):
    
    CP_G = const_pressure_heat_capacity_gas(parameters, temperature, pressure)
    
    A = parameters["heat_of_vaporisation_constant_A"]#get_component_values(heat_of_vaporisation_constant_A_cache, components)#pd.to_numeric(db.loc[components]["heat_of_vaporisation_constant_A"])
    B = parameters["heat_of_vaporisation_constant_B"]#get_component_values(heat_of_vaporisation_constant_B_cache, components)#pd.to_numeric(db.loc[components]["heat_of_vaporisation_constant_B"])
    C = parameters["heat_of_vaporisation_constant_C"]#get_component_values(heat_of_vaporisation_constant_C_cache, components)#pd.to_numeric(db.loc[components]["heat_of_vaporisation_constant_C"])
    D = parameters["heat_of_vaporisation_constant_D"]#get_component_values(heat_of_vaporisation_constant_D_cache, components)#pd.to_numeric(db.loc[components]["heat_of_vaporisation_constant_D"])
    Tc = parameters["temperature_critical"]#get_component_values(temperature_critical_cache, components)#pd.to_numeric(db.loc[components]["temperature_critical"])
    
    A = np.array(A, ndmin=1)[None, :]
    B = np.array(B, ndmin=1)[None, :]
    C = np.array(C, ndmin=1)[None, :]
    D = np.array(D, ndmin=1)[None, :]
    Tc = np.array(Tc, ndmin=1)[None, :]
    
    T = np.array(temperature, ndmin=1)[:, None]
    
    Tr = T/Tc
    
    exponent = B + C*Tr + D*Tr**2
    dCP = A*(1-Tr)**(exponent)*(
        (C/Tc + 2*D*T/Tc**2)*np.log(1-Tr) - (B + C*Tr + D*Tr**2)/(Tc*(1-Tr))
        )
    CP_L = CP_G - dCP/1000
    
    const_pressure_heat_capacity_liquid_components = CP_L
    if const_pressure_heat_capacity_liquid_components.shape == (1, 1):
        return const_pressure_heat_capacity_liquid_components.item()
    return const_pressure_heat_capacity_liquid_components


