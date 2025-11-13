# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 19:02:17 2025

@author: biss3
"""

import numpy as np
from sep2p import utils, pure_comp_prop, activity_coefficient_models

def molecular_mass_mixture(parameters, composition):
    
    MW = pure_comp_prop.molecular_mass(parameters)
    
    MW = np.array(MW, ndmin=1)[None, :]
    x = np.atleast_2d(composition)
    
    molecular_mass_mixture_components = np.sum(x*MW, axis=1)
    
    if molecular_mass_mixture_components.size == 1:
        return molecular_mass_mixture_components.item()
    return molecular_mass_mixture_components


def density_liquid_mixture(parameters, temperature, composition):
    rho0_L = pure_comp_prop.density_liquid(parameters, temperature)
    MW0 = pure_comp_prop.molecular_mass(parameters)
    MW_avg = molecular_mass_mixture(parameters, composition)
    
    x = np.atleast_2d(composition)

    rhom0_L = rho0_L/MW0
    density_liquid_mixture_components = MW_avg * (1/np.sum(x/rhom0_L, axis=1))

    if density_liquid_mixture_components.size == 1:
        return density_liquid_mixture_components.item()
    return density_liquid_mixture_components



def density_gas_mixture(parameters, temperature, pressure, composition):
    rho0_G = pure_comp_prop.density_gas(parameters, temperature, pressure)
    MW0 = pure_comp_prop.molecular_mass(parameters)
    MW_avg = molecular_mass_mixture(parameters, composition)
    
    x = np.atleast_2d(composition)

    rhom0_G = rho0_G/MW0
    density_gas_mixture_components = MW_avg * (1/np.sum(x/rhom0_G, axis=1))
    
    if density_gas_mixture_components.size == 1:
        return density_gas_mixture_components.item()
    return density_gas_mixture_components



def enthalpy_gas_mixture(parameters, temperature, pressure, composition):
    h_G = pure_comp_prop.enthalpy_gas(parameters, temperature, pressure)
    
    x = np.atleast_2d(composition)
    
    enthalpy_gas_mixture_components = np.sum(x*h_G, axis=1)
    
    if enthalpy_gas_mixture_components.size == 1:
        return enthalpy_gas_mixture_components.item()
    return enthalpy_gas_mixture_components
    

def enthalpy_liquid_mixture(parameters, temperature, pressure, composition):
    h_L = pure_comp_prop.enthalpy_liquid(parameters, temperature, pressure)
    
    x = np.atleast_2d(composition)
    
    enthalpy_liquid_mixture_components = np.sum(x*h_L, axis=1)
    
    if enthalpy_liquid_mixture_components.size == 1:
        return enthalpy_liquid_mixture_components.item()
    return enthalpy_liquid_mixture_components



def vapour_liquid_equilibrium_constant(parameters, temperature, pressure, composition):
    Psat = pure_comp_prop.saturation_pressure(parameters, temperature)
    gamma = activity_coefficient(parameters, temperature, pressure, composition)
    
    P = np.array(pressure, ndmin=1)[:, None]
    x = np.atleast_2d(composition)
    
    K = gamma * Psat / P
    y = K * x
    return K, y


def activity_coefficient(parameters, temperature, pressure, composition):
    if (parameters["model_liquid"].lower()=="unifac1p"):
        gamma = activity_coefficient_models.unifac(parameters, temperature, composition)
    else:
        gamma = np.ones_like(composition)  # Ideal
    return gamma


def const_pressure_heat_capacity_liquid_mixture(parameters, temperature, pressure, composition):
    CP_L = pure_comp_prop.const_pressure_heat_capacity_liquid(parameters, temperature, pressure)
    
    x = np.atleast_2d(composition)
    
    const_pressure_heat_capacity_liquid_mixture_components = np.sum(x*CP_L, axis=1)
    
    if const_pressure_heat_capacity_liquid_mixture_components.size == 1:
        return const_pressure_heat_capacity_liquid_mixture_components.item()
    return const_pressure_heat_capacity_liquid_mixture_components


#     case 's'
#         %% Entropy
#         % Refence pressure
#         Pref = f_vapor_pressure(PhysicalData,Tref*ones(NS,1));
        
#         % Vapor and liquid excess enthalpy
#         sL_E = -R*sum(x.*log(x),2);% Ideal entropy of mixing  [J/kmol/K]
#         sV_E = -R*sum(x.*log(x),2);
#         if any(PhysicalData.VaporLiquidEquilibrium.modelType==1)
#             [~,~,~,~,sL_EPhi,sV_EPhi] = f_fugacity_coefficient_SRK(PhysicalData,T,P,x);
#             sL_E = sL_E+sL_EPhi;
#             sV_E = sV_E+sV_EPhi;
#         elseif PhysicalData.VaporLiquidEquilibrium.modelType(1)==2
#             gammaL = f_activity_coefficient(PhysicalData,T,x);
#             sL_EgammaL = -R*sum(x.*log(gammaL),2);
#             sL_E = sL_E+sL_EgammaL;
#         end
#         sL_E(isnan(sL_E)) = 0;
#         sV_E(isnan(sV_E)) = 0;


        
#         % Vaporization
#         dhvap = sum(x.*((A_dHvap).*(1-Tr) ...% Enthalpy of vaporization[kJ/kmol]
#             .^(B_dHvap+C_dHvap.*Tr+D_dHvap.*Tr.^2)).*(Tr<1),2);
#         dsvap = dhvap./T;                   % Entropy of vaporization       [J/kmol/K]
        
#         % Vapor
#         s_T_ub = A_igCP.*log(T_)+2*C_igCP.*B_igCP./T_ ...
#             +2*C_igCP.*B_igCP./(T_.*(exp(2*C_igCP./T_)-1)) ...
#             -B_igCP.*log(exp(2*C_igCP./T_)-1) ...
#             -2*E_igCP.*D_igCP./T_+2*E_igCP.*D_igCP ...
#             ./(T_.*(exp(2*E_igCP./T_)+1)) ...
#             +D_igCP.*log(exp(2*E_igCP./T_)+1);% Upper int. bound.  [J/kmol/K]
#         s_T_lb = A_igCP.*log(Tref_)+2*C_igCP.*B_igCP./Tref_ ...
#             +2*C_igCP.*B_igCP./(Tref_.*(exp(2*C_igCP./Tref_)-1)) ...
#             -B_igCP.*log(exp(2*C_igCP./Tref_)-1) ...
#             -2*E_igCP.*D_igCP./Tref_+2*E_igCP.*D_igCP ...
#             ./(Tref_.*(exp(2*E_igCP./Tref_)+1)) ...
#             +D_igCP.*log(exp(2*E_igCP./Tref_)+1); % Low.int.bound.[J/kmol/K]
#         sV_f = sum(x.*igSf,2);          % Entropy of formation          [J/kmol/K]
#         sV_T = sum(x.*(s_T_ub-s_T_lb),2);% Temperature dependency       [J/kmol/K]
#         sV_P = -sum(R*log(P_./Pref),2); % Pressure dependency           [J/kmol/K]
#         sV = sV_f+sV_T+sV_P+sV_E;       % Total vapor entropy           [J/kmol/K]
        
#         % Liquid
#         sL = sV-sV_E+sL_E-dsvap;        % Liquid entropy                [J/kmol/K]
        
#         % Output [kJ/mol]
#         PL = real(sL)*1e-6;
#         PV = real(sV)*1e-6;
#         dPvap = dsvap*1e-6;


#     case 'c'
#         %% Constant pressure heat capacity
#         % Vaporization
#         dhvapdT = sum(x.*(A_dHvap.*(1-Tr) ...
#             .^(B_dHvap+C_dHvap.*Tr+D_dHvap.*Tr.^2) ...
#             .*((C_dHvap./Tc+2*D_dHvap.*T_./Tc.^2).*log(1-Tr) ...
#             -(B_dHvap+C_dHvap.*Tr+D_dHvap.*Tr.^2)./(Tc.*(1-Tr)))),2);%[J/kmol/K]
        
#         % Vapor
#         CPV = sum(x.*(A_igCP+B_igCP.*(C_igCP./T_./sinh(C_igCP./T_)).^2 ...
#             +D_igCP.*(E_igCP./T_./cosh(E_igCP./T_)).^2),2); %      [J/kmol/K]
        
#         % Liquid
#         CPL = CPV-dhvapdT;
        
#         % Output [kJ/mol]
#         PL = CPL*1e-6;
#         PV = CPV*1e-6;
#         dPvap = dhvapdT*1e-6;
# end
