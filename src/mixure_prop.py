# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 13:17:44 2025

@author: biss3
"""

def density_liquid_mixture():
    NS = T.shape[0]
    NC = par_physical["number_of_components"]
    R = par_physical["gas_constant"]            # [J/mol/K]
    
    
    A = repmat(PhysicalData.PureComponent.liquidDensity(1,:),NS,1);
    B = repmat(PhysicalData.PureComponent.liquidDensity(2,:),NS,1);
    C = repmat(PhysicalData.PureComponent.liquidDensity(3,:),NS,1);
    D = repmat(PhysicalData.PureComponent.liquidDensity(4,:),NS,1);
    
    
    MW_0 = PhysicalData.PureComponent.molecularWeight/1000;% [kg/mol]
    rhoL_0 = A_rho./(B_rho.^(1+(1-repmat(T,1,NC)./C_rho).^D_rho))*1000; % [mol/m3]
    MW = x*(MW_0');
    
    
    rhoL = MW.*sum(x./rhoL_0,2).^(-1); % kg/m3
    if PhysicalData.VaporLiquidEquilibrium.modelType(2)==1
        [~,~,~,~,~,~,~,ZV] = f_fugacity_coefficient_SRK(PhysicalData,T,P,x);
    else
        ZV = 1;
    end
    
    return density_liquid


   # NS = size(T,1);
   # NC = PhysicalData.PureComponent.number;
   
   
   
   
   # T_c = repmat(PhysicalData.PureComponent.criticalTemperature,NS,1);
   # A_sig = repmat(PhysicalData.PureComponent.surfaceTension(1,:),NS,1);
   # B_sig = repmat(PhysicalData.PureComponent.surfaceTension(2,:),NS,1);
   # C_sig = repmat(PhysicalData.PureComponent.surfaceTension(3,:),NS,1);
   # D_sig = repmat(PhysicalData.PureComponent.surfaceTension(4,:),NS,1);
   
   
   # T_ = repmat(T,1,NC);
   # T_r = T_./T_c;
   
   
   
   
   # sigma = sum(x.*A_sig.*(1-T_r).^(B_sig+C_sig.*T_r+D_sig.*T_r.^2),2);