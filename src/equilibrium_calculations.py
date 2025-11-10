# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 09:52:07 2025

@author: biss3
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from phase_sep.mixture_prop import vapour_liquid_equilibrium_constant
from phase_sep.pure_comp_prop import saturation_pressure

def solve_flash(parameters, temperature=0, pressure=0, composition_liquid=0, composition_gas=0, composition_overall=0):
    if ((temperature==0) and (composition_gas==0)):
        # Given pressure and liquid composition
        # Calculate temperature and gas compositions
        P = np.atleast_1d(pressure)
        x = np.atleast_2d(composition_liquid)
        NS, NC = x.shape
        
        Tbp = parameters["temperature_normal_boiling"][None, :]
        T_init = np.sum(Tbp*x, axis=1)
        X_init = T_init
                
        def solve_flash_objective(X, parameters, P, x):
            T_sol = X[0:NS]
            
            K_calc, y_calc = vapour_liquid_equilibrium_constant(parameters, T_sol, P, x)

            objective = 1 - np.sum(K_calc*x, axis=1)
            
            return objective
        
        sol = fsolve(solve_flash_objective, X_init, args=(parameters, P, x))
        T = np.array(sol[0:NS])

        y = vapour_liquid_equilibrium_constant(parameters, T, P, x)[1]
        return T, y
    
    elif ((pressure==0) and (composition_gas==0)):
       # Given temperature and liquid composition
       # Calculate pressure and gas compositions
       T = np.atleast_1d(temperature)
       x = np.atleast_2d(composition_liquid)
       NS, NC = x.shape
       
       Tbp = parameters["temperature_normal_boiling"][None, :]
       Pbp = saturation_pressure(parameters, T)
       P_init = np.sum(Pbp*x, axis=1)
       X_init = P_init
               
       def solve_flash_objective(X, parameters, T, x):
           P_sol = X[0:NS]
           
           K_calc, y_calc = vapour_liquid_equilibrium_constant(parameters, T, P_sol, x)

           objective = 1 - np.sum(K_calc*x, axis=1)
           
           return objective
       
       sol = fsolve(solve_flash_objective, X_init, args=(parameters, T, x))
       P = np.array(sol[0:NS])

       y = vapour_liquid_equilibrium_constant(parameters, T, P, x)[1]
       return P, y
    return 1


def binary_phase_diagram_Ty(parameters, pressure):
    x1 = np.linspace(0.0, 1.0, num=150)
    x = np.hstack((x1[:, None], 1 - x1[:, None]))
    P = np.atleast_1d(pressure)
    
    T = np.zeros(150)
    y = np.zeros_like(x)
    for i in range(150):
        Ti, yi = solve_flash(parameters, pressure=P, composition_liquid=x[i,:])
        T[i] = Ti
        y[i, :] = yi
        
    y1 = y[:,0]
    
    plt.figure(1)
    plt.plot(x1, y1, linestyle='-', color='b', label='Data')
    plt.plot(x1, x1, linestyle='-', color='b', label='1')
    plt.show()
    
    plt.figure(2)
    plt.plot(x1, T, linestyle='-', color='b', label='Liquid')
    plt.plot(y1, T, linestyle='-', color='r', label='Vapour')
    plt.legend()
    plt.show()

def binary_phase_diagram_Py(parameters, temperature):
    x1 = np.linspace(0.0, 1.0, num=150)
    x = np.hstack((x1[:, None], 1 - x1[:, None]))
    T = np.atleast_1d(temperature)
    
    P = np.zeros(150)
    y = np.zeros_like(x)
    
    for i in range(150):
        Pi, yi = solve_flash(parameters, temperature=T, composition_liquid=x[i,:])
        P[i] = Pi
        y[i, :] = yi
        
    y1 = y[:,0]
    
    plt.figure(1)
    plt.plot(x1, y1, linestyle='-', color='b', label='Data')
    plt.plot(x1, x1, linestyle='-', color='b', label='1')
    plt.show()
    
    plt.figure(2, figsize=(5, 4), dpi=80)
    plt.plot(x1, P/100000, linestyle='-', color='b', label='Liquid', linewidth=3)
    plt.plot(y1, P/100000, linestyle='-', color='r', label='Vapour', linewidth=3)
    plt.legend()
    plt.ylabel("Pressure [bar]")
    plt.xlabel("Molefraction Acetonitrile")
    plt.title("Temperature 318K")
    plt.show()