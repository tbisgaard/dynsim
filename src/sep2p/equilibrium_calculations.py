import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from sep2p import mixture_prop
from sep2p import pure_comp_prop

def solve_flash(parameters, temperature=0, pressure=0, composition_liquid=0, composition_gas=0, composition_overall=0):
    """
    Solves flash calculations for binary and multicomponent mixtures.
    
    Depending on which variables are set to zero, the function calculates the missing
    variables based on the provided inputs.

    UNDER DEVELOPMENT....

    Parameters
    ----------
    temperature : array-like (NS,)
        Temperatures
    pressure : array-like (NS,)
        Pressures
    composition_liquid : ndarray (NS x NC) 
        Mole fractions in liquid phase
    composition_gas : ndarray (NS x NC)
        Mole fractions in gas phase
    composition_overall : ndarray (NS x NC)
        Overall mole fractions
    Returns
    -------
    Depending on the input variables set to zero, returns the calculated variables.
    - If temperature and composition_gas are zero: returns temperature and composition_gas
    - If pressure and composition_gas are zero: returns pressure and composition_gas
    """
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
            
            K_calc, y_calc = mixture_prop.vapour_liquid_equilibrium_constant(parameters, T_sol, P, x)

            objective = 1 - np.sum(K_calc*x, axis=1)
            
            return objective
        
        sol = fsolve(solve_flash_objective, X_init, args=(parameters, P, x))
        T = np.array(sol[0:NS])

        y = mixture_prop.vapour_liquid_equilibrium_constant(parameters, T, P, x)[1]
        return T, y
    
    elif ((pressure==0) and (composition_gas==0)):
       # Given temperature and liquid composition
       # Calculate pressure and gas compositions
       T = np.atleast_1d(temperature)
       x = np.atleast_2d(composition_liquid)
       NS, NC = x.shape
       
       Tbp = parameters["temperature_normal_boiling"][None, :]
       Pbp = pure_comp_prop.saturation_pressure(parameters, T)
       P_init = np.sum(Pbp*x, axis=1)
       X_init = P_init
               
       def solve_flash_objective(X, parameters, T, x):
           P_sol = X[0:NS]
           
           K_calc, y_calc = mixture_prop.vapour_liquid_equilibrium_constant(parameters, T, P_sol, x)

           objective = 1 - np.sum(K_calc*x, axis=1)
           
           return objective
       
       sol = fsolve(solve_flash_objective, X_init, args=(parameters, T, x))
       P = np.array(sol[0:NS])

       y = mixture_prop.vapour_liquid_equilibrium_constant(parameters, T, P, x)[1]
       return P, y
    return 1

def generate_binary_phase_equilibrium_data(parameters, pressure, component1=0, component2=1, num_points=150):
    """
    Generates binary phase equilibrium data for given components at specified pressure.
    Parameters
    ----------
    parameters : dict
        System parameters including component information.
    pressure : float
        Pressure at which to generate the phase equilibrium data.
    component1 : int or str, optional
        Index or name of the first component. Default is 0.
    component2 : int or str, optional
        Index or name of the second component. Default is 1.
    num_points : int, optional
        Number of data points to generate. Default is 150.
    Returns
    -------
    x1 : ndarray (num_points,)
        Mole fractions of component1 in the liquid phase.
    y1 : ndarray (num_points,)
        Mole fractions of component1 in the vapour phase.
    T : ndarray (num_points,)
        Temperatures corresponding to the phase equilibrium data.   
    """
    if isinstance(component1, str):
        key_light = parameters["components"].index(component1)
        key_heavy = parameters["components"].index(component2)
    elif isinstance(component1, int):
        key_light = component1
        key_heavy = component2
    else:
        raise ValueError("component1 and component2 must be str or int")
        key_light = int(0)
        key_heavy = int(1)
    NC = len(parameters["components"])

    x1 = np.linspace(0.0, 1.0, num=num_points)
    x = np.zeros((num_points, NC))#np.hstack((x1[:, None], 1 - x1[:, None]))
    x[:, key_light] = x1
    x[:, key_heavy] = 1 - x1
    P = np.atleast_1d(pressure)
    T = np.zeros(num_points)
    y = np.zeros_like(x)
    for i in range(150):
        Ti, yi = solve_flash(parameters, pressure=P, composition_liquid=x[i,:])
        T[i] = Ti
        y[i, :] = yi
        
    y1 = y[:,0]
    return x1, y1, T