# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 19:55:54 2025

@author: biss3
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import cumulative_trapezoid
from scipy.integrate import solve_ivp

from system_parameters import SystemParameters
from mixture_prop import enthalpy_liquid_mixture
from column_models import column_ode, sample_column


system = SystemParameters(
    ["Ethanol", "Water", "2-Methyl-1-propanol", "1-Butanol"], 
    model_liquid="unifac1p"
    )

system.column = {
    'number_of_stages':         40,  #[-]
    'stage_feed':               [11],  #[-]
    'constant_liquid':          1/0.2**1.5,  #[1/s]
    'constant_vapour':          1/(700/101325)**0.5,  #[1/s]
    'weir_height':              0.05,  #[m]
    'weir_crossarea':           (1.0)**2*np.pi/4,  #[m^2]
    'stage_height':             0.4  #[m]
}



u_test = {
     'stage_valve_opening_liquid':np.ones(40),
     'stage_valve_opening_vapour':np.ones(40),
     'condenser_reflux_ratio':  4,
     'condenser_pressure':      3*1e5,          # [Pa]
     'reboiler_duty':           160e3,  #[W]
     'feed_flow_rate':          np.array([3]),  # [mol/s]
     'feed_composition':        np.array([0.05,0.95-0*1e-4-0*1e-4, 0*1e-4, 0*1e-4]),
     'feed_enthalpy':           enthalpy_liquid_mixture(system.parameters, 365, 301325, np.array([0.05,0.95-1e-4-1e-4, 1e-4, 1e-4])),
     }


M_L_test = 0.5*(system.column["weir_crossarea"]*0.05*900/0.018)*np.ones((40, 4))
T_test = np.ones(40)*400
X_test = np.concat((np.concat((np.reshape(M_L_test, (40*4,)), T_test)), np.array([0.])))



t_final = 5*60*60
sol = solve_ivp(column_ode, [0, t_final], X_test, args=(u_test, system), method='BDF')

sol_ss = sol
P_ss, T_ss, x_ss, y_ss, L_ss, V_ss, MT_L_ss, D_ss, B_ss, e_PID1_ss = sample_column(sol_ss, u_test, system)


Nt = len(sol_ss.t)
NC = system.parameters['number_of_components']
NS = system.column['number_of_stages']
nF = system.column['stage_feed'][0]
NSC = NS*NC

trays = [str(i) for i in range(NS)]
trays[0] = 'cnd'
trays[nF] = 'nF'
trays[NS-1] = 'rbl'

plt.figure()
for i in range(NS):
    plt.plot(sol_ss.t, P_ss[i, :], label = trays[i])
plt.title('Pressure')
plt.ylabel('[bar]')
plt.legend()


plt.figure()
plt.plot(x_ss[:,-1,0], y_ss[:,-1,0],'o')
plt.show()


from equilibrium_calculations import solve_flash
x1 = np.linspace(0.0, 1.0, num=150)
x = np.hstack((x1[:, None], 1 - x1[:, None]))
P = np.atleast_1d(P_ss[0,-1])
system2 = SystemParameters(
    ["Ethanol", "Water"], 
    model_liquid="unifac1p"
    )
T=np.zeros(150)
y=np.zeros_like(x)
for i in range(150):
    Ti, yi = solve_flash(system2.parameters, pressure=P, composition_liquid=x[i,:])
    T[i]=Ti
    y[i,:]=yi
    
y1 = y[:,0]

plt.figure(1)
plt.plot(x1, y1, linestyle='-', color='b', label='Data')
plt.plot(x1, x1, linestyle='-', color='b', label='1')
plt.plot(x_ss[:,-1,0], y_ss[:,-1,0],'o', label="column")
plt.legend()
plt.show()
