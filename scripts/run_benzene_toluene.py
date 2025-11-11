# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 13:54:49 2025

@author: biss3
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import cumulative_trapezoid
from scipy.integrate import solve_ivp

from .mixture_prop import enthalpy_liquid_mixture
from .column_models import column_ode, sample_column




par_column = {
    'number_of_components':     2,
    'components':               ['Benzene', 'Toluene'],
    'number_of_stages':         20,  #[-]
    'stage_feed':               [10],  #[-]
    'constant_liquid':          1/0.2**1.5,  #[1/s]
    'constant_vapour':          1/(700/101325)**0.5,  #[1/s]
    'weir_height':              0.02,  #[m]
    'weir_crossarea':           (0.2)**2*np.pi/4,  #[m^2]
    'stage_height':             0.4  #[m]
}



u_test = {
     'stage_valve_opening_liquid':np.ones(20),
     'stage_valve_opening_vapour':np.ones(20),
     'condenser_reflux_ratio':  2.5,
     'condenser_pressure':      1*1e5,          # [Pa]
     'reboiler_duty':           5000,  #[W]
     'feed_flow_rate':          np.array([0.1]),  # [mol/s]
     'feed_composition':        np.array([0.5, 0.5]),
     'feed_enthalpy':           enthalpy_liquid_mixture(365, 101325, np.array([0.5, 0.5]),par_column["components"]),
     }


M_L_test = 0.5*(0.2**2*np.pi/4*0.02*800/0.08442353)*np.ones((20, 2))
T_test = np.ones(20)*350
X_test = np.concat((np.concat((np.reshape(M_L_test, (20*2,)), T_test)), np.array([0.])))



t_final = 2*60*60
sol = solve_ivp(column_ode, [0, t_final], X_test, args=(u_test, par_column), method='BDF')

sol_ss = sol
P_ss, T_ss, x_ss, y_ss, L_ss, V_ss, MT_L_ss, D_ss, B_ss, e_PID1_ss = sample_column(sol_ss, u_test, par_column)


Nt = len(sol_ss.t)
NC = par_column['number_of_components']
NS = par_column['number_of_stages']
nF = par_column['stage_feed'][0]
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
# f_visualise_profile(sol_ss,u_test,par_column, par_phys,plot_types=['V','L'])