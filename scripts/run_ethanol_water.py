import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import fsolve
#from scipy.integrate import cumulative_trapezoid
from scipy.integrate import solve_ivp

from sep2p import system_parameters
from sep2p import mixture_prop
from sep2p import column_models
from sep2p import plots

def run_ethanol_water_simulation() -> None:
    system = system_parameters.SystemParameters(
        ["Ethanol", "Water", "2-Methyl-1-propanol", "2-Butanol"],
        model_liquid="unifac1p"
        )

    system.column = {
        'number_of_stages':         30,  #[-]
        'stage_feed':               [22],  #[-]
        'constant_liquid':          1*3/0.1**1.5,  #[1/s]
        'constant_vapour':          1.5*3*1/(700/101325)**0.5,  #[1/s]
        'weir_height':              0.05,  #[m]
        'weir_crossarea':           (1.0)**2*np.pi/4  #[m^2]
    }


    NS = system.column['number_of_stages']
    u_test = {
        'stage_valve_opening_liquid':np.ones(NS),
        'stage_valve_opening_vapour':np.ones(NS),
        'condenser_reflux_ratio':   4.0,
        'condenser_pressure':       3*1e5,          # [Pa]
        'reboiler_duty':            0.15*160e3,  #[W]
        'feed_flow_rate':           np.array([3]),  # [mol/s]
        'feed_composition':         np.array([0.05, 0.95, 0, 0]),
        'feed_enthalpy':            mixture_prop.enthalpy_liquid_mixture(system.parameters, 405, 301325, np.array([0.05, 0.95, 0, 0])),
        }

    X_test = column_models.generate_initial_column_state(system, u_test)

    #X_test = utils.run_until_steady_state(column_models.column_ode, X_test, u_test, system, t_final=30000)
    t_final = 10*60*60
    sol_ss = solve_ivp(column_models.column_ode, [0, t_final], X_test, args=(u_test, system), method='BDF')
    #fsolve(lambda x: column_models.column_ode(0, x, u_test, system), sol_ss.y[:,-1])
    P_ss, T_ss, x_ss, y_ss, L_ss, V_ss, MT_L_ss, D_ss, B_ss, e_PID1_ss = column_models.sample_column(sol_ss, u_test, system)

    u_new = u_test.copy()
    u_new['feed_composition'] = np.array([0.05,0.95-1e-4-1e-4, 1e-4, 1e-4]),
    sol_new = solve_ivp(column_models.column_ode, [0, 10*60*60], sol_ss.y[:, -1], args=(u_new, system), method='BDF')
    P_new, T_new, x_new, y_new, L_new, V_new, MT_L_new, D_new, B_new, e_PID1_new = column_models.sample_column(sol_new, u_new, system)

    
    NC = system.parameters['number_of_components']
    nF = system.column['stage_feed'][0]

    trays = [str(i) for i in range(NS)]
    trays[0] = 'cnd'
    trays[nF] = 'nF'
    trays[NS-1] = 'rbl'

    plt.figure()
    for i in range(NS):
        plt.plot(sol_ss.t, T_ss[i, :], label = trays[i])
    plt.title('Temperature [K]')
    plt.ylabel('Temperature [K]')
    plt.xlabel('Time [s]')
    plt.legend()
    plt.show()
    
    fig, ax = plt.subplots(figsize=(6, 4))
    key_comp = range(2, NC)
    for i in range(NS):
        ax.plot(sol_new.t/3600, np.sum(x_new[i, :, 2:NC], axis=1), label = trays[i])
        if any(np.sum(x_new[i, :, key_comp], axis=0) > 0.2*np.max(np.sum(x_new[:, :, key_comp], axis=2))):
            ax.text(sol_new.t[-3]/3600, np.sum(x_new[i, -3, key_comp], axis=0), 'Tray '+str(trays[i]), 
                    va='center', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none')
            )
    #ax.set_title('Title')
    #ax.set_xlim([0, 1])
    #ax.set_ylim([0, 1])
    ax.set_xlabel('Time [h]')
    ax.set_ylabel('Molefraction fusel [-]')
    #ax.grid(True)
    plt.tight_layout()
    plt.show()

    plots.plot_working_lines(system, u_test['condenser_pressure'], x_ss[:, -1, :], y_ss[:, -1, :])



if __name__ == "__main__":
    run_ethanol_water_simulation()