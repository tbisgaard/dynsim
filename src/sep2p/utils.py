import numpy as np

def step_function(t):
    """
    Smooth step function from 0 to 1 in time interval from 0 to 1
    
    Parameters
    ----------
    t : float or array-like
        Time
    Returns
    -------
    f : float or array-like
        Function value
    """
    # Smooth transition from 0 to 1 after in t inteval from 0 to 1
    t = (t>0)*t
    f = ((1-np.exp(-(t+0)/0.2))**3) * ((1-np.exp(-(t+0.95)/0.1))**2)
    return f

def run_until_steady_state(odefun, states_initial, inputs_unit, system, t_final=1e3):
    """
    Run ODE function until steady state is reached
    Parameters
    ----------
    odefun : function
        ODE function
    states_initial : ndarray
        Initial states
    inputs_unit : dict
        Unit inputs
    system : SystemParameters
    t_final : float, optional
        Final time for each integration step. The default is 1e3.
    Returns
    -------
    states_final : ndarray
        Final states at steady state
    """
    from scipy.integrate import solve_ivp
    atol_ss = 1e-2
    states_start = states_initial
    states_final = states_initial + 2*atol_ss
    objective_ss = 2*atol_ss*np.ones_like(states_initial)
    timesteps = 0
    while (any(objective_ss > atol_ss) and (timesteps < 20)) :
        sol = solve_ivp(odefun, [0, t_final], states_start, args=(inputs_unit, system), method='BDF', rtol=1e-6, atol=1e-8)
        states_final = sol.y[:, -1]
        objective_ss = abs(states_final - states_start)
        states_start = states_final
        timesteps+=1
        print('Step='+str(timesteps))
    return states_final