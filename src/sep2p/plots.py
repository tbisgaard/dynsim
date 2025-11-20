import numpy as np
import matplotlib.pyplot as plt

from sep2p import equilibrium_calculations
from sep2p import mixture_prop

def plot_working_lines(system, pressure, composition_liquid, composition_gas, component1=0, component2=1):
    if isinstance(component1, str):
        key_light = system.parameters["components"].index(component1)
        key_heavy = system.parameters["components"].index(component2)
    elif isinstance(component1, int):
        key_light = component1
        key_heavy = component2
    else:
        raise ValueError("component1 and component2 must be str or int")
        key_light = int(0)
        key_heavy = int(1)
    x = composition_liquid
    y = composition_gas
    P = pressure
    NS = x.shape[0]
    x1_wl = np.hstack([x[0, key_light], x[0:-1, key_light], x[-1, key_light]])
    y1_wl = np.hstack([x[0, key_light], y[1:, key_light], x[-1, key_light]])
    # Equilibrium
    x1_eq, y1_eq = equilibrium_calculations.generate_binary_phase_equilibrium_data(system.parameters, P, component1=component1, component2=component2, num_points=150)[0:2]
    # Stepping
    x1_step = np.zeros(2*NS + 1)
    y1_step = np.zeros(2*NS + 1)
    for i in range(0, NS):
        j = 2*i
        x1_step[j] = x1_wl[i]
        x1_step[j+1] = x1_wl[i+1]
        y1_step[j] = np.interp(x1_wl[i+1], x1_eq, y1_eq)
        y1_step[j+1] = y1_step[j]
    x1_step[-1] = x[-1, key_light]
    y1_step[-1] = x[-1, key_light]
    # Plot
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x1_step, y1_step, linewidth=2)
    ax.plot(x1_wl, y1_wl, linewidth=2)
    ax.plot(x1_eq, y1_eq, linewidth=2)
    #ax.set_title('Title')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    plt.tight_layout()
    plt.show()

def plot_binary_phase_diagram_Ty(parameters, pressure):
    x1 = np.linspace(0.0, 1.0, num=150)
    x = np.hstack((x1[:, None], 1 - x1[:, None]))
    P = np.atleast_1d(pressure)
    
    T = np.zeros(150)
    y = np.zeros_like(x)
    for i in range(150):
        Ti, yi = equilibrium_calculations.solve_flash(parameters, pressure=P, composition_liquid=x[i,:])
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

def plot_binary_phase_diagram_Py(parameters, temperature):
    x1 = np.linspace(0.0, 1.0, num=150)
    x = np.hstack((x1[:, None], 1 - x1[:, None]))
    T = np.atleast_1d(temperature)
    
    P = np.zeros(150)
    y = np.zeros_like(x)
    
    for i in range(150):
        Pi, yi = equilibrium_calculations.solve_flash(parameters, temperature=T, composition_liquid=x[i,:])
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

def plot_trace_component_K_factor(parameters, pressure, component1=0, component2=1):
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

    num_points = 150
    x1, y1, T = equilibrium_calculations.generate_binary_phase_equilibrium_data(parameters, pressure, component1=component1, component2=component2, num_points=num_points)
    x_trace = 1e-5
    x = x_trace*np.ones((num_points, NC))
    number_of_trace_components = NC - 2
    x[:, key_light] = x1
    x[:, key_heavy] = 1 - x1 - number_of_trace_components*x_trace

    Keq = mixture_prop.vapour_liquid_equilibrium_constant(parameters, T, pressure, x)[0]
    alpha_heavy = Keq/Keq[:, key_heavy][:, None]
    
    fig, ax = plt.subplots(figsize=(6, 4))
    for i in range(NC):
        ax.plot(x1, np.log(Keq[:, i]), label = str(parameters['components'][i]))
    #ax.set_title('Title')
    ax.set_xlim([0, 1])
    #ax.set_ylim([0, 1])
    ax.legend()
    ax.set_xlabel('Molefraction ' + parameters['components'][key_light])
    ax.set_ylabel('Log. K-factor')
    #ax.grid(True)
    plt.tight_layout()
    plt.show()