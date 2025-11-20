import numpy as np

#from util_general import gas_constant, acceleration_constant, step_function
#from mixture_prop import molecular_mass_mixture, density_liquid_mixture, density_gas_mixture, enthalpy_gas_mixture, enthalpy_liquid_mixture, vapour_liquid_equilibrium_constant, const_pressure_heat_capacity_liquid_mixture
from sep2p import utils
from sep2p import phys_constants
from sep2p import mixture_prop
from sep2p import pure_comp_prop

def column(time, states, inputs_unit, system):
    """


    Column model ODEs for distillation column.
    NS is number of stages, NC is number of components.

    Parameters
    ----------
    ttime : float
        Time [s].
    states : Numpy array of length (NS*NC + NS + 1).
        Column state.
    inputs_unit : Dictionary
        Unit operation inputs.
    system : SystemParameters object
        System parameters.
    Returns
    -------
    dXdt : Numpy array of length (NS*NC + NS + 1).
        Column state derivatives.
    X : Numpy array of length (NS*NC + NS + 1).
        Column state. (For convenience)
    P : Numpy array of length NS.
        Stage pressures [Pa].
    T : Numpy array of length NS.
        Stage temperatures [K].
    x : Numpy array of shape (NS, NC).
        Stage liquid mole fractions [-].
    y : Numpy array of shape (NS, NC).
        Stage vapour mole fractions [-].
    L : Numpy array of shape (NS, NC).
        Stage liquid flow rates [mol/s].
    V : Numpy array of shape (NS, NC).
        Stage vapour flow rates [mol/s].
    MT_L : Numpy array of length NS.
        Stage total liquid holdups [mol].
    D : float
        Distillate flow rate [mol/s].
    B : float
        Bottoms flow rate [mol/s].
    """
    t = time
    X = states
    u = inputs_unit

    NS = system.column['number_of_stages']      # Number of stages
    nF = system.column['stage_feed']            # Feed stage number
    cL = system.column['constant_liquid']       # Liquid hydraulic constant [mol/s]
    cV = system.column['constant_vapour']       # Vapour hydraulic constant [mol/s]
    Hwr = system.column['weir_height']          # Weir height [m]
    A_T = system.column['weir_crossarea']       # Tray cross sectional area  [m^2]
    NC = system.parameters['number_of_components']# Number of components
    NSC = NS*NC

    alphaL = u['stage_valve_opening_liquid']    # Valve opening, 0..1 [-]
    alphaV = u['stage_valve_opening_vapour']    # Valve opening, 0..1 [-]
    RR = u['condenser_reflux_ratio']            # Reflux ratio [-]
    Pset_cnd = u['condenser_pressure']          # Condenser pressure set point [Pa]
    Q_rbl = u['reboiler_duty']                  # Reboiler duty [W]
    F_feed = u['feed_flow_rate']                # Feed flow rate [mol/s]
    hF_feed = u['feed_enthalpy']                # Feed enthalpy [J/mol]
    z_feed = u['feed_composition']              # Feed composition, j=1,..,NC
    # Pcnd = u.Pcnd
    # Pdiff = u.Pdiff
    
    # State variables
    M_L = np.reshape(X[0:NSC], [NS, NC])        # Liquid holdup NS x NC [mol]
    T = X[NSC:NSC+NS]                           # Temperature NS [K]
    e_int_PID1 = X[NSC+NS]
    # e_int_PID2 = X[NSC+NS+1]
    # e_int_PID3 = X[NSC+NS+2]
    
    # Simple calculations
    # Mole fractions
    MT_L = np.sum(M_L, 1)                       # Total stage liquid holdup [mol]
    x = M_L / np.tile(np.reshape(MT_L, (NS, 1)), (1, NC))# Liquid molefraction[-]
    
    # Physical properties
    rho_L = mixture_prop.density_liquid_mixture(system.parameters, T, x)
    # rho_V = density_gas_mixture(T, np.ones_like(T), x, components)
    MW_L = mixture_prop.molecular_mass_mixture(system.parameters, x)
    
    
    # VLE
    P_1bar = 101325*np.ones_like(T)             # Assume 1 bar for K-values [Pa]
    Keq_1bar = mixture_prop.vapour_liquid_equilibrium_constant(system.parameters, T, P_1bar, x)[0]
    P = P_1bar*np.sum(Keq_1bar * x, 1)          # Pressure [Pa]
    Keq = Keq_1bar*(P_1bar/P)[:, None]          # K-values at actual pressure
    y = Keq*x                                   # Vapour molefraction[-]
    y[0, :] = x[0, :]                           # Total condenser
    
    # Control loops
    # PID1: Qcnd->P0
    e_PID1 = (Pset_cnd - P[0])/101325
    P_PID1 = -0.0008*1.5
    I_PID1 = -0.0005*1.5
    u0_PID1 = -1e2
    u_PID1 = u0_PID1 + (P_PID1*e_PID1 + I_PID1*e_int_PID1)*101325*u0_PID1
    Q_cnd = u_PID1
    
    # PID2
    e_PID2 = 0
    
    #PID3
    e_PID3 = 0
    
    F = np.zeros(NS)
    F[nF] = F_feed
    z = np.zeros((NS, NC))
    z[nF, :] = z_feed
    hF = np.zeros(NS)
    hF[nF] = hF_feed
    
    # Tray hydraulics
    Li = np.zeros(NS)
    Vi = np.zeros(NS)
    H = np.zeros(NS)
    dP = np.zeros(NS)
    L = np.zeros((NS, NC))
    V = np.zeros((NS, NC))
    Pref = 101325  #[Pa]
    H[0] = 0
    H[1:NS-1] = MT_L[1:NS-1] * MW_L[1:NS-1] / rho_L[1:NS-1] / A_T  #[m]
    H[NS-1] = 0
    dHrel = H/Hwr - 1                       # Rel. driving force, weir hight [-]
    dP_loss = 0.5*rho_L*phys_constants.ACCELERATION_CONSTANT*H# Tray pressure drop [Pa]
    dP[0] = 0
    dP[1:NS] = P[1:NS] - P[0:NS-1] - dP_loss[0:NS-1]*(dP_loss[0:NS-1]>0)
    dPrel = dP/Pref                         # Rel. driving force, pressure [-]
    Li = cL * utils.step_function(50*dHrel) * abs(dHrel)**1.5 * alphaL**0.5
    Vi = cV * abs(dPrel)**0.5 * np.sign(dPrel) * alphaV**0.5 * utils.step_function(50*dPrel)
    
    # Perfect level control - Condenser Reflux drum
    Li[0] = (Vi[1] - Vi[0]) * RR / (RR + 1) * utils.step_function((Vi[1] - Vi[0])*100)
    D = Vi[1] - Vi[0] - Li[0]
    # Perfect level control - Reboiler
    Li[NS-1] = 0
    B = Li[NS-2] - Li[NS-1] - Vi[NS-1]
    
    # Heat transfer
    Q = np.zeros(NS)
    Q[0] = Q[0] + Q_cnd
    Q[NS-1] = Q[NS-1] + Q_rbl
    
    # Enthalpy
    h_L = mixture_prop.enthalpy_liquid_mixture(system.parameters, T, P, x)
    h_V = mixture_prop.enthalpy_gas_mixture(system.parameters, T, P, x)
    
    # ========== Conservation equations ==========
    dM_Ldt = np.zeros((NS, NC))
    dMhTdt = np.zeros(NS)
    
    L = np.tile(np.reshape(Li, (NS, 1)), NC)
    V = np.tile(np.reshape(Vi, (NS, 1)), NC)
    
    # Top stage
    i = 0
    dM_Ldt[i, :] = - L[i, :] * x[i, :] \
        - D * x[i, :] \
        + V[i+1, :] * y[i+1, :] \
        - V[i, :] * y[i, :]
    dMhTdt[i] = - Li[i] * h_L[i] \
        - D * h_L[i] \
        + Vi[i+1] * h_V[i+1] \
        - Vi[i] * h_V[i] \
        + Q[i]

    # Middle stages
    dM_Ldt[1:NS-1, :] = L[0:NS-2, :] * x[0:NS-2, :] \
        - L[1:NS-1, :] * x[1:NS-1, :] \
        + V[2:NS, :] * y[2:NS, :] \
        - V[1:NS-1, :] * y[1:NS-1, :]
    dMhTdt[1:NS-1] = Li[0:NS-2] * h_L[0:NS-2] \
        - Li[1:NS-1] * h_L[1:NS-1] \
        + Vi[2:NS] * h_V[2:NS] \
        - Vi[1:NS-1] * h_V[1:NS-1] \
        + Q[1:NS-1]
    
    # Bottom stages
    i = NS-1
    dM_Ldt[i, :] = L[i-1, :] * x[i-1, :] \
        - B * x[i, :] \
        - V[i, :] * y[i, :]
    dMhTdt[i] = Li[i-1] * h_L[i-1] \
        - B * h_L[i] \
        - Vi[i] * h_V[i] \
        + Q[i]
    
    # Feed stage correction
    i = nF
    dM_Ldt[i, :] = dM_Ldt[i, :] + F[i] * z[i, :]
    dMhTdt[i] = dMhTdt[i] + F[i] * hF[i]
    
    # Convert energy balance to temperature derivative
    hLdM_Ldt = mixture_prop.enthalpy_liquid_mixture(system.parameters, T, P, dM_Ldt)
    M_LCPL = mixture_prop.const_pressure_heat_capacity_liquid_mixture(system.parameters, T, P, M_L)
    dTdt = (dMhTdt - hLdM_Ldt) / (M_LCPL)
    
    
    dXdt = np.concat((
        np.concatenate( (
        np.reshape(dM_Ldt, NSC),
        dTdt
        ) ),
        np.array([e_PID1])))#, e_PID2, e_PID3
    
    return dXdt, X, P, T, x, y, Li, Vi, MT_L, D, B



def column_ode(time, states, inputs_unit, system):
    """
    Wrapper for column ODEs for use with ODE solvers.
    Parameters
    
    ----------
    time : float
        Time [s].
    states : Numpy array of length (NS*NC + NS + 1).
        Column state.
    inputs_unit : Dictionary
        Unit operation inputs.
    system : SystemParameters object
        System parameters.
    Returns
    -------
    dXdt : Numpy array of length (NS*NC + NS + 1).
        Column state derivatives.
    """
    return column(time, states, inputs_unit, system)[0]


def sample_column(sol, inputs_unit, system, idx=None):
    """


    Sample column solution from ODE solver output.
    Parameters
    ----------
    sol : ODE solver solution object
        ODE solver solution.
    inputs_unit : Dictionary
        Unit operation inputs.
    system : SystemParameters object
        System parameters.
    idx : List or array, optional
        Indices of time points to sample. The default is None, which samples all time points.
    Returns
    -------
    P : Numpy array of shape (NS, Nt).
        Stage pressures [Pa].
    T : Numpy array of shape (NS, Nt).
        Stage temperatures [K].
    x : Numpy array of shape (NS, Nt, NC).
        Stage liquid mole fractions [-].
    y : Numpy array of shape (NS, Nt, NC).
        Stage vapour mole fractions [-].
    L : Numpy array of shape (NS, Nt).
        Stage liquid flow rates [mol/s].
    V : Numpy array of shape (NS, Nt).
        Stage vapour flow rates [mol/s].
    MT_L : Numpy array of shape (NS, Nt).
        Stage total liquid holdups [mol].
    D : Numpy array of length Nt.
        Distillate flow rate [mol/s].
    B : Numpy array of length Nt.
        Bottoms flow rate [mol/s].
    e_PID1 : Numpy array of length Nt.
        PID1 error integral.
    """
    u = inputs_unit
    NC = system.parameters['number_of_components']
    NS = system.column['number_of_stages']
    nF = system.column['stage_feed']
    NSC = NS*NC

    Nt = len(sol.t)
    if idx==None:
        idx = range(Nt)
    else:
        Nt = len(idx)
    
    P = np.zeros([NS, Nt])
    T = np.zeros([NS, Nt])
    L = np.zeros([NS, Nt])
    V = np.zeros([NS, Nt])
    MT_L = np.zeros([NS, Nt])
    x = np.zeros([NS, Nt, NC])
    y  = np.zeros([NS, Nt, NC])
    D = np.zeros(Nt)
    B = np.zeros(Nt)
    e_PID1 = np.zeros(Nt)

    for k in idx:
        dXdtk, Xk, Pk, Tk, xk, yk, Lk, Vk, MT_Lk, Dk, Bk = \
            column(sol.t[k], sol.y[:, k], u, system)
        e_PID1[k] = sol.y[NSC+NS, k]
        P[:, k] = Pk
        T[:, k] = Tk
        L[:, k] = Lk
        V[:, k] = Vk
        MT_L[:, k] = MT_Lk
        D[k] = Dk
        B[k] = Bk
        x[:, k, :] = np.reshape(xk, (NS, NC))
        y[:, k, :] = np.reshape(yk, (NS, NC))
            
    return P, T, x, y, L, V, MT_L, D, B, e_PID1

def generate_initial_column_state(system, inputs_unit):
    """ 
    
    
    Generate an initial column state based on feed conditions and column parameters 
    Parameters
    ----------
    system : SystemParameters object
        System parameters.

    inputs_unit : Dictionary 
        Unit operation inputs.
    Returns
    -------
    X : Numpy array of length (NS*NC + NS + 1).
        Initial column state.

    """
    NC = system.parameters['number_of_components']
    NS = system.column['number_of_stages']
    NSC = NS*NC
    Hw = system.column['weir_height']
    Across = system.column['weir_crossarea']
    H = 0.03 # Height above weir [m]
    u = inputs_unit
    P = u['condenser_pressure']
    z = u['feed_composition']
    T_const = 373.15
    dHvap_const = pure_comp_prop.heat_of_vaporisation(system.parameters, T_const)
    Tnbp = system.parameters['temperature_normal_boiling']
    Tbp = (1/Tnbp - phys_constants.GAS_CONSTANT/dHvap_const*np.log(P/101325))**(-1)

    Tavg = np.sum(z*Tbp)

    rho_L = mixture_prop.density_liquid_mixture(system.parameters, Tavg, z)
    MW_L = mixture_prop.molecular_mass_mixture(system.parameters, z)
    MT_L = Hw*Across*rho_L/MW_L
    M_L = MT_L*z[None, :]*np.ones((NS, NC))
    T = np.ones(NS)*Tavg
    X = np.concat((np.concat((np.reshape(M_L, (NSC,)), T)), np.array([0.])))
    return X