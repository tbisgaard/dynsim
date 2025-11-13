# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 10:22:58 2025

@author: biss3
"""

import numpy as np

def unifac(parameters, temperature, composition):
    """
    Original UNIFAC (UNIfied Functional-group Activity Coefficient) model.
    References:
        1-parameter: Aa. Fredenslund, R. Jones, J.M. Prausnitz, AIChE J., 1975,
        21, 1086
        2-parameter: H.K. Hansen, B. Coto, B. Kuhlmann, UNIFAC with linearly
        temperature-dependent group-interaction parameters. IVC-SEP International
        Report SEP 9212, Institut for Kemiteknik, Technical University of
        Denmark, 1992.

    Parameters
    ----------
    T : array-like (NS,)
        Temperatures
    x : ndarray (NS x NC)
        Mole fractions
    PhysicalData : object or dict
        Must contain VaporLiquidEquilibrium with attributes/dicts:
        - nu  : (NG x NC) subgroup coefficient matrix
        - RQ  : (NG x 2) array, RQ[:,0]=R (group volume), RQ[:,1]=Q (group surface area)
        - aij : (NG x NG) interaction parameter matrix
        - bij : (NG x NG) temperature-dependent interaction parameter matrix
    """

    nu = parameters["stoichiometry"]            # (NG, NC)
    R = parameters["size_parameter_r"]          # (NG,)
    Q = parameters["size_parameter_q"]          # (NG,)
    a = parameters["interaction_parameter"]     # (NG, NG)
    b = np.zeros_like(a)                        # (NG, NG)

    x = np.atleast_2d(np.array(composition, dtype=float)) # (NS, NC)
    T = np.atleast_1d(temperature).astype(float)          # (NS,)
    
    NS, NC = x.shape
    NG = nu.shape[0]

    # Expand T
    T2 = np.concatenate([T, np.tile(T, NC*NS//len(T))])  # (NS+NSC,)
    NSC = NS*NC

    Z = 10.0  # Coordination number

    # --- Combinatorial term (entropic) ---
    r = R @ nu                                  # (NC,)
    q = Q @ nu                                  # (NC,)
    r_tot = x @ r                               # (NS,)
    q_tot = x @ q                               # (NS,)
    phi_over_x = np.outer(1.0/r_tot, r)         # (NS, NC)
    theta_over_x = np.outer(1.0/q_tot, q)       # (NS, NC)
    theta_over_phi = theta_over_x/phi_over_x    # (NS, NC)
    l = Z/2*(r - q) - (r - 1)                   # (NC,)
    sum_xl = x @ l                              # (NS,)
    
    lngamma_comb = (                            # (NS, NC)
        np.log(phi_over_x) +
        Z/2*np.tile(q, (NS,1))*np.log(theta_over_phi) +
        np.tile(l, (NS,1)) -
        phi_over_x*sum_xl[:, None]
        )
    
    # --- Residual term (enthalpic) ---
    # All new
    def ln_gamma_calc(temperature, composition):
        T = temperature
        x = composition
        
        X = (x @ nu.T) / np.sum(x @ nu.T, keepdims=True, axis=1)   # (NS, NG)
        Theta = X*Q / np.sum(X*Q, keepdims=True, axis=1)# (NS, NG)
        
        T0 = 298.15
        delta_a = a - np.diag(a)[:, None]           # (NG, NG)
        delta_b = b - np.diag(b)[:, None]           # (NG, NG)
        Psi = np.exp(-(                             # (NS, NG, NG)
            delta_a[None, :, :] +
            delta_b[None, :, :] * (T[:, None, None] - T0)
        ) / T[:, None, None])
            
        term_1 = np.einsum('ij,ijk->ik', Theta, Psi) # (NS, NG)
        term_2 = np.sum(                            # (NS, NG)
            (Theta[:, None, :] * Psi) /
            np.einsum('in, inm -> im', Theta, Psi)[:, None, :],
            axis=2
            )
        lnGamma = Q*(1 - np.log(term_1) - term_2)
        return lnGamma
    
    
    x_ref = np.tile(np.diag(np.ones(NC)), (NS,1))
    T_ref = np.tile(T,NC)
    
    ln_gamma_group = ln_gamma_calc(T, x)    # (NS, NG)
    ln_gamma_ref = ln_gamma_calc(T_ref, x_ref).reshape((NS, NC, NG))
    ln_gamma_res = np.einsum('nj, in -> ij', nu, ln_gamma_group) - np.einsum('nj, ijn -> ij', nu, ln_gamma_ref)#nu*(ln_gamma_group - ln_gamma_ref)
    gamma = np.exp(lngamma_comb + ln_gamma_res)
    
    
    return gamma