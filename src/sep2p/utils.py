# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 21:08:00 2025

@author: biss3
"""
import numpy as np

def step_function(t):
    # From 0 to 1 after in t inteval from 0 to 1
    t = (t>0)*t
    f = ((1-np.exp(-(t+0)/0.2))**3) * ((1-np.exp(-(t+0.95)/0.1))**2)
    return f
