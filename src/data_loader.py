# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 09:02:32 2025

@author: biss3
"""

import os
import pandas as pd


def load_pure_comp_properties(components):
    def try_data_clean_up(col):
        try:
            # pd.to_numeric will raise if conversion fails (no errors arg)
            return (pd.to_numeric(col))
        except (ValueError, TypeError):
            # return original (as string-ish) if conversion can't be done
            cleaned_string = col.astype(str).str.strip() if hasattr(col, "astype") else pd.Series(col)
            return cleaned_string
    
    root_folder = "C:\\Users\\biss3\\OneDrive\\Python\\phase_sep_project\\"
    data_file = os.path.join(root_folder, "phase_sep\\data","pure_comp_data.csv")
    db = pd.read_csv(data_file).set_index("component")
    db_components = db.loc[components]

    for col in db_components.columns:
        db_components[col] = try_data_clean_up(db_components[col])
    
    return {col: db_components[col].to_numpy() for col in db_components.columns}


def load_liquid_mixture_model(components, model_liquid="ideal", model_gas="ideal"):
    param_liquid = {"model_liquid": "ideal"}
    
    if (model_liquid.lower()=="unifac1p"):
        root_folder = "C:\\Users\\biss3\\OneDrive\\Python\\phase_sep_project\\"
        
        # Select only relevant components in "_nu"
        data_file = os.path.join(root_folder, "phase_sep\\data","mixture_liquid_unifac1p_nu.csv")
        db = pd.read_csv(data_file).set_index("component")
        db = db.apply(pd.to_numeric, errors="coerce").astype("Int64")
        db = db.fillna(0).astype(int)
        db_components_nu = db.loc[components]
        # mask = (db_components_nu != 0)
        # locations = mask.stack()[mask.stack()].index
        
        subgroups = (db_components_nu != 0).any(axis=0)
        nu = (db_components_nu.loc[:, subgroups].to_numpy()).T # shape(NG,NC)
        
        # Select only relevant subgroups in "_rq"
        data_file = os.path.join(root_folder, "phase_sep\\data","mixture_liquid_unifac1p_rq.csv")
        db = pd.read_csv(data_file, dtype={'subgroup': str}).set_index("subgroup")
        rq = (db.loc[subgroups].to_numpy())
        r = rq[:, 0]
        q = rq[:, 1]
                
        #  Relevant subgroups gives relevant main groups in "_main2sub"
        data_file = os.path.join(root_folder, "phase_sep\\data","mixture_liquid_unifac1p_main2sub.csv")
        db = pd.read_csv(data_file, dtype={'subgroup': str}).set_index("subgroup")
        maingroups = (db.loc[subgroups]["maingroup"]).astype(str)
        
        # Select relevant maingroups in "_aij"
        data_file = os.path.join(root_folder, "phase_sep\\data","mixture_liquid_unifac1p_aij.csv")
        db = pd.read_csv(data_file, dtype={'maingroup': str}).set_index("maingroup")
        aij = db.loc[maingroups, maingroups].to_numpy()
        
        param_liquid = {
            "model_liquid":             "UNIFAC1p",
            "stoichiometry":            nu,
            "size_parameter_r":         r,
            "size_parameter_q":         q,
            "interaction_parameter":    aij
            }

    return param_liquid

