import pandas as pd

def load_pure_comp_properties(components):
    """
    Load pure component properties from CSV database file.
    The CSV file should have a column "component" with component names,
    and other columns with properties.
    Parameters
    ----------
    components : list of str
        List of component names to load properties for.
    Returns
    -------
    dict
        Dictionary with component properties arrays.
    """
    def try_data_clean_up(col):
        try:
            # pd.to_numeric will raise if conversion fails (no errors arg)
            return (pd.to_numeric(col))
        except (ValueError, TypeError):
            # return original (as string-ish) if conversion can't be done
            cleaned_string = col.astype(str).str.strip() if hasattr(col, "astype") else pd.Series(col)
            return cleaned_string
    
    data_file = "src/sep2p/data/pure_comp_data.csv" # So far, hardcoded path is used here
    db = pd.read_csv(data_file).set_index("component")
    db_components = db.loc[components]

    for col in db_components.columns:
        db_components[col] = try_data_clean_up(db_components[col])
    
    return {
        "components": db_components.index.to_numpy(),
        **{col: db_components[col].to_numpy() for col in db_components.columns}
    }


def load_liquid_mixture_model(components, model_liquid="ideal", model_gas="ideal"):
    """
    Load liquid mixture model parameters from CSV database files.
    Currently supports "unifac1p" model or "ideal".
    Parameters
    ----------
    components : list of str
        List of component names in the mixture.
    model_liquid : str, optional
        Liquid mixture model to use. Default is "ideal".
    model_gas : str, optional
        Gas mixture model to use. Default is "ideal".
    Returns
    -------
    dict
        Dictionary with liquid mixture model parameters.
    """
    param_liquid = {"model_liquid": "ideal"}
    
    if (model_liquid.lower()=="unifac1p"):        
        # Select only relevant components in "_nu"
        data_file = "src/sep2p/data/mixture_liquid_unifac1p_nu.csv" # So far, hardcoded path is used here
        db = pd.read_csv(data_file).set_index("component")
        db = db.apply(pd.to_numeric, errors="coerce").astype("Int64")
        db = db.fillna(0).astype(int)
        db_components_nu = db.loc[components]
        # mask = (db_components_nu != 0)
        # locations = mask.stack()[mask.stack()].index
        
        subgroups = (db_components_nu != 0).any(axis=0)
        nu = (db_components_nu.loc[:, subgroups].to_numpy()).T # shape(NG,NC)
        
        # Select only relevant subgroups in "_rq"
        data_file = "src/sep2p/data/mixture_liquid_unifac1p_rq.csv" # So far, hardcoded path is used here
        db = pd.read_csv(data_file, dtype={'subgroup': str}).set_index("subgroup")
        rq = (db.loc[subgroups].to_numpy())
        r = rq[:, 0]
        q = rq[:, 1]
                
        #  Relevant subgroups gives relevant main groups in "_main2sub"
        data_file = "src/sep2p/data/mixture_liquid_unifac1p_main2sub.csv" # So far, hardcoded path is used here
        db = pd.read_csv(data_file, dtype={'subgroup': str}).set_index("subgroup")
        maingroups = (db.loc[subgroups]["maingroup"]).astype(str)
        
        # Select relevant maingroups in "_aij"
        data_file = "src/sep2p/data/mixture_liquid_unifac1p_aij.csv" # So far, hardcoded path is used here
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

