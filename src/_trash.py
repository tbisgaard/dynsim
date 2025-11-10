# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 16:01:30 2025

@author: biss3
"""


# Load database
# filename = "C:\\Users\\biss3\\OneDrive\\Python\\phase_sep_project\\"
# data_file = os.path.join(filename, "phase_sep\\data","pure_comp_data.csv")
# db = pd.read_csv(data_file).set_index("component")

# component_short_cache = ((db["component_short"].iloc[1:])).to_dict()
# molecular_mass_cache = (pd.to_numeric(db["molecular_mass"].iloc[1:])).to_dict()
# omega_cache = (pd.to_numeric(db["molecular_mass"].iloc[1:])).to_dict()
# temperature_critical_cache = (pd.to_numeric(db["temperature_critical"].iloc[1:])).to_dict()
# pressure_critical_cache = (pd.to_numeric(db["pressure_critical"].iloc[1:])).to_dict()
# density_critical_cache = (pd.to_numeric(db["density_critical"].iloc[1:])).to_dict()
# temperature_normal_boiling_cache = (pd.to_numeric(db["temperature_normal_boiling"].iloc[1:])).to_dict()
# rg_cache = (pd.to_numeric(db["rg"].iloc[1:])).to_dict()
# dm_cache = (pd.to_numeric(db["dm"].iloc[1:])).to_dict()
# solubility_parameter_cache = (pd.to_numeric(db["solubility_parameter"].iloc[1:])).to_dict()
# vanderWals_volume_cache = (pd.to_numeric(db["vanderWals_volume"].iloc[1:])).to_dict()
# epsilon_cache = (pd.to_numeric(db["epsilon"].iloc[1:])).to_dict()
# enthalpy_of_formation_ideal_gas_cache = (pd.to_numeric(db["enthalpy_of_formation_ideal_gas"].iloc[1:])).to_dict()
# entropy_of_formation_ideal_gas_cache = (pd.to_numeric(db["entropy_of_formation_ideal_gas"].iloc[1:])).to_dict()
# saturation_pressure_constant_A_cache = (pd.to_numeric(db["saturation_pressure_constant_A"].iloc[1:])).to_dict()
# saturation_pressure_constant_B_cache = (pd.to_numeric(db["saturation_pressure_constant_B"].iloc[1:])).to_dict()
# saturation_pressure_constant_C_cache = (pd.to_numeric(db["saturation_pressure_constant_C"].iloc[1:])).to_dict()
# saturation_pressure_constant_D_cache = (pd.to_numeric(db["saturation_pressure_constant_D"].iloc[1:])).to_dict()
# saturation_pressure_constant_E_cache = (pd.to_numeric(db["saturation_pressure_constant_E"].iloc[1:])).to_dict()
# heat_of_vaporisation_constant_A_cache = (pd.to_numeric(db["heat_of_vaporisation_constant_A"].iloc[1:])).to_dict()
# heat_of_vaporisation_constant_B_cache = (pd.to_numeric(db["heat_of_vaporisation_constant_B"].iloc[1:])).to_dict()
# heat_of_vaporisation_constant_C_cache = (pd.to_numeric(db["heat_of_vaporisation_constant_C"].iloc[1:])).to_dict()
# heat_of_vaporisation_constant_D_cache = (pd.to_numeric(db["heat_of_vaporisation_constant_D"].iloc[1:])).to_dict()
# heat_capacity_ideal_gas_constant_A_cache = (pd.to_numeric(db["heat_capacity_ideal_gas_constant_A"].iloc[1:])).to_dict()
# heat_capacity_ideal_gas_constant_B_cache = (pd.to_numeric(db["heat_capacity_ideal_gas_constant_B"].iloc[1:])).to_dict()
# heat_capacity_ideal_gas_constant_C_cache = (pd.to_numeric(db["heat_capacity_ideal_gas_constant_C"].iloc[1:])).to_dict()
# heat_capacity_ideal_gas_constant_D_cache = (pd.to_numeric(db["heat_capacity_ideal_gas_constant_D"].iloc[1:])).to_dict()
# heat_capacity_ideal_gas_constant_E_cache = (pd.to_numeric(db["heat_capacity_ideal_gas_constant_E"].iloc[1:])).to_dict()
# molar_density_constant_A_cache = (pd.to_numeric(db["molar_density_constant_A"].iloc[1:])).to_dict()
# molar_density_constant_B_cache = (pd.to_numeric(db["molar_density_constant_B"].iloc[1:])).to_dict()
# molar_density_constant_C_cache = (pd.to_numeric(db["molar_density_constant_C"].iloc[1:])).to_dict()
# molar_density_constant_D_cache = (pd.to_numeric(db["molar_density_constant_D"].iloc[1:])).to_dict()
# surface_tension_constant_A_cache = (pd.to_numeric(db["surface_tension_constant_A"].iloc[1:])).to_dict()
# surface_tension_constant_B_cache = (pd.to_numeric(db["surface_tension_constant_B"].iloc[1:])).to_dict()
# surface_tension_constant_C_cache = (pd.to_numeric(db["surface_tension_constant_C"].iloc[1:])).to_dict()
# surface_tension_constant_D_cache = (pd.to_numeric(db["surface_tension_constant_D"].iloc[1:])).to_dict()





# db = db.apply(pd.to_numeric, errors="coerce")
# db_dict = db.to_dict(orient="index")
# db = db.apply(lambda s: pd.to_numeric(s.astype(str).str.replace(",", "").str.strip(), errors="ignore"))
# def get_component_values(data_cache, components):
#     if isinstance(components, str):
#         components = [components]
#     return np.array([data_cache[c] for c in components])

# FROM UNIFAC() in activity_coefficient_models

    # --- Residual term (enthalpic) ---
    # # Group mole fraction
    
    # X_group = (x @ nu.T) / np.sum(x @ nu.T, keepdims=True, axis=1)   # (NS, NG)
    # Theta_group = X_group*Q / np.sum(X_group*Q, keepdims=True, axis=1)# (NS, NG)
    
    # X_ref = nu / np.sum(nu, axis=0, keepdims=True)   # (NG, NC)
    # Theta_ref = Q[:, None]*X_ref / (Q @ X_ref)  # (NG, NC) 

    # Theta = np.vstack([Theta_group, np.tile(Theta_ref.T, (NS, 1))])  # (NS+NSC, NG)

    # sum_1 = np.zeros((NS+NSC, NG))
    # sum_2 = np.zeros((NS+NSC, NG))
    
    # # Old Psi
    # def Psi(k, l):
    #     T0 = 298.15
    #     return np.exp(-(a[k, l] - a[l, l] + (b[k, l] - b[l, l]) * (T2 - T0)) / T2)
    # # old Psi end

    # # Loop over groups
    # for k in range(NG):
    #     sum_3 = np.zeros(NS+NSC)
    #     prod_1 = np.zeros((NS+NSC, NG))
    #     prod_2 = np.zeros((NS+NSC, NG))

    #     for p in range(NG):
    #         # Term 1
    #         prod_1[:, p] = Theta[:, k] * Psi(k, p)

    #         # Term 2
    #         sum_3 += Theta[:, p] * Psi(p, k)
    #         prod_2[:, p] = Theta[:, k] * Psi(p, k)

    #     sum_1 += prod_1
    #     sum_2 += prod_2 / sum_3[:, None]
    

    # lnGamma = np.tile(Q, (NS+NSC, 1)) * (1 - np.log(sum_1) - sum_2)
    # lnGamma[np.isinf(lnGamma)] = 0
    # lnGamma[np.isnan(lnGamma)] = 0

    # lngammaRes = np.zeros((NS, NC))
    # for k in range(NG):
    #     lnGammaGroup = np.tile(lnGamma[0:NS, k].reshape(-1, 1), (1, NC))
    #     lnGammaPure = lnGamma[NS:NS+NSC, k].reshape(NS, NC) #lnGamma[NS:NS+NSC, k].reshape(NC, NS).T
    #     lngammaRes += np.tile(nu[k, :], (NS, 1)) * (lnGammaGroup - lnGammaPure)

    # # --- Total activity coefficient ---
    # gamma = np.exp(lngamma_comb + lngammaRes)
    