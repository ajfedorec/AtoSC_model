import numpy as np
from model import *
from scipy.integrate import odeint
import pandas as pd
import seaborn as sns


def select_model(model_name):
    """
    This function gets the function and inputs to run an ODE model

    :param model_name: string identifying which model to use: 'original', 'original_dimer', 'simple'
    :return: a tuple containing the model ODE function and dictionaries for parameters and initial conditions
    """
    model_func = None
    model_params = None
    model_species = None

    full_param_dict = {'k_b1': 0.5,  # AtoC binds AtoSP
                       'k_d1': 0.5,  # AtoC unbinds AtoSP
                       'k_b2': 0.05,  # AtoCP binds AtoS
                       'k_d2': 0.5,  # AtoCP unbinds AtoS
                       'k_b3': 0.5,  # AtoC binds AtoS
                       'k_d3': 0.5,  # AtoC unbinds AtoS
                       'k_ap': 0.15,  # Acetoacetate phosphorylates AtoS
                       'k_ad': 0.001,  # Dephosphorylation of AtoSP
                       'k_pt': 1.5,  # Phosphorylation of AtoC
                       'k_ph': 0.00,  # Dephosphorylation of AtoC
                       'k_dim': 0.0083,  # Dimerisation of AtoCP, from Merk et al., Biorxiv, 2020
                       'k_mon': 0.5,  # Monomeristaion of AtoCP, from Merk et al., Biorxiv, 2020
                       'k_dbnd': 0.5,  # Promoter binding - phosphorylated dimer
                       'k_dunbnd': 0.05,  # Promoter unbinding - phosphorylated dimer
                       'k_b4': 0.5,  # Binding of alternative dephosphatase
                       'k_d4': 0.5,  # Unbinding of alternative dephosphatase
                       'k_cat': 0.05,  # Alternative dephosphorylation
                       'k_exp': 0.08,  # GFP expression
                       'k_lexp': 0.00008,  # Leaky GFP expression
                       'k_deg': 0.00057,  # GFP degradation
                       'k_bnd': 0.5,  # Promoter binding - phosphorylated monomer
                       'k_unbnd': 0.05,  # Promoter unbinding - phosphorylated monomer
                       'k_nsbnd': 0.01,  # Promoter binding - unphosphorylated monomer
                       'k_nsunbnd': 0.001,  # Promoter unbinding - unphosphorylated monomer
                       'k_pmgexp': 0.07,  # GFP transcription - phosphorylated monomer
                       'k_npmgexp': 0.007,  # GFP transcription - unphosphorylated monomer
                       'k_mgbnd': 0.5,  # Ribosome binding mRNA
                       'k_gexp': 0.08,  # GFP translation
                       'k_mat': 0.00042,  # GFP maturation
                       'k_mgdeg': 0.00223,  # mRNA degradation
                       'k_dmgexp': 0.07,  # GFP transcription - phosphorylated dimer
                       'k_lgexp': 0.00007}  # leaky GFP transcription

    full_init_dict = {'C': 6.0,
                      'S': 0.17,
                      'Cp': 0.00,
                      'Sp': 0.00,
                      'pato': 10.0,
                      'Ph': 0.17,
                      'Cpd': 0.00,
                      'CS': 0.00,
                      'CSp': 0.00,
                      'CpS': 0.00,
                      'patoCpd': 0.00,
                      'PhCp': 0.00,
                      'patoCp': 0.00,
                      'patoC': 0.00,
                      'GFP': 0.00,
                      'mGFP': 0.00,
                      'R': 10.00,
                      'mGFPR': 0.00,
                      'uGFP': 0.00
                      }

    if model_name == 'original':
        model_params = ['k_b1', 'k_d1', 'k_b2', 'k_d2', 'k_b3', 'k_d3', 'k_ap', 'k_ad', 'k_pt', 'k_ph', 'k_b4', 'k_d4',
                        'k_cat', 'k_deg', 'k_bnd', 'k_unbnd', 'k_nsbnd', 'k_nsunbnd', 'k_pmgexp', 'k_npmgexp',
                        'k_mgbnd', 'k_gexp', 'k_mat', 'k_mgdeg']

        model_species = ['C', 'S', 'Cp', 'Sp', 'pato', 'Ph', 'CS', 'CSp', 'CpS', 'PhCp', 'patoCp', 'patoC', 'GFP',
                         'mGFP', 'R', 'mGFPR', 'uGFP']

        model_func = AtoSC_ori_dot

    elif model_name == 'original_dimer':
        model_params = ['k_b1', 'k_d1', 'k_b2', 'k_d2', 'k_b3', 'k_d3', 'k_ap', 'k_ad', 'k_pt', 'k_ph', 'k_dim',
                        'k_mon', 'k_b4', 'k_d4', 'k_cat', 'k_deg', 'k_dbnd', 'k_dunbnd', 'k_dmgexp', 'k_mgbnd',
                        'k_gexp', 'k_mat', 'k_mgdeg', 'k_lgexp']

        model_species = ['C', 'S', 'Cp', 'Sp', 'pato', 'Ph', 'Cpd', 'CS', 'CSp', 'CpS', 'patoCpd', 'PhCp', 'GFP',
                         'mGFP', 'R', 'mGFPR', 'uGFP']

        model_func = AtoSC_ori_d_dot

    elif model_name == 'simple':
        model_params = ['k_b1', 'k_d1', 'k_b2', 'k_d2', 'k_b3', 'k_d3', 'k_ap', 'k_ad', 'k_pt', 'k_ph', 'k_dim',
                        'k_mon', 'k_dbnd', 'k_dunbnd', 'k_b4', 'k_d4', 'k_cat', 'k_exp', 'k_lexp', 'k_deg']

        model_species = ['C', 'S', 'Cp', 'Sp', 'pato', 'Ph', 'Cpd', 'CS', 'CSp', 'CpS', 'patoCpd', 'PhCp', 'GFP']

        model_func = AtoSC_dot

    else:
        exit(-1)

    model_param_dict = {k: full_param_dict[k] for k in model_params if k in full_param_dict}

    model_init_dict = {k: full_init_dict[k] for k in model_species if k in full_init_dict}

    return model_func, model_param_dict, model_init_dict


if __name__ == '__main__':
    # choose ode_model (original, original_dimer, simple)
    model_names = ['original', 'original_dimer', 'simple']
    # for MODEL in model_names:
    MODEL = 'original_dimer'
    model, param_dict, init_dict = select_model(MODEL)

    # define parameter search space
    k_ap_space = np.logspace(-6, 1, 100)
    C0_space = np.logspace(-1, 1, 5)
    S0_space = np.logspace(-1, 1, 5)
    pato0_space = np.linspace(1, 9, 3)

    # get ready to store output data
    C0s = []
    S0s = []
    pato0s = []
    k_aps = []
    GFPs = []
    models = []

    for C0 in C0_space:
        for S0 in S0_space:
            for pato0 in pato0_space:
                for k_ap in k_ap_space:
                    # set parameters and initial conditions for this experiment
                    temp_param_dict = param_dict
                    temp_param_dict['k_ap'] = k_ap
                    temp_init_dict = init_dict
                    temp_init_dict['C'] = C0
                    temp_init_dict['S'] = S0
                    temp_init_dict['pato'] = pato0
                    params = list(temp_param_dict.values())
                    s0 = list(temp_init_dict.values())

                    tmax = 10 * 60 * 60  # Should reach steady state --> can be increased to be more sure but will take longer.
                    t_obs = np.linspace(0, tmax, 1001)

                    atoSC_obs = odeint(model, s0, t_obs, args=(params,))

                    GFPs.append(atoSC_obs[-1, 12])  # Get the final GFP observation. Can we remove the hardcoded position?
                    C0s.append(C0)
                    S0s.append(S0)
                    pato0s.append(pato0)
                    k_aps.append(k_ap)
                    models.append(MODEL)

    # stick data into Pandas dataframe
    all_data = pd.DataFrame({
        "model": models,
        "C0": C0s,
        "S0": S0s,
        "pato0": pato0s,
        "k_ap": k_aps,
        "GFP": GFPs
    })

    # save results to .csv
    all_data.to_csv('data_' + MODEL + '.csv')

    # # PLOT
    # p = sns.FacetGrid(data=all_data,
    #                   col="C0",
    #                   row="S0",
    #                   hue="pato0")
    #
    # p.map(sns.lineplot,
    #       "k_ap",
    #       "GFP")
    #
    # p.set(xscale='log')
    # p.add_legend()
    # p.fig.subplots_adjust(top=0.9)
    # p.fig.suptitle('Model = ' + MODEL)
    # p.fig.show()
