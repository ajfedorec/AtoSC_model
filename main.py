import numpy as np
from model import *
from scipy.integrate import odeint
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from SALib.sample import morris
from SALib.analyze import morris as morris_analysis
import multiprocessing
from joblib import Parallel, delayed
from tqdm import tqdm
import math
import random


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
                       'k_ph': 0.001,  # Dephosphorylation of AtoC
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

    elif model_name == 'original_simple':
        model_params = ['k_b1', 'k_d1', 'k_b2', 'k_d2', 'k_b3', 'k_d3', 'k_ap', 'k_ad', 'k_pt', 'k_ph', 'k_deg',
                        'k_bnd', 'k_unbnd', 'k_pmgexp', 'k_mgbnd', 'k_gexp', 'k_mat', 'k_mgdeg', 'k_lgexp']


        model_species = ['C', 'S', 'Cp', 'Sp', 'pato', 'CS', 'CSp', 'CpS', 'patoCp', 'GFP',
                         'mGFP', 'R', 'mGFPR', 'uGFP']

        model_func = AtoSC_ori_s_dot

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


def leaky_hill(x, ymin, ymax, K, n):
    return ymin + (ymax - ymin) * x ** n / (K ** n + x ** n)


def fit_hill(y, x):
    # try:
    pars, cov = curve_fit(f=leaky_hill,
                          xdata=x,
                          ydata=y,
                          p0=[min(y), max(y), 1E-3, 1.],
                          bounds=(-np.inf, np.inf))
    # except:
    #     pars = np.empty(5)
    #     pars[:] = np.nan
    #     return pars

    return np.append(pars, pars[1] / pars[0])


def morris_problem(eval_params, var_dict):
    bounds = [[var_dict[x] / 10, var_dict[x] * 10] for x in eval_params]

    return {'num_vars': len(eval_params),
            'names': eval_params,
            'bounds': bounds}


def single_dose(model, param_dict, init_dict, k_ap):
    # set parameters and initial conditions for this experiment
    temp_param_dict = param_dict.copy()
    temp_init_dict = init_dict.copy()
    temp_param_dict['k_ap'] = k_ap
    params = list(temp_param_dict.values())
    s0 = list(temp_init_dict.values())

    tmax = 10 * 60 * 60  # Should reach steady state --> can be increased to be more sure but will take longer.
    t_obs = np.linspace(0, tmax, 1001)

    atoSC_obs = odeint(model, s0, t_obs, args=(params,))

    obs_dict = dict(zip(temp_init_dict.keys(), atoSC_obs[-1, :]))
    out_dict = {**obs_dict, **temp_param_dict}  # concatenate parameter and species outputs
    return [out_dict["GFP"], k_ap]


def parallel_dose_response(model, param_dict, init_dict, k_ap_space):
    num_cores = multiprocessing.cpu_count()

    # inputs = tqdm(k_ap_space)

    # store dose response data
    dose_response_data = Parallel(n_jobs=num_cores)(delayed(single_dose)(model, param_dict, init_dict, k_ap) for k_ap in k_ap_space)

    return np.array(dose_response_data)


def sim_dose_response(model, param_dict, init_dict, k_ap_space):
    # store dose response data
    dose_response_data = []

    for k_ap in k_ap_space:
        # set parameters and initial conditions for this experiment
        param_dict['k_ap'] = k_ap
        params = list(param_dict.values())
        s0 = list(init_dict.values())

        # simulate model
        tmax = 10 * 60 * 60  # Should reach steady state --> can be increased to be more sure but will take longer.
        t_obs = np.linspace(0, tmax, 1001)

        atoSC_obs = odeint(model, s0, t_obs, args=(params,))

        obs_dict = dict(zip(init_dict.keys(), atoSC_obs[-1, :]))
        out_dict = {**obs_dict, **param_dict}  # concatenate parameter and species outputs

        out_series = pd.Series(out_dict)

        # update dose response curve
        dose_response_data.append(out_series["GFP"])

    return dose_response_data


def run_morris():
    problem = morris_problem(eval_params, {**param_dict, **init_dict})

    if(1):
        # create samples in morris problem space
        morris_param_values = morris.sample(problem, N=1000, optimal_trajectories=50, seed=123)
        np.savetxt('morris_sample.csv', morris_param_values, delimiter=",")

        # create morris trajectories
        summary_output = []
        # k_ap_space = np.logspace(-8, 1, 64)  # define parameter search space
        for i, X in enumerate(tqdm(morris_param_values)):
            # print("%i of %i" % (i, len(morris_param_values)))
            morris_dict = dict(zip(eval_params, X))
            temp_param_dict = param_dict.copy()
            temp_init_dict = init_dict.copy()
            for j, Y in enumerate(morris_dict):
                if(Y in temp_param_dict):
                    temp_param_dict[Y] = morris_dict[Y]
                elif(Y in temp_init_dict):
                    temp_init_dict[Y] = morris_dict[Y]

            sim = parallel_dose_response(model, temp_param_dict, temp_init_dict, k_ap_space)
            summary_output.append(fit_hill(sim[:, 0], sim[:, 1]))

        hill_params = ['ymin', 'ymax', 'K', 'n', 'fold_change']
        df = pd.DataFrame(summary_output, columns=hill_params)
        df.to_csv('hill_summary.csv')
    # morris_param_values = np.genfromtxt('morris_problem.csv', delimiter=",")
    # hill_params = ['ymin', 'ymax', 'K', 'n', 'fold_change']
    # summary_output = np.genfromtxt('hill_summary.csv', delimiter=",", skip_header=1, usecols=(1,2,3,4,5))

    for i, p in enumerate(hill_params):
        # crappy fix for nan values from Hill model
        Y = np.array(summary_output)[:, i]
        meanY = Y[~np.isnan(Y)].mean()
        Y[np.isnan(Y)] = meanY

        Si = morris_analysis.analyze(problem,
                                     X=morris_param_values,
                                     Y=Y)
        Si.to_df().to_csv('%s_morris.csv' % p)


def sim_weird_morris():
    hill_summary = np.genfromtxt('hill_summary.csv', delimiter=",", skip_header=1, usecols=(1, 2, 3, 4, 5))

    weird_indices = []

    for i in range(len(hill_summary)):
        # find ymin >= ymax
        if hill_summary[i, 0] >= hill_summary[i, 1]:
            weird_indices.append(i)
            break

        # find negative params
        for j in range(5):
            if hill_summary[i, j] < 0:
                weird_indices.append(i)
                break

    # read in samples in morris problem space
    morris_param_values = np.genfromtxt('morris_sample.csv', delimiter=",")

    weird_morris_param_values = morris_param_values[weird_indices, :]

    # create morris trajectories
    summary_output = []
    for i, X in enumerate(weird_morris_param_values):
        morris_dict = dict(zip(eval_params, X))
        temp_param_dict = param_dict.copy()
        temp_init_dict = init_dict.copy()
        for j, Y in enumerate(morris_dict):
            if (Y in temp_param_dict):
                temp_param_dict[Y] = morris_dict[Y]
            elif (Y in temp_init_dict):
                temp_init_dict[Y] = morris_dict[Y]

        sim = sim_dose_response(model, temp_param_dict, temp_init_dict, k_ap_space)
        hill = leaky_hill(k_ap_space,
                          hill_summary[weird_indices[i], 0],
                          hill_summary[weird_indices[i], 1],
                          hill_summary[weird_indices[i], 2],
                          hill_summary[weird_indices[i], 3])

        fig, ax = plt.subplots()
        ax.scatter(k_ap_space, sim)
        ax.plot(k_ap_space, hill)
        ax.set_xscale('log')
        fig.show()


def sim_rand_morris():
    hill_summary = np.genfromtxt('hill_summary.csv', delimiter=",", skip_header=1, usecols=(1, 2, 3, 4, 5))

    weird_indices = random.sample(range(len(hill_summary)), k=5)

    # read in samples in morris problem space
    morris_param_values = np.genfromtxt('morris_sample.csv', delimiter=",")

    weird_morris_param_values = morris_param_values[weird_indices, :]

    # create morris trajectories

    summary_output = []
    k_ap_space = np.logspace(-8, 1, 64)  # define parameter search space
    for i, X in enumerate(weird_morris_param_values):
        morris_dict = dict(zip(eval_params, X))
        temp_param_dict = param_dict.copy()
        temp_init_dict = init_dict.copy()
        for j, Y in enumerate(morris_dict):
            if (Y in temp_param_dict):
                temp_param_dict[Y] = morris_dict[Y]
            elif (Y in temp_init_dict):
                temp_init_dict[Y] = morris_dict[Y]

        sim = sim_dose_response(model, temp_param_dict, temp_init_dict, k_ap_space)
        hill = leaky_hill(k_ap_space,
                          hill_summary[weird_indices[i], 0],
                          hill_summary[weird_indices[i], 1],
                          hill_summary[weird_indices[i], 2],
                          hill_summary[weird_indices[i], 3])

        fig, ax = plt.subplots()
        ax.scatter(k_ap_space, sim)
        ax.plot(k_ap_space, hill)
        ax.set_xscale('log')
        fig.show()


def sim_best_morris():
    hill_summary = np.genfromtxt('hill_summary.csv', delimiter=",", skip_header=1, usecols=(1, 2, 3, 4, 5))

    weird_indices = []

    for i in range(len(hill_summary)):
        # find best 10 fold_change
        if hill_summary[i, 0] >= hill_summary[i, 1]:
            weird_indices.append(i)
            break

        # find negative params
        for j in range(5):
            if hill_summary[i, j] < 0:
                weird_indices.append(i)
                break
    # find best 10 fold_change

    # find lowest 10 K


if __name__ == '__main__':
    # choose ode_model (original, original_dimer, simple)
    model_names = ['original_simple']
    # for MODEL in model_names:
    MODEL = 'original_simple'
    model, param_dict, init_dict = select_model(MODEL)

    # define morris problem space
    eval_params = ['k_b1', 'k_d1', 'k_b2', 'k_d2', 'k_b3', 'k_d3', 'k_ad', 'k_pt', 'k_ph',
                   'k_bnd', 'k_unbnd', 'k_pmgexp', 'k_lgexp',
                   'C', 'S', 'pato']

    k_ap_space = np.logspace(-9, 1, 64)  # define parameter search space
    # run_morris()
    # sim_weird_morris()
    sim_rand_morris()



    # ratio_space = np.linspace(0, 1, 5)
    # C0_space = np.logspace(-1, 1, 5)
    # S0_space = np.logspace(-1, 1, 5)
    # pato0_space = np.linspace(1, 9, 3)
    #
    # # get ready to store output data
    # all_timeseries_data = pd.DataFrame()
    #
    # for C0 in C0_space:
    #     # for ratio in ratio_space:
    #     #     S0 = C0 * ratio
    #     for S0 in S0_space:
    #         for pato0 in pato0_space:
    #
    #             # store dose response data
    #             dose_response_data = []
    #
    #             for k_ap in k_ap_space:
    #                 # set parameters and initial conditions for this experiment
    #                 temp_param_dict = param_dict.copy()
    #                 temp_param_dict['k_ap'] = k_ap
    #                 temp_init_dict = init_dict.copy()
    #                 temp_init_dict['C'] = C0
    #                 temp_init_dict['S'] = S0
    #                 temp_init_dict['pato'] = pato0
    #                 params = list(temp_param_dict.values())
    #                 s0 = list(temp_init_dict.values())
    #
    #                 # simulate model
    #                 out_series = simulate_model(model, s0, params)
    #
    #                 # update dose response curve
    #                 dose_response_data.append(out_series["GFP"])
    #
    #                 all_timeseries_data = all_timeseries_data.append(out_series, ignore_index=True)
    #
    #             # fit hill function to dose response curve
    #             pars = fit_hill(dose_response_data, k_ap_space)
    #
    #             print(pars)
    #
    # # # save results to .csv
    # # all_data.to_csv('data_' + MODEL + '.csv')
    #
    # # PLOT
    # p = sns.FacetGrid(data=all_timeseries_data,
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
