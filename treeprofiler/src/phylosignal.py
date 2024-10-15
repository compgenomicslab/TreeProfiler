#!/usr/bin/env python3
import os, math, re
from multiprocessing.pool import ThreadPool
import numpy as np
from scipy.stats import entropy
import math
#from numba import njit, float64, int64

from pastml.tree import read_tree, name_tree, annotate_dates, DATE, read_forest, DATE_CI, resolve_trees, IS_POLYTOMY, \
    unresolve_trees
from pastml.acr import acr, _serialize_acr
from pastml.annotation import preannotate_forest
from pastml import col_name2cat
from collections import defaultdict, Counter

from treeprofiler.src.utils import add_suffix
from treeprofiler.src.acr_continuous import ml_acr, by_acr

''' ADDITIONAL INFORMATION

[1] Borges, R. et al. (2019). Measuring phylogenetic signal between categorical traits and phylogenies. Bioinformatics, 35, 1862-1869.

[2] Ribeiro, D. et al. (2023). Testing phylogenetic signal with categorical traits and tree uncertainty, Bioinformatics, Volume 39, Issue 7, July 2023

[3] Ishikawa SA, Zhukova A, Iwasaki W, Gascuel O. (2019). A Fast Likelihood Method to Reconstruct and Visualize Ancestral Scenarios. Molecular Biology and Evolution, msz131.
'''

# Source Script of delta method from https://github.com/diogo-s-ribeiro/delta-statistic/blob/master/Delta-Python/delta_functs.py
def mhalpha(a, b, x, l0, se):
    """
    Metropolis-Hastings step for alpha parameter (no numba).
    
    Parameters:
    - a: The current value of the alpha parameter.
    - b: The current value of the beta parameter.
    - x: An array of data points used in the acceptance ratio computations.
    - l0: A constant value used in the acceptance ratio computations.
    - se: The standard deviation for the random walk in the Metropolis-Hastings algorithm.
    
    Returns:
    - The updated alpha value.
    """
    a1 = np.exp(np.random.normal(np.log(a), se))
    lp_a = np.exp(len(x) * (math.lgamma(a1 + b) - math.lgamma(a1)) - a1 * (l0 - np.sum(np.log(x))) - 
                  (len(x) * (math.lgamma(a + b) - math.lgamma(a)) - a * (l0 - np.sum(np.log(x)))))
    
    r = min(1, lp_a)

    while np.isnan(lp_a):
        a1 = np.exp(np.random.normal(np.log(a), se))
        lp_a = np.exp(len(x) * (math.lgamma(a1 + b) - math.lgamma(a1)) - a1 * (l0 - np.sum(np.log(x))) - 
                      (len(x) * (math.lgamma(a + b) - math.lgamma(a)) - a * (l0 - np.sum(np.log(x)))))
        r = min(1, lp_a)

    # Acceptance
    if np.random.uniform(0, 1) < r:
        return a1
    else:
        return a


# Metropolis-Hastings step for beta parameter

def mhbeta(a, b, x, l0, se):
    """
    Metropolis-Hastings step for beta parameter (no numba).
    
    Parameters:
    - a: The current value of the alpha parameter.
    - b: The current value of the beta parameter.
    - x: An array of data points used in the acceptance ratio computations.
    - l0: A constant value used in the acceptance ratio computations.
    - se: The standard deviation for the random walk in the Metropolis-Hastings algorithm.
    
    Returns:
    - The updated beta value.
    """
    b1 = np.exp(np.random.normal(np.log(b), se))
    lp_b = np.exp(len(x) * (math.lgamma(a + b1) - math.lgamma(b1)) - b1 * (l0 - np.sum(np.log(1 - x))) - 
                  (len(x) * (math.lgamma(a + b) - math.lgamma(b)) - b * (l0 - np.sum(np.log(1 - x)))))

    r = min(1, lp_b)

    while np.isnan(lp_b):
        b1 = np.exp(np.random.normal(np.log(b), se))
        lp_b = np.exp(len(x) * (math.lgamma(a + b1) - math.lgamma(b1)) - b1 * (l0 - np.sum(np.log(1 - x))) - 
                      (len(x) * (math.lgamma(a + b) - math.lgamma(b)) - b * (l0 - np.sum(np.log(1 - x)))))
        r = min(1, lp_b)

    # Acceptance
    if np.random.uniform(0, 1) < r:
        return b1
    else:
        return b


# Metropolis-Hastings algorithm using alpha and beta
def emcmc(params):
    """
    Metropolis-Hastings algorithm for alpha and beta parameters (no numba).
    
    Parameters:
    - params: tuple (alpha, beta, x, l0, se, sim, thin, burn)
    
    Returns:
    - Gibbs samples for alpha and beta.
    """
    alpha, beta, x, l0, se, sim, thin, burn = params
    n_size = np.linspace(burn, sim, int((sim - burn) / thin + 1))
    usim = np.round(n_size, 0).astype(int)
    gibbs = []
    p = 0

    for i in range(sim + 1):
        alpha = mhalpha(alpha, beta, x, l0, se)
        beta = mhbeta(alpha, beta, x, l0, se)
        
        if i == usim[p]:
            gibbs.append((alpha, beta))
            p += 1
            
    return np.asarray(gibbs)

# def parallel_emcmc(threads, alpha, beta, x, l0, se, sim, thin, burn):
#     params = [(alpha, beta, x, l0, se, sim, thin, burn) for _ in range(threads)]
#     with ThreadPool(processes=threads-1) as pool:
#         results = pool.map(func=emcmc, iterable=params)
#     return results

# Calculate uncertainty using different types
def entropy_type(prob, ent_type):
    '''prob  = A matrix of ancestral probabilities.
    ent_type = A string indicating the type of entropy calculation. (options: 'LSE', 'SE', or any other value for Gini impurity).'''
    
    # Linear Shannon Entropy
    if ent_type == 'LSE':
        k    = np.shape(prob)[1]
        prob = np.asarray(np.where(prob<=(1/k), prob, prob/(1-k) - 1/(1-k)))
        tent = np.sum(prob, 1)
        
        # Ensure absolutes
        tent = np.asarray(np.where(tent != 0, tent, tent + np.random.uniform(0,1,1)/10000))
        tent = np.asarray(np.where(tent != 1, tent, tent - np.random.uniform(0,1,1)/10000))
        
        return tent

    # Shannon Entropy
    elif ent_type == 'SE':
        k    = np.shape(prob)[1]
        tent = entropy(prob, base=k, axis=1)
        
        # Ensure absolutes
        tent = np.asarray(np.where(tent != 0, tent, tent + np.random.uniform(0,1,1)/10000))
        tent = np.asarray(np.where(tent != 1, tent, tent - np.random.uniform(0,1,1)/10000))

        return tent

    # Ginni Impurity
    else:
        k    = np.shape(prob)[1]
        tent = ((1 - np.sum(prob**2, axis=1))*k)/ (k - 1)
        
        # Ensure absolutes
        tent = np.asarray(np.where(tent != 0, tent, tent + np.random.uniform(0,1,1)/10000))
        tent = np.asarray(np.where(tent != 1, tent, tent - np.random.uniform(0,1,1)/10000))

        return tent


# Calculate delta-statistic after an MCMC step
def delta(x,lambda0,se,sim,thin,burn,ent_type, threads=1):
    '''x     = A matrix of ancestral probabilities.
    lambda0  = A constant value used in the acceptance ratio computations.
    se       = The standard deviation used for the random walk in the Metropolis-Hastings algorithm.
    sim      = The number of total iterations in the Markov Chain Monte Carlo (MCMC) simulation.
    thin     = The thinning parameter, i.e., the number of iterations to discard between saved samples.
    burn     = The number of burn-in iterations to discard at the beginning of the simulation.
    ent_type = A string specifying the type of entropy calculation (options: 'LSE', 'SE', or any other value for Gini impurity).'''
    params = (np.random.exponential(),np.random.exponential(),entropy_type(x, ent_type),lambda0,se,sim,thin,burn)
    # if threads > 1:
    #     mc1    = parallel_emcmc(threads,np.random.exponential(),np.random.exponential(),entropy_type(x, ent_type),lambda0,se,sim,thin,burn)
    #     mc2    = parallel_emcmc(threads,np.random.exponential(),np.random.exponential(),entropy_type(x, ent_type),lambda0,se,sim,thin,burn)
    # else:
    mc1    = emcmc(params)
    mc2    = emcmc(params)
    mchain = np.concatenate((mc1,mc2), axis=0)
    
    deltaA = (np.mean(mchain[:,1]))/(np.mean(mchain[:,0]))
    
    return deltaA

# Calculate the marginal probabilities for each discrete trait
def run_acr_discrete(tree, columns, prediction_method="MPPA", model="F81", threads=1, outdir="./"):
    prop2acr = {}
    column2states = {c: np.array(sorted(list(set(states)))) for c, states in columns.items()}
    features = list(column2states.keys())
    forest = [tree]

    for key in column2states.keys():
        single_column = {key:columns[key]}
        single_column2states = {key:column2states[key]}
        # Run ACR
        acr_result = acr(forest=forest, columns=single_column.keys(), column2states=single_column2states, prediction_method=prediction_method, model=model, threads=threads)
        prop2acr[key] = acr_result
        if outdir:
            _serialize_acr((acr_result[0], outdir))
    return prop2acr, forest[0]

# Calculate the marginal probabilities for each continuous trait
def run_acr_continuous(tree, transformed_dict, model="BM", prediction_method="ML", threads=1, outdir="./"):
    acr_results = {}

    # Parameters for OU model
    alpha = 1.0
    theta = 30.0
    sigma = 5.0
    sigma_prior = 10
    sigma_drift = 5.0

    for key, observed_traits in transformed_dict.items():
        # Run ACR
        if prediction_method == 'ML':
            tree, acr_result = ml_acr(tree, key, observed_traits, model=model, sigma=sigma, alpha=alpha, theta=theta)
        elif prediction_method == 'BAYESIAN':
            tree, acr_result = by_acr(tree, key, observed_traits, model=model, sigma_prior=sigma_prior, sigma_drift=sigma_drift, alpha=alpha, theta=theta)
        acr_results[key] = acr_result

    return acr_result, tree

# Calculate delta-statistic of marginal probabilities each discrete trait
def run_delta(acr_results, tree, run_whole_tree=False, ent_type='LSE', lambda0=0.1, se=0.5, sim=10000, burn=100, thin=10, threads=1):
    prop2delta = {}
    prop2marginals = {}
    leafnames = tree.leaf_names()
    node2leaves = tree.get_cached_content(leaves_only=False)
    # extract the marginal probabilities for each discrete trait
    if run_whole_tree:
        for node in tree.traverse():
            if node.is_leaf:
                continue

            for prop, acr_result in acr_results.items():
                # Get the marginal probabilities for each node
                children_data = acr_result[0]['marginal_probabilities'].loc[[child.name for child in node2leaves[node] if not child.is_leaf]]
                
                #internal_nodes_data[node.name] = children_data
                if node.is_root:
                    marginal_probs = np.asarray(acr_result[0]['marginal_probabilities'].drop(leafnames))
                else:
                    marginal_probs = np.asarray(children_data)
                # run delta for each discrete trait
                # load annotations to leaves
                delta_result = delta(marginal_probs, lambda0, se, sim, burn, thin, ent_type, threads)
                node.add_prop(add_suffix(prop, "delta"), delta_result)
    else:
        # this is the case when we only want to calculate delta for the root
        for prop, acr_result in acr_results.items():
            # Get the marginal probabilities for each node
            marginal_probs = np.asarray(acr_result[0]['marginal_probabilities'].drop(leafnames))
            # run delta for each discrete trait
            # load annotations to leaves
            delta_result = delta(marginal_probs, lambda0, se, sim, burn, thin, ent_type, threads)
            #tree.add_prop(add_suffix(prop, "delta"), delta_result)
            prop2delta[prop] = delta_result
        return prop2delta

# Calculate Pagel's lambda statistic for each continuous trait
def run_lambda():
    return

# Calculate Blomberg's kappa statistic for each continuous trait
def run_kappa():
    return