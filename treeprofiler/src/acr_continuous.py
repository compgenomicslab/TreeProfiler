# methods.py
import numpy as np
import pymc as pm
import arviz as az
import aesara.tensor as at

def build_variance_covariance_matrix(tree, species, sigma, alpha=None, model='BM'):
    """
    Build the variance-covariance matrix for BM or OU models.
    
    Parameters:
    - tree: Phylogenetic tree
    - species: List of species names
    - sigma: Drift rate (variance scaling)
    - alpha: Selection strength (OU model only)
    - model: 'BM' or 'OU' to select the model type
    
    Returns:
    - Variance-covariance matrix (V)
    """
    n = len(species)
    V = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            leaf_i = next(tree.search_leaves_by_name(name=species[i]))
            leaf_j = next(tree.search_leaves_by_name(name=species[j]))
            mrca = tree.common_ancestor([leaf_i, leaf_j])
            shared_time = tree.get_distance(tree, mrca)
            
            if model == 'BM':
                V[i, j] = V[j, i] = sigma ** 2 * shared_time
            elif model == 'OU':
                V[i, j] = V[j, i] = (sigma ** 2 / (2 * alpha)) * (1 - np.exp(-2 * alpha * shared_time))
    return V

def bm_model(V, Y, sigma):
    """
    Likelihood function for Brownian Motion model.
    
    Parameters:
    - V: Variance-covariance matrix
    - Y: Observed trait values vector
    - sigma: Drift rate
    
    Returns:
    - Likelihood of the observed data
    """
    # Placeholder function structure; for Bayesian, the actual likelihood computation happens in PyMC3
    return V, Y

def ou_model(V, Y, sigma, alpha, theta):
    """
    Likelihood function for Ornstein-Uhlenbeck model.
    
    Parameters:
    - V: Variance-covariance matrix
    - Y: Observed trait values vector
    - sigma: Drift rate
    - alpha: Selection strength
    - theta: Optimal trait value
    
    Returns:
    - Likelihood of the observed data
    """
    # Placeholder function structure; for Bayesian, the actual likelihood computation happens in PyMC3
    return V, Y, alpha, theta

def ml_acr(tree, prop, observed_traits, model='BM', sigma=1.0, alpha=None, theta=None):
    """
    Maximum Likelihood Ancestral Character Reconstruction.
    
    Parameters:
    - tree: Phylogenetic tree
    - observed_traits: Observed trait values
    - model: 'BM' or 'OU'
    - sigma: Drift rate
    - alpha: Selection strength (OU model only)
    - theta: Optimal trait value (OU model only)
    
    Returns:
    - Annotated tree with estimated traits
    - Results with node values and credible intervals
    """
    species = list(observed_traits.keys())
    V = build_variance_covariance_matrix(tree, species, sigma, alpha, model)
    Y = np.array([observed_traits[sp] for sp in species])

    # Compute ML estimate for the root value
    V_inv = np.linalg.inv(V)
    ones = np.ones(len(species))
    root_value = (ones @ V_inv @ Y) / (ones @ V_inv @ ones)
    
    results = {'root': {prop: root_value}}
    print(f"Estimated ancestral {prop} value at the root ({model}-ML): {root_value:.2f}")

    # Recursive estimation for internal nodes
    def estimate_internal_states(node, parent_value):
        if node.is_leaf:
            node.add_prop(prop, observed_traits[node.name])
            results[node.name] = {prop: observed_traits[node.name]}
        else:
            left_child, right_child = node.children
            branch_length = node.dist if node.dist is not None else 1  # Ensure branch length is not None
            
            # Calculate expected value based on the model
            expected_value = parent_value * np.exp(-alpha * branch_length) + theta * (1 - np.exp(-alpha * branch_length)) if model == 'OU' else parent_value
            
            # Estimate values for child nodes with noise based on branch length
            left_child_value = expected_value + np.random.normal(0, np.sqrt(branch_length) if branch_length > 0 else 0)
            right_child_value = expected_value + np.random.normal(0, np.sqrt(branch_length) if branch_length > 0 else 0)

            # Assign the estimated trait to the current node
            node.add_prop(prop, expected_value)
            results[node.name or 'Unnamed'] = {prop: expected_value}

            # Recurse to estimate states for children nodes
            estimate_internal_states(left_child, left_child_value)
            estimate_internal_states(right_child, right_child_value)

    estimate_internal_states(tree, root_value)
    return tree, results



def by_acr(tree, prop, observed_traits, model='BM', sigma_prior=10, sigma_drift=5.0, alpha=None, theta=None):
    """
    Bayesian Inference Ancestral Character Reconstruction using PyMC.
    
    Parameters:
    - tree: Phylogenetic tree
    - observed_traits: Observed trait values
    - model: 'BM' or 'OU'
    - sigma_prior: Prior standard deviation
    - sigma_drift: Drift rate
    - alpha: Selection strength (OU model only)
    - theta: Optimal trait value (OU model only)
    
    Returns:
    - Annotated tree with estimated traits and credible intervals
    - Results with node values and credible intervals
    """
    species = list(observed_traits.keys())
    Y = np.array([observed_traits[sp] for sp in species])
    
    # Build the variance-covariance matrix
    V = build_variance_covariance_matrix(tree, species, sigma_drift, alpha, model)

    # Ensure V is a proper float64 NumPy array
    V = np.array(V, dtype=np.float64)

    with pm.Model() as bayesian_model:
        # Define priors for root_value
        root_value = pm.Normal('root_value', mu=theta if model == 'OU' else 30, sigma=sigma_prior)
        
        # Use Data within the model context
        V_tensor = pm.Data("V", V)

        # Define the likelihood based on the multivariate normal distribution
        pm.MvNormal('Y', mu=root_value, cov=V_tensor, observed=Y)

        # Perform sampling
        trace = pm.sample(2000, tune=1000)
    
    # Extract the root_value samples
    root_samples = trace.posterior['root_value'].values.flatten()

    mean_value = np.mean(root_samples)
    lower_bound = np.percentile(root_samples, 2.5)
    upper_bound = np.percentile(root_samples, 97.5)
    
    results = {'root': {prop: mean_value, 'credible_interval': (lower_bound, upper_bound)}}
    root_node = tree
    root_node.add_prop(prop, mean_value)

    print(f"Root node ({model}-Bayesian): Estimated Trait = {mean_value:.2f}, 95% CI = [{lower_bound:.2f}, {upper_bound:.2f}]")

    # Recursive estimation for internal nodes
    def infer_ancestral_state(node, parent_samples, sigma):
        if node.is_leaf:
            node.add_prop(prop, observed_traits[node.name])
            results[node.name] = {prop: observed_traits[node.name]}
            return np.array([observed_traits[node.name]] * len(parent_samples))
        else:
            branch_length = node.dist
            node_samples = parent_samples + np.random.normal(0, np.sqrt(sigma ** 2 * branch_length), len(parent_samples))
            mean_value = np.mean(node_samples)
            lower_bound = np.percentile(node_samples, 2.5)
            upper_bound = np.percentile(node_samples, 97.5)

            node.add_prop(prop, mean_value)
            results[node.name] = {prop: mean_value, 'credible_interval': (lower_bound, upper_bound)}

            print(f"Internal node {node.name or 'Unnamed'}: Estimated Trait = {mean_value:.2f}, 95% CI = [{lower_bound:.2f}, {upper_bound:.2f}]")

            for child in node.children:
                infer_ancestral_state(child, node_samples, sigma)

    for child in root_node.children:
        infer_ancestral_state(child, root_samples, sigma_drift)

    return tree, results
