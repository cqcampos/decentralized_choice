import os
import numpy as np
import pandas as pd
import torch
from scipy.optimize import minimize  # For L-BFGS-B
from torch.quasirandom import SobolEngine
from scipy.stats import qmc


def process_dta(ROOT, K, R, lmbda, homog, eta_out):

    """
     Args:
         ROOT (string): Path where the subfolder decentralized_choice is located
         Other args   : Arguments related to MLE that'll be loaded into the output dictionary

     Returns:
         dictionary: Contains matrices of covariate and admission probabilities for the MLE
     """
    root_norm = os.path.normpath(ROOT)
    data_path = os.path.join(root_norm, "data")
    raw_data = pd.read_stata(os.path.join(data_path, "structural_data_2004_2008_sample.dta"))
    N = raw_data.shape[0]

    # Create indicator matrix of applications and enrollments
    E = raw_data.filter(regex='^D_').to_numpy(dtype=float)
    A = raw_data.filter(regex='^A_').to_numpy(dtype=float)
    J = A.shape[1]
    Z = raw_data.filter(regex='^Z_').to_numpy(dtype=float)
    d = raw_data.filter(regex='^d_').to_numpy(dtype=float)

    # Helper to find column index of max value per row, or -9999 for all-zero rows
    def find_col_index(mat):
        row_sums = mat.sum(axis=1)
        result = np.argmax(mat, axis=1)
        result[row_sums == 0] = -9999
        return result

    choice_col = find_col_index(A)
    E_col = find_col_index(E)

    # Extract covariates as a numpy array
    x_names = [
        'female', 'black', 'hispanic', 'white', 'poverty', 'el', 'eng_home',
        'lag_ela', 'lag_math', 'median_income', 'current_in_mag', 'current_peer_quality'
    ]

    # Extract as NumPy array and column-stack
    X = np.column_stack([raw_data[name].to_numpy(dtype=float) for name in x_names])

    # Demean X
    X -= X.mean(axis=0)
    Q_X = X.shape[1]

    # Add an intercept term
    intercept = np.ones((X.shape[0], 1))
    X_intercept = np.hstack((intercept, X))
    Q_W = X_intercept.shape[1]

    # Process distance matrix
    Q_D =  Q_W + 1
    D_X = np.zeros((N, Q_D, J))

    for j in range(J):
        # Broadcast d[:, j] to (N, Q_D) and multiply by X
        D_X[:, :Q_W, j] = X_intercept * d[:, j].reshape(-1, 1)
        # Last column is squared distance
        D_X[:, Q_W, j] = d[:, j]**2

    # Get admission probabilities and the admission matrix:
    pi_i_df = pd.read_stata(os.path.join(data_path, "student_p_i_2004_2008_sample.dta"))
    pi_i = pi_i_df.iloc[:, 2:].to_numpy(dtype=float)

    # Y matrix (test scores)
    Y_m_6 = raw_data['F1math'].to_numpy(dtype=float)
    Y_e_6 = raw_data['F1ela'].to_numpy(dtype=float)
    Y_m_7 = raw_data['F2math'].to_numpy(dtype=float)
    Y_e_7 = raw_data['F2ela'].to_numpy(dtype=float)
    Y_m_8 = raw_data['F3math'].to_numpy(dtype=float)
    Y_e_8 = raw_data['F3ela'].to_numpy(dtype=float)
    Y_m_6_8 = raw_data['avg_Fmath'].to_numpy(dtype=float)
    Y_e_6_8 = raw_data['avg_Fela'].to_numpy(dtype=float)

    Y = np.column_stack([
        Y_m_6, Y_e_6, Y_m_7, Y_e_7, Y_m_8, Y_e_8, Y_m_6_8, Y_e_6_8
    ])

    # Cohort dummies: yearapp == 2004 ~ 2008
    yearapp    = raw_data['endyear'].to_numpy(dtype=float)

    cohort = np.column_stack([
        (yearapp == 2004).astype(int),
        (yearapp == 2005).astype(int),
        (yearapp == 2006).astype(int),
        (yearapp == 2007).astype(int),
        (yearapp == 2008).astype(int)
    ])

    lat = raw_data['block_lat'].to_numpy(dtype=float)
    lon = raw_data['block_lon'].to_numpy(dtype=float)

    return {
        'N': N, 'J': J,                   # Number of students and schools
        'A': A, 'choice_idx': choice_col, # Application dummy and column idx of applied school (N X J)
        'pi_i' : pi_i,                    # Applicaiton probabilities (N X J)
        'E': E, 'enroll_idx': E_col,      # Enrollment dummy and column idx of enrolled school (N X J)
        'Z': Z,                           # Offer dummy matrix
        'X': X, 'W':X_intercept, 'D': D_X,         # Outcome, covariates, and X interacted with dist
        'Q_X': Q_X, 'Q_D': Q_D, 'Q_W' :Q_W,
        'x_names': x_names,
        'K' : K, 'R' : R, 'homog' : homog, 'lmbda' : lmbda, 'eta_out' : eta_out
    }


def convert_to_tensor(data_dict, device='cpu'):
    """
    Convert numpy arrays in data_dict to PyTorch tensors.

    Args:
        data_dict (dict): Output from process_dta(), including arrays and scalars.
        device (str or torch.device): Device for tensors ('cpu' or 'cuda')
    Returns:
        dict: Same keys as data_dict, where all numpy.ndarray values are converted to
              torch.Tensor on the specified device, with appropriate dtype.
              Other Python types (int, float, list, str) are left unchanged.
    """
    torch_dict = {}
    for key, value in data_dict.items():
        if isinstance(value, np.ndarray):
            # Choose dtype: integers → long, floats → float32
            if np.issubdtype(value.dtype, np.integer):
                dtype = torch.long
            else:
                dtype = torch.float32
            tensor = torch.tensor(value, dtype=dtype, device=device)
            torch_dict[key] = tensor
        else:
            # Leave scalars, lists (e.g., x_names), strings, etc., as-is
            torch_dict[key] = value
    return torch_dict


# Draws from normal for eta/theta
def draw_normal(N, R, J, method='mc', device='cpu', dtype=torch.float32):
    """
    Generate Monte Carlo, Sobol or (scrambled) lattice quasi‑random normal draws.

    Args:
        N (int):    number of “individuals”
        R (int):    draws per individual
        J (int):    dimension of each draw
        method (str):
            'mc'      → torch.randn
            'sobol'   → SobolEngine → inverse‑CDF
            'lattice' → rank‑1 lattice (via SciPy)
        device:     torch device
        dtype:      torch dtype
    Returns:
        Tensor (N, R, J) of ~N(0,1) draws.
    """
    total = N * R

    if method == 'mc':
        return torch.randn(N, R, J, dtype=dtype, device=device)

    elif method == 'sobol':
        engine = SobolEngine(dimension=J, scramble=True, seed=None)
        u = engine.draw(total).to(device=device, dtype=dtype)        # (total, J)
        normals = torch.erfinv(2 * u - 1) * (2 ** 0.5)
        return normals.view(N, R, J)

    elif method == 'lattice':
        # build a Korobov lattice rule
        engine = qmc.Lattice(d=J, seed=None)
        engine.randomize()
        u = engine.random(n=total)                                     # (total, J) in numpy
        u = torch.from_numpy(u).to(device=device, dtype=dtype)
        normals = torch.erfinv(2 * u - 1) * (2 ** 0.5)
        return normals.view(N, R, J)

    else:
        raise ValueError(f"Unknown method '{method}'. Use 'mc','sobol' or 'lattice'.")


# Set the starting values for MLE optimization
def set_start_params(mle_input_dict, perturbation = False, seed_value = 123456):
    """
    Args:
        mle_input_dict (dict): Output from process_dta(), containing data and model matrices for MLE.
        perturbation (bool): If True, apply small random noise to the starting values. Useful for testing the
                             sensitivity of the MLE solution to initial parameter values.
        seed_value (int, optional): Random seed used when perturbation is enabled. Defaults to 123456.
    Returns:
        omega_names (list) : A list of string containing variable labels corresponding to each omega
        omega       (list) : A list of starting value floats
    """
    # Create variable labels
    omega_names = []

    if mle_input_dict["homog"]:
        mean_util = np.array([-1])
        omega_names.append("Mean Utility")
    else:
        # Set j-specific mean based on share of students applying to j
        mean_util = 10*mle_input_dict["A"].mean(axis=0) - 1
        J = mle_input_dict["A"].shape[1]
        omega_names.extend([f"Mean Utility {j+1}" for j in range(J)])

    # Demeaned covariate effects on utility
    beta_x = np.zeros(mle_input_dict["Q_X"])
    omega_names.extend([f"U x {name}" for name in mle_input_dict["x_names"]])

    # Distance effects
    psi_d = np.zeros(mle_input_dict["Q_D"])
    # Set the starting value of baseline distance effects to -1
    psi_d[0] = -1

    omega_names.append("Distance Cost")  # First element
    omega_names.extend([f"Distance Cost x {name}" for name in mle_input_dict["x_names"]])
    omega_names.append("Distance Cost Sq")

    # Mean effects inside exp() cost function
    mu_x_cost = np.zeros(mle_input_dict["Q_W"])

    omega_names.append("App. Cost")
    omega_names.extend([f"App. Cost {name}" for name in mle_input_dict["x_names"]])

    # ln(sd(eta))
    # Exponentiated during the likelihood calculations
    # This ensures that the sd(theta) > 0
    sd_eta = np.array([-0.5])

    # Parameters related to theta
    # In case when K<2, we only need to estimate ln(sd(theta))
    K = mle_input_dict["K"]

    # Always include ln(sd_theta)
    ln_sd_theta = np.zeros(K)

    if K >= 2:
        # Components of soft-max for type shares
        alpha_theta = np.zeros(K - 1)

        # Mean utility value for K-1 types
        mu_theta    = np.zeros(K - 1)
        theta_parms = np.concatenate([ln_sd_theta, alpha_theta, mu_theta])

        omega_names.extend([f"log(sd theta {k+1})" for k in range(K)])
        omega_names.extend([f"alpha theta {k+1}" for k in range(K - 1)])
        omega_names.extend([f"mu theta {k+1}" for k in range(K - 1)])
    else:
        theta_parms = ln_sd_theta

    omega = np.concatenate([mean_util, beta_x, psi_d, mu_x_cost, sd_eta, theta_parms])

    if (perturbation):
        np.random.seed(seed_value)
        omega += np.random.uniform(low=-1, high=1, size=omega.shape)

    return omega, omega_names
