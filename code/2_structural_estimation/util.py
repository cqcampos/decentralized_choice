import os
import numpy as np
import pandas as pd
import torch                         
from scipy.optimize import minimize  # For L-BFGS-B
from torch.quasirandom import SobolEngine # For normal draws
from scipy.stats import qmc  # For normal draws (Lattice)
import psutil
import gc
torch.set_grad_enabled(False) 

def print_cpu_memory():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / (1024 ** 2)  # in MB
    print(f"CPU Memory Usage: {mem:.2f} MB")


# Loads in datasets used for MLE
def process_dta(ROOT, K, R, lmbda, homog, eta_out, in_mag_ind, device, full_samp = False, last_yr=2008):

    """
     Args:
         ROOT (string): Path where the subfolder decentralized_choice is located
         Other args   : Arguments related to MLE that'll be loaded into the output dictionary

     Returns:
         dictionary: Contains matrices of covariate and admission probabilities for the MLE
     """
    root_norm = os.path.normpath(ROOT)
    data_path = os.path.join(root_norm, "data")
    sample_str = "_sample" if not full_samp else ""
    raw_data = pd.read_stata(os.path.join(data_path, f"structural_data_2004_{last_yr}{sample_str}.dta"))
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

    if not in_mag_ind:
        [x for x in x_names if x != "current_in_mag"] 

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
    pi_i_df = pd.read_stata(os.path.join(data_path, f"student_p_i_2004_{last_yr}{sample_str}.dta"))
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

    # figure out the min/max years
    years_range = np.arange(yearapp.min(), yearapp.max() + 1)

    # build dummy columns with list comprehension
    cohort = np.column_stack([(yearapp == y).astype(int) for y in years_range])
    

    lat = raw_data['block_lat'].to_numpy(dtype=float)
    lon = raw_data['block_lon'].to_numpy(dtype=float)

    # Note: choice_col tells which of 40 schools students applied
    #       In contrast, E_col tells which of the 41 options students chose. 
    #       (The 0 is for the outside option for E_col)
    chose_app_sch = (choice_col == (E_col-1)).astype(int)
    has_offer     =  Z.sum(axis=1).astype(int)
    return {
        'N': N, 'J': J,                   # Number of students and schools
        'A': A, 'choice_idx': choice_col, # Application dummy and column idx of applied school (N X J)
        'pi_i' : pi_i,                    # Applicaiton probabilities (N X J)
        'E': E, 'enroll_idx': E_col,      # Enrollment dummy and column idx of enrolled school (N X J)
        'Z': Z,                           # Offer dummy matrix
        'X': X, 'W':X_intercept, 'D': D_X,         # Outcome, covariates, and X interacted with dist
       'Q_X': Q_X, 'Q_D': Q_D, 'Q_W' :Q_W,
        'x_names': x_names,
        'chose_app_sch': chose_app_sch, # Did student enroll into school they applied to?
        'has_offer'    : has_offer,   # Does student has a offer from a magnet school?
        'K' : K, 'R' : R, 'homog' : homog, 'lmbda' : lmbda, 'eta_out' : eta_out, 'in_mag_ind':in_mag_ind, 'device' : device
    }


# Convert MLE inputs which are np.array to PyTorch tensor
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


def convert_to_tensor_blocks(data_dict, n_blocks, sims_theta, sims_c):
    """
    Partition observation-specific numpy arrays in data_dict into n_blocks and assign each to a separate GPU.

    Args:
        data_dict (dict): Output from process_dta(), including arrays and scalars.
        n_blocks (int): Number of partitions (assumes 1 GPU per block)

    Returns:
        list of dicts: Each dict corresponds to a partition with tensors on cuda:b.
    """

    N = data_dict['N']
    block_size = (N + n_blocks - 1) // n_blocks  # ceiling division

    # Keys that vary across individuals (dim[0] == N) and need partitioning
    variable_keys = [
        'A', 'choice_idx', 'pi_i', 'E', 'Z', 'X', 'W', 'D', 
        'chose_app_sch', 'has_offer'
    ]

    # List to store dict for each block
    block_dicts = []

    for b in range(n_blocks):
        start = b * block_size
        end = min((b + 1) * block_size, N)

        if data_dict["device"] == torch.device("cuda") and torch.cuda.device_count() == 1:
            device = torch.device("cuda:0") 
        elif data_dict["device"] == torch.device("cuda") and (b < torch.cuda.device_count()): 
            device = torch.device(f"cuda:{b}")
        else:
            device = torch.device("cpu")
        
        print(f"Loading block {b} (row idx {start}:{end}) to {device}...")

        block_dict = {}
        for key, val in data_dict.items():
            if key in variable_keys:
                sliced_val = val[start:end] if val.ndim == 1 else val[start:end, ...]
                if np.issubdtype(sliced_val.dtype, np.integer):
                    dtype = torch.long
                else:
                    dtype = torch.float32
                block_dict[key] = torch.tensor(sliced_val, dtype=dtype, device=device)
            elif key == 'N':
                block_dict['N'] = end - start  # block-specific N
            else:
                block_dict[key] = val  # shared values (e.g., K, J, device, lmbda)

        block_dict['device'] = device

        #  Slice sims_theta and sims_c for this block
        block_dict["sims_theta"] = sims_theta[start:end].to(device)
        block_dict["sims_c"]     = sims_c[start:end].to(device)
        block_dicts.append(block_dict)

    

    return block_dicts





# Simulate normal draws to later transform into eta and theta
def draw_normal(N, J, R, method='mc', device='cpu', dtype=torch.float32, seed=12345):
    """
    Generate Monte Carlo, Sobol, or (scrambled) lattice quasi‑random normal draws.

    Args:
        N (int):    number of individuals
        J (int):    dimension of each draw (set to None or 1 for scalar draw per individual)
        R (int):    draws per individual
        method (str): 'mc', 'sobol', or 'lattice'. if 'debug', it returns N by J by R array filled with scalars
        device:     torch device
        dtype:      torch dtype
        seed:       seed to use for random draws

    Returns:
        Tensor:
            If J is None or 1: shape (N, R)
            Else: shape (N, J, R)
    """
    torch.manual_seed(seed)

    scalar_draw = J is None or J == 1
    dim = 1 if scalar_draw else J
    total = N * R

    if method == 'mc':
        normals = torch.randn(N, R, dtype=dtype, device=device) if scalar_draw \
            else torch.randn(N, dim, R, dtype=dtype, device=device)

    elif method == 'sobol':
        engine = SobolEngine(dimension=dim, scramble=True, seed=seed)
        u = engine.draw(total).to(device=device, dtype=dtype)          # (total, dim)
        eps = torch.finfo(dtype).eps  # Smallest representable number > 0
        u = u.clamp(min=eps, max=1 - eps)
        print("Min u:", u.min().item())
        print("Max u:", u.max().item())
        normals = torch.erfinv(2 * u - 1) * (2 ** 0.5)
        if scalar_draw:
            normals = normals.view(N, R)
        else:
            normals = normals.view(N, R, dim).transpose(1, 2)          # (N, J, R)

    elif method == 'lattice':
        gens = torch.arange(1, dim+1, dtype=dtype, device=device)
        idx = torch.arange(total, dtype=dtype, device=device).unsqueeze(1)
        u = (idx * gens / total) % 1.0
        shift = torch.rand(1, dim, dtype=dtype, device=device)
        u = (u + shift) % 1.0
        u = u.clamp(min=1e-6, max=1 - 1e-6)  # clamp just in case
        normals = torch.erfinv(2 * u - 1) * (2 ** 0.5)

        if scalar_draw:
            normals = normals.view(N, R)
        else:
            normals = normals.view(N, R, dim).transpose(1, 2)  # (N, J, R)

    elif method == 'debug':
        if scalar_draw:
            # Vary across n and r: use broadcasting of outer sum
            base = torch.arange(N, dtype=dtype, device=device).view(N, 1)
            step = torch.arange(R, dtype=dtype, device=device).view(1, R)
            normals = (torch.log(base+1) + step) * 0.01
        else:
            # Efficient broadcasting without storing large intermediate arrays
            n = torch.arange(N, dtype=dtype, device=device).reshape(N, 1, 1)
            j = torch.arange(J, dtype=dtype, device=device).reshape(1, J, 1)
            r = torch.arange(R, dtype=dtype, device=device).reshape(1, 1, R)
            normals = (torch.log(n+1) + j + r) * 0.01
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'mc', 'sobol' or 'lattice'.")

    return normals


# Set the starting values for MLE optimization
def set_start_params(mle_input_dict, last_yr = 2008, perturbation = False, seed_value = 123456, hot_start=False):
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
        mean_util = 10*mle_input_dict["A"].mean(axis=0)
        J = mle_input_dict["A"].shape[1]
        omega_names.extend([f"Mean Utility {j+1}" for j in range(J)])

    # Demeaned covariate effects on utility
    beta_x = np.zeros(mle_input_dict["Q_X"])
    omega_names.extend([f"U x {name}" for name in mle_input_dict["x_names"]])

    # Distance effects
    psi_d = np.zeros(mle_input_dict["Q_D"])
    # Set the starting value of baseline distance effects to -1
    psi_d[0] = -0.3

    omega_names.append("Distance Cost")  # First element
    omega_names.extend([f"Distance Cost x {name}" for name in mle_input_dict["x_names"]])
    omega_names.append("Distance Cost Sq")

    # Mean effects inside exp() cost function
    mu_x_cost = np.zeros(mle_input_dict["Q_W"])
    mu_x_cost[0] = 0.3

    omega_names.append("App. Cost")
    omega_names.extend([f"App. Cost {name}" for name in mle_input_dict["x_names"]])

    # ln(sd(eta))
    # Exponentiated during the likelihood calculations
    # This ensures that the sd(theta) > 0
    sd_eta = np.array([-0.5])
    omega_names.append("log(sd eta)")

    # Parameters related to theta
    # In case when K<2, we only need to estimate ln(sd(theta))
    K = mle_input_dict["K"]

    # Always include ln(sd_theta)
    ln_sd_theta = np.zeros(K)
    omega_names.extend([f"log(sd theta {k+1})" for k in range(K)])
    
    if K >= 2:
        # Components of soft-max for type shares
        alpha_theta =  np.arange(K - 1, 0, -1)
        # Mean utility value for K-1 types
        mu_theta    = -0.5*np.arange(K - 1, 0, -1)
        theta_parms = np.concatenate([ln_sd_theta, alpha_theta, mu_theta])

        omega_names.extend([f"alpha theta {k+1}" for k in range(K - 1)])
        omega_names.extend([f"mu theta {k+1}" for k in range(K - 1)])
    else:
        theta_parms = ln_sd_theta

    omega = np.concatenate([mean_util, beta_x, psi_d, mu_x_cost, sd_eta, theta_parms])

    if (perturbation):
        np.random.seed(seed_value)
        omega += np.random.uniform(low=-1, high=1, size=omega.shape)

    # If needed, supply starting values manually 
    if hot_start:    
        if K==3 and mle_input_dict["eta_out"]:
            if last_yr ==2013:
                last_omega_0 = np.array([
                    -0.082344, -0.812008, -1.657717, -2.702595, 1.672002, 1.428475, -1.928099, -0.531554,
                    -1.535448, -1.792133, -2.338426, -1.606545, 0.223483, -0.901974, -0.633352, -0.598536,
                    1.634750, -2.540916, -0.963694, -1.282817, 0.230879, 0.861140, -1.913761, -1.024508,
                    -1.102711, 0.682461, -1.704199, -1.273383, -1.352369, -0.040991, -1.561058, 0.095407,
                    -2.344869, -0.751074, -2.311019, -0.977392, -1.581657, 2.189675, -1.302953, -1.507649,
                    -0.402384, 0.044038, 2.522300, -1.180310, -1.631329, -0.945015, -0.984468, -0.544250,
                    -0.554778, -1.059589, -0.362605, 3.226497, 4.136221,
                    0.006077, -0.114490, -0.109517, -1.107454, 0.060185, -0.219201, -0.148606, 0.166927,
                    0.201657, -0.058748, 1.270226, 0.130848,
                    -0.326071, -0.001849, 0.018942, -0.010787, -0.127129, -0.016193, 0.016478, 0.022929,
                    0.000121, -0.006810, -0.001514, 0.075448, 0.026666, 0.003646,
                    0.225312, 0.002595, -0.036312, 0.064333, -0.083029, 0.012232, 0.011003, -0.011522,
                    -0.033674, -0.018631, 0.005244, -0.096653, -0.025149,
                    -0.995783, -1.562508, -2.038194, -1.038199,
                    1.904746, 2.156522,
                    -0.922100, 0.390567
                ])
                

                last_omega_0 = np.array([
                        -1.32217771e-01, -8.70095008e-01, -1.75852438e+00, -2.81430749e+00,
                        1.61847704e+00,  1.20576946e+00, -1.97187257e+00, -5.66734038e-01,
                        -1.57768659e+00, -1.82657872e+00, -2.40363309e+00, -1.60682441e+00,
                        1.84788175e-01, -9.47621889e-01, -6.90232235e-01, -6.43365084e-01,
                        1.58744340e+00, -2.60851116e+00, -1.00006462e+00, -1.31719313e+00,
                        1.84236503e-01,  8.13016871e-01, -1.97331620e+00, -1.07029127e+00,
                        -1.14715816e+00,  6.38674904e-01, -1.80634231e+00, -1.32336576e+00,
                        -1.37252942e+00, -8.42440580e-02, -1.59526859e+00,  4.63683567e-02,
                        -2.41214178e+00, -8.37722269e-01, -2.35689698e+00, -1.02148455e+00,
                        -1.62999020e+00,  2.12030204e+00, -1.34037548e+00, -1.54687694e+00,
                        -4.55410954e-01,  2.23096716e-03,  2.47269254e+00, -1.22447991e+00,
                        -1.67872296e+00, -9.01994684e-01, -1.04378849e+00, -5.85888172e-01,
                        -6.50642696e-01, -1.10803072e+00, -4.12909939e-01,  3.16316865e+00,
                        4.04649520e+00,  1.01825145e-03, -1.05177078e-01, -1.06782020e-01,
                        -1.09837379e+00,  7.05140220e-02, -2.75192258e-01, -1.35925029e-01,
                        1.74469204e-01,  2.01377423e-01, -5.89721464e-02,  1.27402269e+00,
                        1.22983414e-01, -3.24293787e-01, -1.75825718e-03,  1.60303061e-02,
                        -1.22532664e-02, -1.30616810e-01, -1.54285542e-02,  1.72439752e-02,
                        2.36264454e-02, -1.80377574e-04, -6.54024083e-03, -1.44372853e-03,
                        7.49707712e-02,  2.64667581e-02,  3.63884650e-03,  2.08410899e-01,
                        2.61186030e-03, -3.53482592e-02,  6.45838147e-02, -8.31671733e-02,
                        1.34950568e-02,  9.18135977e-03, -1.02168664e-02, -3.29008249e-02,
                        -1.85985142e-02,  5.15625180e-03, -9.73609280e-02, -2.60463578e-02,
                        -1.01376227e+00, -1.66609549e+00, -2.07630755e+00, -6.13482919e-01,
                        1.91335463e+00,  2.00155779e+00, -8.75870673e-01,  4.30343607e-01
                    ])
            elif last_yr == 2008:
                last_omega_0 = np.array([
                -1.3738813582, -1.8131448801,  0.3322918695,  0.1828422365, -2.9219479949,
                -1.3880158025, -2.5299264494, -2.8505220889, -3.3007381730, -1.0030499825,
                -2.1859844829, -1.6317770542, -1.7693628546,  0.4295998814, -3.3995154961,
                -1.8030147856, -2.4889515661, -0.7013622166, -0.3677407217, -2.9851618048,
                -1.9410918911, -1.8070196962, -0.4380884385, -2.4578770378, -1.2841329768,
                -2.6699177098, -1.2459118087, -3.3366222028, -3.2755558878,  0.7652032472,
                -2.1392219084, -2.3187437385, -1.3419951396, -0.8062325298,  1.2990470138,
                -2.1655967849, -2.6838978663, -1.3618228345,  1.5420457954,  7.6592284467,
                -0.0260878015, -0.1145484008, -0.2155931525, -1.1573838908,  0.0990846220,
                -0.3364281073, -0.0769051059,  0.2393047341,  0.1973153443, -0.0646135518,
                1.5339229697,  0.1585593693, -0.2805871747, -0.0054473725,  0.0102690190,
                -0.0106919164, -0.1664611367, -0.0208464096,  0.0032503056,  0.0081994382,
                0.0005325848, -0.0080815157, -0.0008007884,  0.0702064611,  0.0244372586,
                0.0027340396, -0.2793898098, -0.0037414934, -0.0594522496,  0.0642134945,
                -0.1154789557, -0.0005690875, -0.0034363702, -0.0185022971, -0.0369477275,
                -0.0252408417,  0.0075722274, -0.0858324749, -0.0353668282, -1.5643172237,
                -0.4127090570,  0.1992331564, -0.0617436418,  2.4684198426,  0.1177352521,
                -0.2854836465,  0.5043194784
                ])  

                last_omega_0 =  np.array([
                    -0.942753, -1.636921,  0.845378,  0.626647, -2.891469,
                    -1.138307, -2.444745, -2.783882, -3.296146, -0.710319,
                    -2.048366, -1.383319, -1.597016,  0.942246, -3.377102,
                    -1.598180, -2.376562, -0.379636, -0.010647, -2.931106,
                    -1.752307, -1.598185, -0.077504, -2.335897, -1.035773,
                    -2.520229, -0.932417, -3.327573, -3.252453,  1.332111,
                    -1.990557, -2.216876, -1.100975, -0.468819,  1.928117,
                    -2.032934, -2.636557, -1.090020,  2.277180,  8.917205,
                    0.008817, -0.084557, -0.244818, -1.264036,  0.139166,
                    -0.348725, -0.058840,  0.265762,  0.250225, -0.068912,
                    1.704066,  0.107186, -0.322247, -0.009288,  0.012943,
                    -0.013015, -0.171456, -0.019942,  0.001446,  0.008834,
                    0.002107, -0.011222, -0.001785,  0.067749,  0.027791,
                    0.003609,  0.044249, -0.001760, -0.048040,  0.058551,
                    -0.097508,  0.004446,  0.000076, -0.018042, -0.034481,
                    -0.022481,  0.007052, -0.105499, -0.036334, -1.174516,
                    -0.406261, -0.214435, -0.204659,  2.959618,  1.272712,
                    -0.468257,  1.492037
                ])

        elif K==3 and not mle_input_dict["eta_out"]:
            last_omega_0 = np.array( [
            -0.407411493, -2.591162950,  4.291152134,  2.405320336, -4.743821164,
            -0.797324402, -3.402076121, -4.083176409, -5.400690808,  0.047455299,
            -2.595707639, -0.954384693, -0.193192209,  3.662625989, -5.976142528,
            -1.364086915, -3.321703586,  0.508138923,  2.224818533, -4.370287966,
            -1.720920630, -1.934678749,  1.350634383, -2.711605816,  0.462392838,
            -3.448579550,  0.898825736, -5.575813486, -4.778151455,  4.565892764,
            -1.805228149, -3.224209754, -0.789500131,  1.160300083,  5.104865990,
            -1.548216672, -4.135140689, -0.149393265,  5.453167699, 13.003500159,
            -0.086603428,  0.606461033,  0.273212237, -0.813645205,  0.466434224,
            -0.780353119,  0.034402790,  0.532882233,  0.474955632, -0.109399081,
            3.609010894, -0.219184249, -0.812417606, -0.027309281, -0.019704053,
            -0.076299924, -0.242018454, -0.009030303,  0.029345192, -0.003242337,
            -0.004621648, -0.005738294, -0.023424982,  0.064491638,  0.051151748,
            0.014459775,  2.632499789, -0.081991544,  0.004148379,  0.454872575,
            0.028032585,  0.105388184,  0.100008993, -0.063598924, -0.155382435,
            -0.091529274,  0.006808806, -0.523505390, -0.305890002,  0.475472536,
            0.400041707,  0.978100996,  0.613058936,  2.217432011,  0.120037841,
            -0.782908227,  0.269805617
            ])  

            last_omega_0 =  np.array([
            0.380921, -2.630633, 4.784158, 2.519803, -4.166717, -0.435462, -3.504696,
            -4.120323, -5.620902, 0.479236, -2.502938, -0.336962, -0.221202, 4.347957,
            -5.392543, -1.270361, -3.029415, 0.955110, 2.387914, -4.176130, -1.715212,
            -1.898262, 1.479264, -2.206856, 0.998908, -3.120860, 1.009419, -4.927807,
            -5.093832, 5.303804, -1.893854, -2.943071, -0.776752, 1.436988, 5.819950,
            -1.492123, -3.928668, 0.261559, 6.057055, 14.084856, -0.125970, 0.361664,
            0.330792, -0.696620, 0.664917, -0.604086, 0.029160, 0.557371, 0.594538,
            -0.112917, 3.812333, -0.399081, -0.888914, -0.018617, -0.015004, -0.099481,
            -0.245605, -0.038646, 0.014712, -0.009988, -0.011418, -0.015030, -0.025217,
            0.052859, 0.071521, 0.015756, 3.045541, -0.015753, -0.104674, 0.510568,
            0.064976, 0.148560, 0.102977, -0.039297, -0.193619, -0.089968, 0.013582,
            -0.711415, -0.323713, 0.551019, 0.375340, 0.362457, 0.083803,
            3.145111, 1.709262, -1.251001, 3.585348
        ])
        
        elif K==2 and not mle_input_dict["eta_out"]:
              last_omega_0 = np.array( [  -0.400833284, -2.611656703, 4.296780516, 2.440372581, -4.833043355,
                    -0.720298664, -3.504644524, -4.121032452, -5.379819642, 0.037375500,
                    -2.681177873, -0.994054573, -0.224430045, 3.707744809, -6.006841140,
                    -1.355538574, -3.385970799, 0.535331396, 2.229492773, -4.317572287,
                    -1.708824395, -1.939851844, 1.336145363, -2.679377833, 0.483351054,
                    -3.480310412, 0.917710381, -5.589586363, -4.734127133, 4.605721938,
                    -1.835093204, -3.223163923, -0.806974079, 1.184202481, 5.135906134,
                    -1.543211928, -4.144221443, -0.172082654, 5.481942473, 12.985857231,
                    -0.158704949, 0.586385366, 0.287528812, -0.822882192, 0.516235025,
                    -0.834260099, -0.037027887, 0.473625444, 0.498625916, -0.101171828,
                    3.639313188, -0.207594670, -0.810588957, -0.025968450, -0.012797736,
                    -0.071835889, -0.242461593, -0.011846078, 0.031484073, -0.005826691,
                    -0.001523217, -0.007531544, -0.023108395, 0.061276280, 0.053326506,
                    0.014439075, 2.638397676, -0.093120825, 0.060986354, 0.491874658,
                    0.046317294, 0.114883119, 0.083615180, -0.102295099, -0.159868642,
                    -0.092787953, 0.009554377, -0.512577930, -0.293276905, 0.476935534,
                    0.480271214, 0.561717453, 2.339803588, -0.666311760
            ])  
              
              last_omega_0 = np.array([
                  -0.101892, -2.868711, 4.453104, 2.371892, -4.635641, -0.794283, -3.467362,
                    -3.803538, -5.492338, 0.126414, -2.478574, -0.780871, -0.478816, 3.788586,
                    -5.931695, -1.475414, -3.256448, 0.751192, 2.080675, -4.420392, -1.887375,
                    -1.868628, 1.269005, -2.746726, 0.747317, -3.434557, 0.690869, -5.217615,
                    -4.997231, 4.689237, -1.822805, -3.112853, -1.100172, 1.085385, 5.113574,
                    -1.624092, -3.977442, -0.019206, 5.462147, 13.086923, 0.086613, 0.743469,
                    0.164377, -0.856201, 0.563157, -0.695737, 0.069964, 0.598348, 0.517248,
                    -0.118687, 3.588028, -0.275730, -0.813530, -0.011791, -0.019829, -0.081370,
                    -0.246876, -0.024815, 0.025082, 0.002675, -0.010033, -0.011905, -0.019480,
                    0.038415, 0.070000, 0.014694, 2.930697, 0.039094, -0.036920, 0.465981,
                    0.020119, 0.133076, 0.061565, -0.020236, -0.156718, -0.119364, 0.006632,
                    -0.679909, -0.309880, 0.535384, 0.553221, 0.461341, 2.358516, -0.668596
              ])
    
        elif K==2:
            last_omega_0 =  np.array([
                -1.3902372805, -1.8526023226, 0.3009050142, 0.1998554291, -2.9830671611,
                -1.4328039414, -2.5842448974, -2.9175638424, -3.3778959143, -1.0457194632,
                -2.2269434869, -1.6815558951, -1.7896966710, 0.3932910323, -3.4813255064,
                -1.8548114027, -2.5444084357, -0.7304388678, -0.3855495703, -3.0635995101,
                -1.9992969869, -1.8647097967, -0.4580688080, -2.5111326680, -1.3014435782,
                -2.7304332511, -1.2855214293, -3.4199849808, -3.3506015071, 0.7493261822,
                -2.1805915388, -2.3581213562, -1.4150999855, -0.8281222195, 1.2769356337,
                -2.2037145715, -2.7442226824, -1.4026167572, 1.5170770727, 7.6319744157,
                -0.0192539800, -0.0957056696, -0.2241712782, -1.1580334482, 0.1020350237,
                -0.3270057895, -0.0661048634, 0.2483187920, 0.1990288276, -0.0660201116,
                1.5564785101, 0.1608986057, -0.2832288274, -0.0056574025, 0.0084258327,
                -0.0129952967, -0.1695916393, -0.0215096999, 0.0043080746, 0.0087984647,
                0.0012780795, -0.0082392517, -0.0008712602, 0.0689626900, 0.0245605940,
                0.0027709174, -0.2784694117, -0.0034392788, -0.0574381418, 0.0644226815,
                -0.1147370397, -0.0012846052, -0.0014528448, -0.0173120966, -0.0363895077,
                -0.0253499420, 0.0074924140, -0.0892312507, -0.0355494821,
                -1.5609555951, -0.2373604358, -0.1076556597, 2.5353127625, -0.2317634763
            ])
            
            last_omega_0 = np.array([
            -0.334677, -1.144551, 1.457349, 1.268105, -2.404609, -0.590302, -1.981747,
            -2.294407, -2.778581, -0.159959, -1.538968, -0.799822, -1.080445, 1.526922,
            -2.899392, -1.045288, -1.872876, 0.156975, 0.524901, -2.432136, -1.215572,
            -1.079824, 0.469828, -1.818183, -0.520124, -2.023893, -0.318745, -2.799358,
            -2.738927, 1.900119, -1.470902, -1.750559, -0.486922, 0.069170, 2.480572,
            -1.548873, -2.155838, -0.495823, 2.936726, 10.881447, -0.002135, -0.114684,
            -0.183708, -1.371948, 0.145011, -0.285355, -0.111175, 0.217547, 0.215976,
            -0.062986, 1.470508, 0.096824, -0.333670, -0.006092, 0.013552, -0.012230,
            -0.148046, -0.020087, 0.008516, 0.016027, 0.001525, -0.012997, -0.001807,
            0.068795, 0.030276, 0.003801, 0.175129, -0.000638, -0.047147, 0.056441,
            -0.098665, 0.004915, 0.003294, -0.019280, -0.034034, -0.023987, 0.007229,
            -0.103237, -0.034005, -1.032338, -0.459127, -0.429350, 2.576849, -0.228445
        ])
                
        elif K==1 and not mle_input_dict["eta_out"]:
            last_omega_0 = np.array([
                2.366074, 0.059213, 5.148037, 3.386800, -0.547653, 1.727686, 0.100660,
                -0.460113, -0.909345, 2.107812, 0.220875, 2.432544, 1.268625, 4.949938,
                -1.142202, 0.997456, 0.032869, 2.460724, 2.818388, -0.471590, 1.117054,
                1.160102, 2.469760, 0.398559, 1.997084, 0.044519, 2.524530, -0.755548,
                -0.876547, 5.152577, 0.439670, -0.323729, 1.670135, 2.518843, 5.825523,
                0.611124, -0.506067, 2.300430, 5.917927, 14.270349, 0.038592, 0.003434,
                0.462015, -0.057329, 0.460140, -0.193103, -0.191469, 0.147477, 0.251984,
                -0.025049, 1.006159, -0.425352, -0.642877, -0.020939, 0.027263, -0.042517,
                -0.185480, -0.028885, -0.014818, 0.021532, 0.007733, -0.011456, -0.012382,
                0.091015, 0.050773, 0.009112, 3.351653, -0.016316, -0.033439, 0.502317,
                -0.022185, 0.104259, -0.007480, -0.128264, -0.172535, -0.091730, 0.022054,
                -0.840821, -0.299139, 0.623408, -0.045544
            ])
             
            last_omega_0 = np.array([
                -0.068975004, -1.447333899,  2.844808111,  1.072838153, -2.661174841,
                -0.573936846, -2.166927622, -2.519570570, -2.937558959, -0.046437078,
                -1.750056736, -0.352152961, -1.116669271,  2.150368872, -3.160588540,
                -1.033052790, -1.916995346,  0.097804001,  0.724322725, -2.731545093,
                -1.202978966, -0.980772541,  0.472389272, -1.787350185, -0.426306389,
                -2.203462459,  0.237086735, -3.050132430, -3.405066597,  2.656325912,
                -1.760090714, -1.969298001, -0.607400016,  0.155634869,  3.277725602,
                -1.513720498, -2.489077770, -0.055862617,  3.448589425, -7.156742251,
                0.046678249, -0.155187878, -0.028811525, -0.523876161,  0.437513255,
                -0.109379964, -0.118051325,  0.298443345,  0.313427166, -0.049280858,
                1.945040069, -0.288322795, -0.539271921, -0.006088469,  0.010106644,
                -0.049649098, -0.219644721, -0.021307183, -0.049060619,  0.025914977,
                0.006311011, -0.024338082, -0.008002563,  0.116228318,  0.059876532,
                0.005761648,   1.532301006,  0.038587742,
                -0.198292777,  0.287125269, -0.130130893,  0.092255119,  0.037081392,
                0.038740851, -0.109157800, -0.073034338,  0.013687356, -0.332951248,
                -0.182517434, 0.235238912,  0.663420768
            ])
            if not mle_input_dict["in_mag_ind"]:
                last_omega_0 = np.array([
                1.757813, -0.275309,  4.468666,  2.625886, -1.180284,
                1.218876, -0.556468, -1.084455, -1.426018,  1.456760,
                -0.284281,  1.649701,  0.629966,  4.179209, -1.639583,
                0.541759, -0.652399,  1.820846,  2.279104, -1.054750,
                0.474384,  0.586384,  1.886237, -0.260888,  1.414933,
                -0.542286,  1.990714, -1.377658, -1.466098,  4.327446,
                -0.080395, -0.805421,  1.026565,  1.915661,  5.131882,
                -0.064497, -1.051424,  1.684564,  5.202625, 12.822899,
                0.002242, -0.057761,  0.392945, -0.211726,  0.401905,
                -0.212774, -0.160073,  0.203410,  0.229029, -0.033225,
                1.149670, -0.289927, -0.612505, -0.018947,  0.037571,
                -0.030789, -0.179427, -0.023423, -0.019785,  0.013797,
                0.005567, -0.007292, -0.010306,  0.106620,  0.046079,
                0.008238,  2.707738, -0.013649, -0.066773,  0.471953,
                -0.017271,  0.095667, -0.032138, -0.107081, -0.145892,
                -0.088672,  0.021824, -0.678287, -0.246365,  0.497195,
                0.216963
                ])
             

             
        elif K==1:
            last_omega_0 =  np.array([
                -1.8020199362, -2.0417599741, 0.0258693350, -0.5148827827, -3.1424253168,
                -1.8092370473, -2.9033689457, -3.1470215892, -3.4217199845, -1.4120645824,
                -2.6575439949, -2.1101033491, -2.3299919827, -0.2781485669, -3.5319666703,
                -2.2227801671, -2.7683200913, -1.2703922788, -0.9448036174, -3.1823943312,
                -2.3371494490, -2.1412835329, -0.9674365497, -2.7076148061, -1.8909599454,
                -2.9604556396, -1.6744139532, -3.5094114587, -3.5532876655, 0.1279092627,
                -2.5986784037, -2.6618871281, -1.8238577293, -1.4039209859, 0.6763729763,
                -2.5826393929, -2.9413521035, -1.7486182158, 0.8125217763, 6.6648124917,
                -0.0693557927, -0.0740581051, -0.2742839662, -0.9044585752, 0.1286865905,
                -0.2625640391, -0.1817705678, 0.2234651322, 0.2258579749, -0.0559432275,
                1.5543512131, -0.0305002934, -0.2637328973, -0.0096651133, 0.0044494216,
                -0.0124639038, -0.1956389171, -0.0154225095, -0.0060645715, 0.0177175039,
                0.0078790818, -0.0128250042, -0.0007910889, 0.0775072183, 0.0298790156,
                0.0017529399,  -0.4932604898, -0.0102590074,
                -0.0823533822, 0.0661900066, -0.1409868811, 0.0059519887, 0.0026665811,
                -0.0035830008, -0.0341578254, -0.0233338644, 0.0061295755, -0.1026614631,
                -0.0377900841, -1.8680058913, 0.3894567027
            ])

         
            last_omega_0 =  np.array([
            -1.494406, -2.171558, -0.094506, -0.365095, -3.205346, -1.734156, -2.853284,
            -3.154746, -3.571787, -1.361581, -2.611712, -2.018439, -2.272089, -0.080715,
            -3.688528, -2.119185, -2.789832, -1.175442, -0.838742, -3.315191, -2.259108,
            -2.111499, -0.906390, -2.707572, -1.781871, -2.879329, -1.642121, -3.633562,
            -3.560114, 0.275445, -2.550220, -2.738969, -1.781161, -1.277361, 0.928003,
            -2.589565, -3.021496, -1.706511, 1.076527, 6.715349, -0.039140, -0.202234,
            -0.243138, -0.925176, 0.150681, -0.233330, -0.013930, 0.229480, 0.251425,
            -0.055356, 1.597221, -0.046047, -0.290902, -0.009173, 0.012544, -0.013812,
            -0.181331, -0.020399, -0.014092, 0.006131, 0.004890, -0.008858, -0.000859,
            0.080516, 0.022203, 0.002647, -0.349280, -0.005592, -0.047915, 0.065488,
            -0.103035, 0.000825, -0.014969, -0.023705, -0.033795, -0.022239, 0.006991,
            -0.120067, -0.034389, -1.640283, 0.474784
            ])

            if not mle_input_dict["in_mag_ind"]:
                last_omega_0 = np.array([
                0.743856,  0.107301,  2.522119,  1.821326, -0.790873,
                0.835469, -0.372329, -0.684404, -0.769669,  1.112745,
                -0.082728,  0.645163,  0.075888,  2.335441, -0.871242,
                0.406338, -0.259175,  1.103318,  1.442234, -0.473009,
                0.234823,  0.420413,  1.364365, -0.179670,  0.493014,
                -0.417971,  0.820018, -0.679722, -0.822955,  2.623647,
                -0.123736, -0.438890,  0.879138,  1.072664,  3.059364,
                -0.234519, -0.578993,  0.857267,  3.784638, 11.665466,
                0.004399, -0.210100, -0.021522, -1.295681,  0.110563,
                -0.061822, -0.107490,  0.107809,  0.079179, -0.025436,
                0.575374, -0.003575, -0.302410, -0.011931,  0.031625,
                0.004526, -0.046994, -0.016697, -0.006849,  0.018486,
                -0.000157, -0.010424, -0.002045,  0.069383,  0.036743,
                0.002773,  0.326570, -0.003099, -0.033519,  0.056313,
                -0.094234,  0.004122,  0.002555, -0.018920, -0.023083,
                -0.023371,  0.005341, -0.085832, -0.009627, -0.864500,
                -3.840447
            ])


        if len(last_omega_0) == len(omega):
            omega = last_omega_0.copy()
        else:
            raise ValueError(
                f"Starting values have the wrong dimension: "
                f"expected {len(omega)}, got {len(last_omega_0)}"
            )


    if len(omega) != len(omega_names):
        print("Warning: len(omega) != len(omega_names). The log file is likely to display incorrect values of omega")
        print(f"len(omega): {len(omega)}")
        print(f"len(omega_names): {len(omega_names)}")

    return omega, omega_names


# Helper function for parsing MLE params from a numpy array
# Convert the params as a PyTorch tensor
def parse_params(omega, mle_inp):

    Q_X = mle_inp["Q_X"] 
    Q_D = mle_inp["Q_D"]
    Q_W = mle_inp["Q_W"]
    J   = mle_inp["J"]
    K   = mle_inp["K"]
    homog = mle_inp["homog"]
    device = mle_inp["device"]

    # Extract MLE params
    n_delta = 1 if homog else J
    cursor = 0
    delta_j = omega[cursor : cursor + n_delta]
    cursor += n_delta

    beta_x = omega[cursor : cursor + Q_X]
    cursor += Q_X

    psi_d = omega[cursor : cursor + Q_D]
    cursor += Q_D

    mu_x_cost = omega[cursor : cursor + Q_W]
    cursor += Q_W

    sigma_eta = np.exp(omega[cursor])  # scalar
    cursor += 1

    # The rest is theta_parms
    ## std dev
    sigma_theta = np.exp(omega[cursor: cursor + K])
    cursor += K

    ## type share
    if K > 1:
        alpha_k = omega[cursor:(cursor+K-1)]
        alpha_k = np.append(alpha_k, 0) # We normalize alpha of last type to 0 
        cursor += (K - 1)
    else:
        alpha_k = 0

    # Use soft-max to convert alphas to type pr
    exp_alpha_k = np.exp(alpha_k)
    
    if K > 1:
        sum_exp_alpha_k = sum(exp_alpha_k)
        p_theta = exp_alpha_k/sum_exp_alpha_k
    else:
        p_theta = 1

    ## mean utility of each type
    if K > 1:
        mu_theta = omega[cursor:]
        mu_last  = -1*np.dot(p_theta[:(K-1)], mu_theta)/p_theta[K-1] # We normalize mu_K = -(sum^{K-1}_{j=1} p_j mu_j)/p_K
        mu_theta = np.append(mu_theta, mu_last)
    else:
        mu_theta = np.array(0)
    
    trace_coefs = False
    if trace_coefs:
        print("J-specific mean utilities:")
        print(delta_j)

        print("Beta for X:")
        print(beta_x)

        print("Distance costs:")
        print(psi_d)

        print("Individual-level costs:")
        print(mu_x_cost)

        print("sd(eta):")
        print(sigma_eta)

        print("sd(theta):")
        print(sigma_theta)

        print("Type share:")
        print(p_theta)

        print("Mean utility for each type:")
        print(mu_theta)

    return {
        "delta_j": torch.tensor(delta_j, dtype=torch.float32, device=device),
        "beta_x": torch.tensor(beta_x, dtype=torch.float32, device=device),
        "psi_d": torch.tensor(psi_d, dtype=torch.float32, device=device),
        "mu_x_cost": torch.tensor(mu_x_cost, dtype=torch.float32, device=device),
        "sigma_eta": torch.tensor(sigma_eta, dtype=torch.float32, device=device),
        "sigma_theta": torch.tensor(sigma_theta, dtype=torch.float32, device=device),
        "p_theta": torch.tensor(p_theta, dtype=torch.float32, device=device),
        "mu_theta": torch.tensor(mu_theta, dtype=torch.float32, device=device),
        "exp_alpha_k": torch.tensor(exp_alpha_k, dtype=torch.float32, device=device) # Later used for gradient calculations
    }


# Tuple of objfn value and gradients for scipy minimization
def objfn_l_and_grad(omega, omega_names, mle_input_dict, sims_theta, sims_c, n_block = 1, grad=True):
    if n_block == 1:
            wgt_l_i, g_i = ind_l_and_grad(omega, mle_input_dict, sims_theta, sims_c, grad)
    else:
        all_l_i = []
        all_g_i = []
        b_cnt = 1
        for block_dict in mle_input_dict:
            device = block_dict['device']
            
            print(f"Processing block {b_cnt} with {device}...")
            with torch.no_grad():
                l_i_block, g_i_block = ind_l_and_grad(omega, block_dict, block_dict["sims_theta"], block_dict["sims_c"], grad)

            # Save log-likelihoods and sum of gradients per block
            all_l_i.append(l_i_block.detach().cpu())                   # (N_block,)
            if grad:
                all_g_i.append(g_i_block.detach().cpu())  # (N_param,)

            b_cnt +=1
        # Combine
        wgt_l_i = torch.cat(all_l_i, dim=0)                  # shape (N_total,)

        if grad:
            g_i = torch.cat(all_g_i, dim=0)  

    # For tracing problematic students
    num_nan_lik = torch.isnan(wgt_l_i).sum().item()
    num_zero_lik = (wgt_l_i <= 0).sum().item()

    print(f"Number of NaNs in likelihood: {num_nan_lik}")
    print(f"Number of zero values in likelihood: {num_zero_lik}")
    
    if num_nan_lik >0:
        nan_lik_idx = torch.nonzero(torch.isnan(wgt_l_i), as_tuple=False).view(-1)
        print("Indices of NaNs in likelihood:", nan_lik_idx.tolist()[0])
    
    if num_zero_lik > 0:
        zero_lik_idx = torch.nonzero(wgt_l_i <= 0, as_tuple=False).view(-1)
        print("Indices of non-positive in likelihood:", zero_lik_idx.tolist()[0])


    if grad:
        num_nan_grad = torch.isnan(g_i).sum().item()
        print(f"Number of NaNs in gradients: {num_nan_grad}")

        if num_nan_grad > 0:
            nan_grad_idx = torch.nonzero(torch.isnan(g_i), as_tuple=False)
            nan_grad_tuples = [tuple(idx.tolist()) for idx in nan_grad_idx]
            print("Indices of NaNs in g_i (row, col):", nan_grad_tuples[0])
    
    # Replace non-positive likelihood with eps before log transformation
    # Also, remove NaNs/Infs if needed
    eps = 1e-4
    wgt_l_i_clipped = torch.where(wgt_l_i <= 0, torch.tensor(eps, dtype=wgt_l_i.dtype, device=wgt_l_i.device), wgt_l_i)
    sum_neg_ll = -torch.log(torch.nan_to_num(wgt_l_i_clipped, nan=eps, posinf=eps, neginf=eps)).sum().item()
    lik = wgt_l_i.detach().cpu()
    lik = lik[torch.isfinite(lik)]

    print("Individual likelihood stats:")
    print(f"Min     : {lik.min().item():.6e}")
    print(f"Max     : {lik.max().item():.6e}")

    quantiles = torch.quantile(
        lik,
        torch.tensor([0.01, 0.25, 0.5, 0.75, 0.99], dtype=lik.dtype)
    )
    print(f"1% quantile  : {quantiles[0].item():.6e}")
    print(f"25% quantile : {quantiles[1].item():.6e}")
    print(f"50% quantile : {quantiles[2].item():.6e}")
    print(f"75% quantile : {quantiles[3].item():.6e}")
    print(f"99% quantile : {quantiles[4].item():.6e}")


    # gradients
    if grad:
        g_val = -1* torch.nan_to_num(g_i, nan=eps, posinf=1e6, neginf=-1e6).sum(axis=0)
    else:
        g_val = None
        
    print("=============================================")
    print(f"Objfn Value: {sum_neg_ll}")
    
    if grad:
        print(f"Gradients: {g_val}")
        # Norm of the average gradients
        grad_norm = torch.norm((g_val).mean(axis=0))
        print(f"Norm of Average Gradient: {grad_norm}" )

    if grad:
        print("Current coefficients:")  
        for name, val in zip(omega_names, omega):
            print(f"  {name}: {val:.6f}")
        print("=============================================")

    if grad:
        return sum_neg_ll, g_val
    else:
        return sum_neg_ll



# Returns N by J matrix, in which (n,j) is the fixed components of
# utility from enrolling into school j by student n.
# By fixed component, I refer to utility components, excluding theta_ik
def calc_fixed_u(omega_dict, mle_inp):
    
    delta_j    = omega_dict["delta_j"]
    beta_x     = omega_dict["beta_x"]
    psi_d      = omega_dict["psi_d"]
    X          = mle_inp["X"]
    D          = mle_inp["D"]

    # Calculate fixed components of enrollment utilities (N by J)
    # u_bar[n,j] = delta_j[j] + \sum_x^Q_X beta_x*X[n,x] + \sum_k^Q_d psi_d[k]*D[n,k,j]
    
    # Delta j : Matrix of J-specific mean utilities
    delta_term = delta_j.unsqueeze(0) # 1 by J

    # X \beta_x: Decision makers' mean utility
    xb_term = X @ beta_x               
    xb_term = xb_term.unsqueeze(1)  # N by 1
    
    # sum_k psi_d[k] * D[n,k,j] : Distance effects
    d_term = torch.tensordot(D, psi_d, dims=([1], [0]))  # (N, J)

    # Sum of fixed components of utility
    # Note that delta_term is 1 by J and xb_term is N by 1 array.
    # When 1 by J array is added to N by 1 array, R creates N by J array, 
    # in which (n, j) th element = delta_term[1, j] + xb_term[n, 1]
    # This process is called broadcasting
    u_bar = delta_term + xb_term + d_term          

    return u_bar # shape (N, J)


# Calculate fixed components of application costs
# i.e., cost excluding effects from eta_ij
def calc_fixed_cost(omega_dict, mle_inp):
    mu_x_cost = omega_dict["mu_x_cost"]
    W         = mle_inp["W"]

    c_bar = torch.exp(torch.clamp((W @ mu_x_cost), max=600.0))

    return c_bar # Shape (N)

# Helper function for ind_l_and_grad(.)
# Returns weighted average of likelihood by type share stored in omega_dict["p_theta"] 
def l_and_grad_by_type(omega_dict, mle_inp, u_bar, c_bar, sims_theta, sims_c, grad):

    # Extract params -----------------------------------------------------------------------S
    sigma_eta   = omega_dict["sigma_eta"]
    mu_theta    = omega_dict["mu_theta"]
    sigma_theta = omega_dict["sigma_theta"]
    p_theta     = omega_dict["p_theta"]
    N           = mle_inp["N"]
    J           = mle_inp["J"]
    K           = mle_inp["K"]
    R           = sims_theta.shape[1]
    eta_out     = mle_inp["eta_out"]
    lmbda       = mle_inp["lmbda"]

    # Transform simualted draws for theta according to its mu and sigma --------------------
    # theta[n, r, k] = mu_theta[k] + sigma_theta[k]*sims_theta[n, r]
    theta = mu_theta.view(1, 1, -1) + sigma_theta.view(1, 1, -1) * sims_theta.unsqueeze(-1)

    # Calculate n's utility of enrolling into school j when draw = r and type = k
    # (during application stage prior to realization of epsilon)
    # u_apply[n, j, r, k] = u_bar[n, j]  + theta[n, r, k]
    # exp_u_apply = exp_U_1 in the original R script
    # u_apply     = u_tilde in the original R script
    u_apply = u_bar[:, :, None, None] + theta[:, None, :, :]
    u_apply = u_apply.to(torch.float64)  # Convert to float64 before exponentiating
    exp_u_apply = torch.exp(torch.clamp(u_apply, max = 600.0))

    # Free up some memory
    del u_apply

    if mle_inp["device"].type == "cuda":
        torch.cuda.set_device(mle_inp["device"].index)
        with torch.no_grad():
            torch.cuda.empty_cache()

    # Similarly, transform eta
    # eta[n, j, r] = sigma_eta * sims_c[n,j, r]
    # i.e., eta[n,j,r] is equal to student n's eta term for school j for simulation round r

    #print(sims_c[868,:,57])

    eta = sigma_eta * sims_c 
    eta = eta.to(torch.float64)
    # Calculate total costs based on the model specification -------------------------------
    if eta_out:
        # Cost = exp(.) + eta
        c_total = c_bar[:, None, None] + eta 
    else: 
        # Cost = exp( . + eta)
        c_total = c_bar[:, None, None] * torch.exp(torch.clamp(eta, max=600.0))

    # Calculate app-stage expected utility of max(portfolio) -------------------------------
    # e_max is a N by J by r by k matrix, where 
    # e_max[n, j, ,] is student n's expected utility of school j+1, conditional on offer from j+1
    # r is for simulation draw, k is for student type

    # e_max[n,j, r, k] = ln(exp(0)  exp_u_apply[n, j, r, k]) + euler
    euler_cons = 0.5772156649
    e_max = torch.log(1 + exp_u_apply) + euler_cons

    # Now calculate net expected utility of each application -------------------------------
    # e_net_util[n, j, k, r] = (1/lmbda)*pi_i[n, j] * e_max[n,j, k, r] +
    #                          (1/lmbda)*(1-pi_i[n, j]) * euler_cons   -
    #                          c_total[n, j , r]/lmbda
    pi = mle_inp["pi_i"].unsqueeze(2).unsqueeze(3)            

    # Note: intentionally distribute the 1/lmbda first, as it helps numerical precision when e_max is really small
    e_net_util = (1 / lmbda) * pi * e_max + (1/lmbda) * (1 - pi) * euler_cons - c_total.unsqueeze(3)/lmbda
    
    # Create the (N, 1, R, K) tensor filled with euler_cons / lmbda
    # Note that E(u(not apply)) = . 57721
    e_not_apply = torch.full((N, 1, R, K), euler_cons / lmbda, dtype=e_net_util.dtype, device=e_net_util.device)

    # Concatenate it to the left of e_net_util along the school dimension (dim=1)
    e_net_util = e_net_util.to(torch.float64)
    e_net_util = torch.clamp(e_net_util, max=600.0)  # 80 is safe for float64
    e_net_util = torch.cat([e_not_apply, e_net_util], dim=1)
    

    # Free up the memory by deleting e_net_utility which will never be used again
    del e_not_apply

    if mle_inp["device"].type == "cuda":
        with torch.no_grad():
            torch.cuda.empty_cache()

    # Calculate app-stage probabilities ----------------------------------------------------
    e_net_util = torch.exp(e_net_util)

    # p_app is a N by J+1 by r by k tensor
    # p_app[n, j, r, k] =  e_net_util[n, j, r, k] / (sum_{j = 0}^{J+1} e_net_util[n, j, r, k])
    p_app = e_net_util / e_net_util.sum(dim=1, keepdim=True)

    
    # Choose p_a chosen by student n
    A_outside = 1 - mle_inp["A"].sum(dim=1, keepdim=True)              # (N, 1)
    A_full = torch.cat([A_outside, mle_inp["A"]], dim=1)               # (N, J+1)

    # Reshape for broadcasting: (N, J+1, 1, 1)
    A_full = A_full.unsqueeze(2).unsqueeze(3)

    # Select chosen application portfolio probability
    # Elementwise multiplication followed by sum over alternatives
    p_app_i = (p_app * A_full).sum(dim=1)                   # (N, K, R)


    # Calculate enrollment-stage probabilities --------------------------------------------
    # i.e., probability of enrolling into option actually enrolled by students 

    # Only options with offers are considered
    # numer_enroll[n, j, r, k] =  Z[n, j] * exp_u_apply[n, j ,r, k]
    Z = mle_inp["Z"]
    Z_expanded = Z.unsqueeze(2).unsqueeze(3)                          # (N, J, 1, 1)
    numer_enroll = Z_expanded * exp_u_apply                          # (N, J, R, K)

    # denom_enroll[n, r, k] = 1 + sum(num_enroll[n,:,r, k])
    denom_enroll = 1 + numer_enroll.sum(dim=1)                       # (N, R, K)

    # p_enroll[n, j, r, k] = numer_enroll[n, j, r, k]/denom_enroll[n, r, k]
    p_enroll = numer_enroll / denom_enroll.unsqueeze(1)             # (N, J, R, K)

    # concatenate one additional dimension in the 2nd axis dimension (for outside option)
    # p_enroll_out is N by 1 by R by K matrix
    # p_enroll_out[n, 0, r, k] = 1 if rowSums(Z) == 0 # i.e., no offer for n 
    #                          = 1/denom_enroll[n,r,k] else # got offer 
    no_offer = (Z.sum(dim=1) == 0).view(-1, 1, 1) # shape: (N, 1, 1)

    # Create scalar 1 as a float tensor with proper shape
    ones_tensor = torch.tensor(1.0, dtype=denom_enroll.dtype, device=denom_enroll.device)

    # Now safely compute p_enroll_out
    p_enroll_out = torch.where(
        no_offer, 
        ones_tensor,              # scalar
        1.0 / denom_enroll        # shape (N, R, K)
    )
    # Now add the additional application choice axis
    p_enroll_out = p_enroll_out.unsqueeze(1)  # shape: (N, 1, R, K)

    # Concatenate outside option to enrollment probabilities
    p_enroll = torch.cat([p_enroll_out, p_enroll], dim=1)  # shape: (N, J+1, R, K)

    # Multiply by indicator of what student actually did and sum over options
    E_expanded = mle_inp["E"].unsqueeze(2).unsqueeze(3)  # shape: (N, J+1, 1, 1)
    p_enroll_i = (p_enroll * E_expanded).sum(dim=1)  # shape: (N, R, K)

    # liklihood_i[n, r , k] is n's likelihood of applying and enrolling for draw = r and type = k
    likelihood_i =  (p_app_i * p_enroll_i)  # shape: (N, K)

    # Now, calculate the weighted average of simulated likelihoods
    # This is equal to average across r = 1,...,R, then weighted average by type share
    if K > 1:
        wgt_l_i =  likelihood_i.mean(dim=1)@p_theta.to(dtype=torch.float64)
    else:
        wgt_l_i = likelihood_i.mean(dim=1)

    if not grad:
        if K==1:
            wgt_l_i = wgt_l_i.squeeze(-1)
        return wgt_l_i, None

    # -------------------------------------------------------------------------------
    # Calculate gradients -----------------------------------------------------------
    #--------------------------------------------------------------------------------

    # Initial preparation (Create objects freqeuntly used across omega) -------------
    # p_e_z[n, j-1, r, k] = Pr(Choose j if face offer from j)
    #                        =  exp_u_apply[n, j, r, k]/(exp_u_apply[n, j, r, k] + 1)
    
    # I don't create Pr(choose outside if face offer from j) as it's just 1 - pr(choose j | z_ij = 1)
    denom = exp_u_apply + 1
    p_e_z = exp_u_apply / denom
    
    # Multiply the pr() by admission probability
    # The product of those two are part of application-stage gradients for most of the omega
    # grad_emax[n, j, r, k] = p_e_z[n, j, r, k]*pi[n, j]
    # Note that this conincides with gradients of log likelihood w.r.t., an intercept term
    grad_emax = p_e_z*pi

    # Obtain the grad_emax corresponding to the school student applied to.
    # grad_emax_i[n, r, k] =  sum(grad_emax[n, j, r, k] * A[n, j])
    A_expanded = mle_inp["A"].unsqueeze(2).unsqueeze(3)
    grad_emax_i = (grad_emax * A_expanded).sum(dim=1)  # shape: (N, R, K)

    # Second, work on the component of enrollment stage gradients
    # i.e., derivative of pr(choosing applied school | offer) w.r.t, intercept term
    # grad_enroll[n,r, k] = p_e_z[n,j,r,k] *A[n, j]
    grad_enroll =  (p_e_z*A_expanded).sum(dim = 1)



    ## Calculate gradients for each parameter in omega (except alpha k) -------------

    # School mean utility (delta) ----------------------------------------------------------
    has_offer = mle_inp["has_offer"]
    if mle_inp["homog"]:

        # N by J by R by K
        print("Ryan:homog specification code is likely to crash due to dimension mismatch")
        print("But, I didn't fix it because it's unlikely that we set homog==True")
        dlog_p_a_delta =  (grad_emax_i/lmbda) - ((grad_emax/lmbda*p_app[:,1:,:]).unsqueeze(2)*x.unsqueeze(1)).sum(axis=1)
        chose_app_sch = mle_inp["chose_app_sch"].view(-1, 1, 1)
        dlog_p_e_delta =  (chose_app_sch - grad_enroll.unsqueeze(1))*x*has_offer.view(-1, 1, 1, 1)

    else:
        # N by J by R by K
        dlog_p_a_delta = (grad_emax/lmbda * A_expanded) - grad_emax/lmbda*p_app[:,1:,:,:]
        chose_app_sch = mle_inp["chose_app_sch"].view(-1, 1, 1, 1)
        dlog_p_e_delta = (A_expanded*(chose_app_sch - p_e_z)*has_offer.view(-1, 1, 1, 1)) 

        # return grad_emax, A_expanded, p_app,  chose_app_sch

    # Take the average across simulation draws. 
    # Doing this ASAP helps memory usage
    dlog_p_type = ((dlog_p_a_delta + dlog_p_e_delta)*likelihood_i.unsqueeze(1)).mean(axis=2) # N by J by K
    

    # free up memory by removing duplicated objects
    del dlog_p_a_delta, dlog_p_e_delta

    if mle_inp["device"].type == "cuda":
        with torch.no_grad():
            torch.cuda.empty_cache()
    gc.collect()

    # Covariate coefficients (beta) -------------------------------------------------------
    x = mle_inp["X"].unsqueeze(2).unsqueeze(3) # N by Q_X by 1 by 1

    # pi_p_a_i has a shape of N by 1 by R by K
    pi_p_a_i = (grad_emax_i/lmbda).unsqueeze(1)

    # pi_p_ae_j is N by J by 1 by R by K
    pi_p_ae_j = (grad_emax/lmbda*p_app[:,1:,:,:]).unsqueeze(2)

    # app stage and enrollment stage gradients
    # Both objects are N by Q_X by R by K
    dlog_p_a_x =  pi_p_a_i*x - (pi_p_ae_j*x.unsqueeze(1)).sum(axis=1) # sum over j > 0
    dlog_p_e_x =  (chose_app_sch - grad_enroll.unsqueeze(1))*x*has_offer.view(-1, 1, 1, 1)

    # Append gradients and free up memory by removing duplicated objects``

    dlog_p_type = torch.cat([dlog_p_type, ((dlog_p_a_x + dlog_p_e_x)*likelihood_i.unsqueeze(1)).mean(axis=2)], dim =1) 

    del x, dlog_p_a_x, dlog_p_e_x

    if mle_inp["device"].type == "cuda":
        with torch.no_grad():
            torch.cuda.empty_cache()
    gc.collect()


    # Distance coefficients (psi) -------------------------------------------------------
    # D = N by Q_D by J
    # rel dist is N x J x Q_D, 1, 1
    rel_dist = mle_inp["D"].permute(0, 2, 1).unsqueeze(-1).unsqueeze(-1)
    
    # rel_dist_chose[n,q]: distance X characteristics term for the school student n applied to
    rel_dist_chose =(mle_inp["D"]*mle_inp["A"].unsqueeze(1)).sum(axis=-1).unsqueeze(-1).unsqueeze(-1)


    # N by Q_D by R by K
    dlog_p_a_d =  pi_p_a_i*rel_dist_chose - (pi_p_ae_j*rel_dist).sum(axis=1)
    dlog_p_e_d =  (chose_app_sch - grad_enroll.unsqueeze(1))*rel_dist_chose*has_offer.view(-1, 1, 1, 1)

    # Append gradients and free up memory by removing duplicated objects
    dlog_p_type = torch.cat([dlog_p_type, ((dlog_p_a_d + dlog_p_e_d)*likelihood_i.unsqueeze(1)).mean(axis=2)], dim =1) 

    del rel_dist, rel_dist_chose, dlog_p_a_d, dlog_p_e_d

    if mle_inp["device"].type == "cuda":
        torch.cuda.empty_cache()
    gc.collect()

    # Fixed cost coefficients (mu_c_x_cost)----------------------------------------------------
    if eta_out:
        # dc[n, q] = c_bar[n] * w[n, q]
        dc = c_bar[:, None]/lmbda * mle_inp["W"]
        
        # If eta is inside the exp(), the dc doesn't vary across j
        # so chosen_dc is just dc
        # Still, we need to multiply rowSums(A), as students who didn't apply incur no cost
        applied = mle_inp["A"].sum(axis=1)
        chosen_dc = dc* applied[:,None]

        # p_app.shape: (N, J+1, R, K)
        # dlog_p_a_w[n, q_w, r, k] = -(chosen_dc[n, q_w, None, None] - sum_j(dc[n, None, qw, None, None] * p_app[n, j+1, None, r, k]))
        dlog_p_a_w = -(chosen_dc[:, :, None, None] - (dc[:, None, :, None, None]*p_app[:,1:,None, :,:]).sum(axis=1) )

    else:
        # If eta is inside the exp(), the dc varies at the j-level
        # dc[n, j, q, r] = c_total[n, j, None, r] * w[n,None, q_w. None]

        dc = c_total[:,:,None,:]/lmbda * mle_inp["W"][:, None,:,None]

        # chosen_dc.shape = (N, Q_W, R)
        # dc of chosen school
        chosen_dc = (dc * mle_inp["A"][:,:,None,None]).sum(axis=1)
    
        # For chosen dc, add one additional axis for types
        # dlog_p_a_w[n, q_w, r, k] = chosen_dc[n, q_w, r]  - (dc[n, j, q_w, r, None]*p_app[n, j+1, None r, k])
        dlog_p_a_w = -(chosen_dc.unsqueeze(-1) - (dc.unsqueeze(-1)*p_app[:,1:,None,:,:]).sum(axis=1) )
    
    dlog_p_e_w =torch.zeros(N, mle_inp["Q_W"], R, K, dtype=torch.float32, device=mle_inp["device"]) # w doesn't affect enrollments


    dlog_p_type = torch.cat([dlog_p_type, ((dlog_p_a_w + dlog_p_e_w)*likelihood_i.unsqueeze(1)).mean(axis=2)], dim = 1) 

    

    del dc, chosen_dc, dlog_p_a_w, dlog_p_e_w
 

    if mle_inp["device"].type == "cuda":
        torch.cuda.empty_cache()
    gc.collect()

    # log(sd(eta_ij)) ------------------------------------------------------------------------
    #print(eta[868,:,57])
    if (eta_out):
        dc = eta/lmbda
    else:
        # eta[n,j,r] is equal to student n's eta term for school j for simulation round r
        dc = c_total*eta/lmbda

    del eta 
 
    if mle_inp["device"].type == "cuda":
        torch.cuda.empty_cache()
    gc.collect()

    # chosen_dc.shape = (N, R)
    chosen_dc = (dc * mle_inp["A"].unsqueeze(-1)).sum(axis=1)

    #print(chosen_dc[868,:])


    # dc.unsqueeze(-1).shape (N, J, R, 1)
    # p_app.shape: (N, J+1, R, K)
    #print(dc[868,:,:])


    dlog_p_a_sd_eta = -(chosen_dc.unsqueeze(-1) - (dc.unsqueeze(-1)*p_app[:,1:,:,:]).sum(axis=1) )
    dlog_p_e_sd_eta = torch.zeros(N, 1, R, K, dtype=torch.float32, device=mle_inp["device"]) # eta doesn't affect enrollments

    dlog_p_type = torch.cat([dlog_p_type, ((dlog_p_a_sd_eta.unsqueeze(1) + dlog_p_e_sd_eta)*likelihood_i.unsqueeze(1)).mean(2)], dim=1) 


    del dlog_p_a_sd_eta, dlog_p_e_sd_eta, chosen_dc, dc 

    if mle_inp["device"].type == "cuda":
        torch.cuda.empty_cache()
    gc.collect()

    # log(sd(theta_iq)) --------------------------------------------------------------
    # IMPORTANT: This part gets confusing as we have two indicies we need to keep track for types.
    #            I use K for types we have to take weighted average for
    #            I use q for types we have to take derivative w.r.t, the std deviation
    #            du[n, q, r] = sigma_theta[q]*sims_theta[n, r]
    du = sigma_theta[None,:,None] * sims_theta.unsqueeze(1)

    # pi_p_a_i.shape = (N, R, K)
    # pi_p_a_i = (grad_emax_i/lmbda).unsqueeze(1)
    # pi_p_ae_j is N by J by 1 by R by K

    # dlog_p_a_x[n, q, r, k] = pi_p_a_i[n, none, r, k]*du[n, q,r, None] - (pi_p_ae_j[n, j, none ,r ,k]*du[n, None, q, r, None])
    dlog_p_a_sd_theta =  pi_p_a_i*du.unsqueeze(-1) - (pi_p_ae_j*du[:, None, :, :, None]).sum(axis=1) # sum over j > 0

    # dlog_p_e_sd_theta[n, q, r, k] = chose_app_sch[n, 0, 0, 0] - grad_enr[n, , r, k] * du[n, q, r, None]
    dlog_p_e_sd_theta =  (chose_app_sch - grad_enroll.unsqueeze(1))*du.unsqueeze(-1)*has_offer.view(-1, 1, 1, 1)

    # Create a [Q, K] mask where only diagonal elements are 1 (q == k)
    mask_q_eq_k = torch.eye(K, K, device=dlog_p_a_sd_theta.device)

    # Expand to match shape: [1, Q, 1, K] so it broadcasts over N and R
    mask_q_eq_k = mask_q_eq_k.view(1, K, 1, K)

    # Zero out entries where q != k
    dlog_p_a_sd_theta = dlog_p_a_sd_theta * mask_q_eq_k
    dlog_p_e_sd_theta = dlog_p_e_sd_theta * mask_q_eq_k

    
    dlog_p_type = torch.cat([dlog_p_type, ((dlog_p_a_sd_theta + dlog_p_e_sd_theta)*likelihood_i.unsqueeze(1)).mean(axis=2)], dim=1) 



    del dlog_p_a_sd_theta, dlog_p_e_sd_theta, du 

    if mle_inp["device"].type == "cuda":
        torch.cuda.empty_cache()
    gc.collect()

    # alpha_q (Type share) ----------------------------------------------------------
    # This part calculates the term that'll be later used for the product rule
    # Note that all gradients are 0, except for the case when type = K 
    if K>1:
        du = -1*omega_dict["exp_alpha_k"][:(K-1)]*mu_theta[:(K-1)]
        

        # dlog_p_a_alpha_K[n, q, r, k] = pi_p_a_i[n, None, r, k]*du[None, q,  None, None] 
        #                                - sum_j(pi_p_ae_j[n, j, None, r, k)] * du[:,:,q,:,None])

        dlog_p_a_alpha_K = pi_p_a_i*du[None, :, None, None] - (pi_p_ae_j*du[None, None, :, None, None]).sum(axis=1) # sum over j > 0
        dlog_p_a_alpha_K[..., :-1] = 0

        # dlog_p_e_alpha_K[n, q, r, k] = (chose_app_sch -  grad_enroll.unsqueeze(1))*du.unsqueeze(-1)*has_offer.view(-1, 1, 1, 1)
        dlog_p_e_alpha_K = (chose_app_sch - grad_enroll.unsqueeze(1))*du[None, : , None, None]*has_offer.view(-1, 1, 1, 1)
        dlog_p_e_alpha_K[..., :-1] = 0

        dlog_p_type = torch.cat([dlog_p_type, ((dlog_p_a_alpha_K + dlog_p_e_alpha_K)*likelihood_i.unsqueeze(1)).mean(axis=2)], dim=1) 



        del dlog_p_a_alpha_K, dlog_p_e_alpha_K, du 
    

        if mle_inp["device"].type == "cuda":
            torch.cuda.empty_cache()
        gc.collect()

    # mu_theta_k (Mean utility of theat_k distribution) -----------------------------
    if K>1:
        de= torch.ones((K - 1, K))
        de[:, -1] =  -p_theta[:K-1]/p_theta[K-1]
        row_idx = torch.arange(K - 1).unsqueeze(1)  # [K-1, 1]
        col_idx = torch.arange(K - 1).unsqueeze(0)  # [1, K-1]
        mask = (row_idx == col_idx)  # [K-1, K-1]

        # Apply mask to left block (columns 0 to K-2)
        # This mask imposes condition if j == k in the original code
        de[:, :K-1] *= mask

        # return de, pi_p_a_i,chose_app_sch, grad_enroll, has_offer , pi_p_ae_j
        de = de.to(pi_p_a_i.device)

        dlog_p_a_mu_theta =  pi_p_a_i*de[None, :, None, :]  - (pi_p_ae_j*de[None, None, :, None, : ]).sum(axis=1)
        dlog_p_e_mu_theta =  (chose_app_sch.view(-1, 1, 1, 1) - grad_enroll.unsqueeze(1))*de[None, :, None, :] * has_offer.view(-1, 1, 1, 1)

        # return dlog_p_a_mu_theta, dlog_p_e_mu_theta
        dlog_p_type = torch.cat([dlog_p_type, ((dlog_p_a_mu_theta + dlog_p_e_mu_theta)*likelihood_i.unsqueeze(1)).mean(axis=2)], dim=1) 


        del dlog_p_a_mu_theta, dlog_p_e_mu_theta, de
  

        if mle_inp["device"].type == "cuda":
            with torch.no_grad():
                torch.cuda.empty_cache()
        gc.collect()


    # Calculate the final gradients -------------------------------------------------


    if K > 1:
        
         # Do product rule for alpha_k ---------------------------------------------------
        dlog_p_wgt = (dlog_p_type@p_theta.to(dtype=torch.float64))
        cursor = J+mle_inp["Q_X"]+mle_inp["Q_W"]+mle_inp["Q_D"] + K + 1

        # gp[n, q] =  p_theta[q]*(1-p_theta[q])*likelihood_i_bar[, q] for q < K
        likelihood_i_bar = likelihood_i.mean(dim=1)
        gp = ((p_theta*(1-p_theta))[None, :]*likelihood_i_bar)[:,:(K-1)]

        p_l_except_q = (likelihood_i_bar@p_theta.unsqueeze(-1).double()- likelihood_i_bar * p_theta)[:,:(K-1)]

        
        # sum_gp_not_q = sum_except_q[n, q]*p_theta[:q] %*% p_theta[-q]
        sum_gp_not_q = (p_theta[:(K-1)].unsqueeze(0)*p_l_except_q)
        
        dlog_p_wgt[:,cursor:(cursor+K-1)] = gp - sum_gp_not_q +  (dlog_p_type[:,cursor:(cursor+K-1),(K-1)]*p_theta[-1].unsqueeze(-1))

        dlog_p_wgt = dlog_p_wgt/wgt_l_i.unsqueeze(1)
    else:
        dlog_p_wgt = dlog_p_type.squeeze(-1)/wgt_l_i
        wgt_l_i = wgt_l_i.squeeze(-1)



    return  wgt_l_i,  dlog_p_wgt


# Matrix of individual likelihood contributions and its gradients w.r.t. params
def ind_l_and_grad(omega, mle_inp, sims_theta, sims_c, grad=True):
    """_summary_

    Args:
        omega (np.array): An array containing coefficients for the model
        mle_inp (dictionary): Outputs from convert_to_tensor(process_dta()). 
                              Contains objects required for the likelihood and gradients calculations
        sims_theta (tensor.float32): Draws from standard normal. Its shape is N by R
        sims_c (tensor.float32): Draws from standard normal. Its shape is N by J by R

    Returns:
        individual level likelihood (Length N tensor)
        Individual level gradient   (N by len(omega) tensor)
    """
    
    # Parse MLE parameters into tensors
    omega_dict = parse_params(omega, mle_inp)

    # Calculate utilities of enrollments excluding random coefficients
    u_bar = calc_fixed_u(omega_dict, mle_inp)

    # Fixed costs
    c_bar = calc_fixed_cost(omega_dict, mle_inp)
    #print(c_bar[868,])
    
    # Calculate weighted avg of likelihood and gradients across type k <= K
    return l_and_grad_by_type(omega_dict, mle_inp, u_bar, c_bar, sims_theta, sims_c, grad)


# Output gradients and final coefficients as a NPZ file for R scripts
def write_npz(out_dir, mle_fit, omega_names, mle_inp,sims_c, sims_theta, n_blocks, last_yr):

    if n_blocks == 1:
        # obtain individual level gradient tensor at final values of omega
        _, g_i = ind_l_and_grad(omega=mle_fit["x"],
                            mle_inp=mle_inp,
                            sims_theta=sims_theta,
                            sims_c=sims_c,
                            grad=True)
    else:
        all_l_i = []
        all_g_i = []
        b_cnt = 1
        for block_dict in mle_inp:
            print(f"Processing block {b_cnt}...")
            l_i_block, g_i_block = ind_l_and_grad(mle_fit["x"], block_dict, block_dict["sims_theta"], block_dict["sims_c"],grad=True)

            # Save log-likelihoods and sum of gradients per block
            all_l_i.append(l_i_block.detach().cpu())                   # (N_block,)
            all_g_i.append(g_i_block.detach().cpu())  # (N_param,)

            b_cnt +=1
        # Combine
        wgt_l_i = torch.cat(all_l_i, dim=0)                  
        g_i = torch.cat(all_g_i, dim=0)  


    # gradients
    g_val = -1* torch.nan_to_num(g_i, nan=1e-4, posinf=1e6, neginf=-1e6).sum(axis=0)
    
    # Output the datasets as NPZ file
    if n_blocks > 1:
        mle_inp = mle_inp[0]

    K = f"K{mle_inp['K']}" 
    homog = "homog" if mle_inp["homog"] else "hetero"
    eta_out = "eta_out" if mle_inp["eta_out"] else "eta_in"
    in_mag_ind  = "in_mag_ind"  if mle_inp["in_mag_ind"] else "no_in_mag_ind"
    npz_fname = f"{in_mag_ind}_{eta_out}_results_choice_{homog}_{K}_{last_yr}_bfgs.npz"

    np.savez(os.path.join(out_dir, npz_fname), 
             omega_choice_hat = mle_fit["x"],
             omega_names      = omega_names,
             Q_hat = mle_fit["fun"],
             G     = g_val.cpu().numpy(),
             G_i   = g_i.cpu().numpy())



def calc_mle_std_err(mle_fit, mle_inp_dict, sims_theta, sims_c):

    # obtain individual level gradient tensor at final values of omega
    _, g_i = ind_l_and_grad(omega=mle_fit["x"],
                         mle_inp=mle_inp_dict,
                         sims_theta=sims_theta,
                         sims_c=sims_c)
    
    # Calculate vcov matrix of the transformed estimates using delta method
    # Since sd(theta) or sd(eta) can't be negative, we estimated ln(sd(eta))
    # We need to transform the estimated ln(sd) back to sd
    inv_fisher = torch.linalg.inv(g_i.T @ g_i / g_i.shape[0])
    
    # Initialize jacobian matrix for delta method
    # Note: +2 for last type's type share and mu
    n_eye = g_i.shape[1] if mle_inp_dict["K"] < 2 else g_i.shape[1] + 2
    j_mat  = torch.eye(g_i.shape[1], dtype=torch.float64, device=mle_inp_dict["device"])

    # Create a dictionary containing final parameters
    mle_dict = parse_params(mle_fit["x"], mle_inp_dict)

    # Create a tensor of final parameters
    omega = torch.cat([
        mle_dict["delta_j"].flatten(),
        mle_dict["beta_x"].flatten(),
        mle_dict["psi_d"].flatten(),
        mle_dict["mu_x_cost"].flatten(),
    ])

    
    # Work on ln(sd())
    cursor = omega.shape[0]
    mle_dict["sigma_eta"]   = torch.exp(mle_dict["sigma_eta"])
    mle_dict["sigma_theta"] = torch.exp(mle_dict["sigma_theta"])

    # Append transformed sigma values to omega
    omega = torch.cat([
        omega,
        mle_dict["sigma_eta"].flatten(),
        mle_dict["sigma_theta"].flatten()
    ])

    # Update j_mat (Jacobian matrix) for the transformed log(sd) parameters
    num_sigma = mle_dict["sigma_eta"].numel() + mle_dict["sigma_theta"].numel()
    for i in range(cursor, cursor + num_sigma):
        j_mat[i, i] = omega[i]  # Because d(exp(x))/dx = exp(x)
    cursor += num_sigma
    
    # If K>1, need to add (p_1,..., p_K), and (mu_1,..,mu_K)
    # Make similar changes for p and mu
    if mle_inp_dict["K"]>1:
        p = mle_dict["p_theta"]
        mu = mle_dict["mu_theta"]
        omega = torch.cat([omega, 
                           p,
                           mu])
        
        # jacobian of p_k, w.r.t, alpha_k
        # dp_k/d_alpha_j = p_k(1-p_k) if j == k
        #                = -p_k*p_j   if j != k
        #                =  exp(alpha_j) /(1+sum_n exp(alpha_n))^2  if k == K

        # Create outer product
       # outer = p * p

        # Fill in the matrix
        #H = -outer  # Off-diagonal terms: -p_k * p_j
        # H.fill_diagonal_(p_sub * (1 - p_sub))  # Diagonal terms: p_k * (1 - p_k)
                


        

        

        


        

    

    