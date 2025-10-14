from util import *

# Basic setup
SERVER  = 0       # Running code from the server or from Github Repo?
R       = 300     # Number of simulation draws for estimation    
R_post  = 1000    # Number of draws for calculating posterior
K       = 2       # Number of types assumed
lmbda   = 0.05    # Smoothing parameter
homog   = 0       # Homogeneous mean utility (0) or j-specfic mean utility (1)?
eta_out = False   # Is cost = exp() + eta or  exp(eta)?
device = "cpu"    # CPU or CUDA (for GPU usage)

# Set input paths
if SERVER:
    ROOT = "Z:/decentralized_choice"
else:
    ROOT = "/project/lausd/decentralized_choice"

if __name__ == "__main__":

    # Read and process raw data for optimization
    mle_input_dict = process_dta(ROOT)

    # Simulate draws for theta and eta
    sims_theta = draw_normal(N=mle_input_dict['N'], R=mle_input_dict['R'], J=1, method='sobol', device=device)
    sims_c = draw_normal(N=mle_input_dict['N'], R=mle_input_dict['R'], J=mle_input_dict['J'], method='Lattice', device=device)

    # Set starting values for MLE
    # Also returns variable labels for omega_0
    omega_0, omega_names = set_start_params(mle_input_dict)

    # Convert the matrices into Pytorch Tensors
    mle_input_dict = convert_to_tensor(mle_input_dict, device=device)

    # TODO:ONCE THE CONVERSION IS DONE, REMOVE THE MATRIX OBJECTS WHICH ARE NO LONGER NEEDED!
    # WE"LL NEED TO SAVE SOME MEMORY

    # Run L-BFGS-B using scipy.minimize()
