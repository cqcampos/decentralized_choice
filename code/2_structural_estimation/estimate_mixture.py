from util import *

# Basic setup
SERVER  = 1       # Running code from the server or from Github Repo?
R       = 300     # Number of simulation draws for estimation    
R_post  = 1000    # Number of draws for calculating posterior
lmbda   = 0.05    # Smoothing parameter
homog   = 0       # Homogeneous mean utility (0) or j-specfic mean utility (1)?
eta_out = False   # Is cost = exp() + eta or  exp(eta)?
in_mag_ind =True # Do we include currently in-magnet indicator into X?
device =  torch.device("cuda" if torch.cuda.is_available() else "cpu") ## specify the GPU id's, GPU id's start from 0    # 'cpu' or 'cuda'(for GPU usage)
draw_norm_method = "sobol" # Method for drawing normal
seed_normal = 12345 # Seed for the normal draws
ftol        = 5e-14  # Conclude convergence when |(Objfn(t)- Objfn(t+1))/max(Objfn(t), Objfn(t+1)))| is smaller than this
gtol        = 1    # Stop when all elements of gradient is smaller than this value
hot_start   = False  # Set the starting values of MLE as the pre-chosen values based on the past optimizations
K           = 2 # Number of types
last_yr   = 2013 # Last year for the application
# At most, 3 GPUs can be used in Mercury
# When using CUDA, split the data so that each GPU gets an equal share,
# and the remaining portion is processed on the CPU
# But use 2 GPUs because you can run 2 jobs using 2GPUs, but 1 job using 3 GPUs
n_blocks = (torch.cuda.device_count() + K)  if device == torch.device("cuda") else 1  

# Set paths
if SERVER:
    ROOT = "/project/lausd/decentralized_choice"
else:
    ROOT = "Z:/decentralized_choice"

OUT = os.path.join(ROOT, "estimates")


if __name__ == "__main__":
   for eta_out in [True]:
        # Show number of threads used
        print(f"Optimization specifications =========================================")
        print(f"Number of threads use for torch operation {torch.get_num_threads()}")
        print(f"Number of simulation draws: {R}")
        print(f"Number of types assumed: {K}")
        print(f"Eta outside of exp(): {eta_out}")
        print(f"X includes currently in magnet indicator: {in_mag_ind}")
        print(f"Method for drawing normal: {draw_norm_method}")
        print(f"Seed value for normal draws: {seed_normal}")
        print(f"Number of blocks for objfn and grad: {n_blocks}")
        print(f"Last year of applications: {last_yr}")
        print("======================================================================")
        
        # Read and process raw data for optimization
        print("Loading raw data...")
        mle_input_dict = process_dta(ROOT, K, R, lmbda, homog, eta_out, in_mag_ind, device,  full_samp=False, last_yr=last_yr)
        print(mle_input_dict["N"])
        print(mle_input_dict["J"])
        print(mle_input_dict["A"].shape)
        # Simulate draws for theta and eta
        print("Drawing from normal...")
        sims_theta = draw_normal(N=mle_input_dict['N'], R=mle_input_dict['R'], J=1, method=draw_norm_method, device=device, seed = seed_normal)
        sims_c = draw_normal(N=mle_input_dict['N'], R=mle_input_dict['R'], J=mle_input_dict['J'], method=draw_norm_method, device=device, seed=seed_normal)

        print("Printing the theta draws for round 1-5 for student 1-5")
        print(sims_theta[0:5, 0:5])

        print("Printing some of the sims_c")
        print(sims_c[868,:,min(57, R)])


        # Set starting values for MLE
        # Also returns variable labels for omega_0
        print("Picking starting values for optimization...")
        omega_0, omega_names = set_start_params(mle_input_dict, hot_start=hot_start, last_yr=last_yr)


        # Convert the matrices into Pytorch Tensors
        print("Converting numpy arrays to tensors...")
        if n_blocks == 1:
            mle_input_dict = convert_to_tensor(mle_input_dict, device=device)
        else:
            mle_input_dict = convert_to_tensor_blocks(mle_input_dict, n_blocks, sims_theta, sims_c)
            sims_c = None
            sims_theta = None
            gc.collect()
                    
        # For developing code for objective function values and gradients
        #sum_neg_ll, g_i = objfn_l_and_grad(omega_0, omega_names, mle_input_dict, sims_theta, sims_c)
        #print(f"Objfn Value: {sum_neg_ll}")
        #print(f"Gradients: {g_i}")
   
        # Run nelder-mead to refine starting values of omega
        # Use subset of simulation draws for refining the start points
        conv_ind = 0
        final_objfn = 999999
        np.random.seed(123456)

        #if K == 3:
        #    eta_str = "eta_out" if eta_out else "eta_in"
        #    prev_optim_path = os.path.join(OUT, f"in_mag_ind_{eta_str}_results_choice_hetero_K2_bfgs.npz")
        #    if os.path.exists(prev_optim_path):
        #        prev_optim = np.load(prev_optim_path)
        #         omega_0 = prev_optim["omega_choice_hat"]

        #        # Insert from the back so earlier insertions don't shift later positions
        #        omega_0 = np.insert(omega_0, 84, -0.1)  # mu theta 2
        #        omega_0 = np.insert(omega_0, 83, -0.1)  # alpha theta 2
        #        omega_0 = np.insert(omega_0, 82, -0.2)  # log(sd(theta 3))

        print("starting values:")
        print(omega_0)


        n_attempt = 0
        # Don't terminate unless convergence has been achieved
        while not conv_ind:

            # When the L-BFGS-B fails, use starting values refined by Nelder-Meads
            if n_attempt > 0:
                print("Starting nelder-mead to refine starting values of omega...")
                max_iter = 1500 if K <=2 else 100
                mle_fit_nm = minimize(
                    fun=objfn_l_and_grad,
                    x0=omega_0,                # initial guess (numpy array)
                    method='Nelder-Mead',
                    args=(omega_names, mle_input_dict,  sims_theta[:,:50], sims_c[:,:, :50], n_blocks, False),
                    options={
                    'disp': True,
                    'maxiter':max_iter, 
                    }
                )

                omega_0 = mle_fit_nm["x"]
            

            # Run L-BFGS-B using scipy.minimize()
            print("Starting L-BFGS-B optimization...")
            mle_fit = minimize(
                fun=objfn_l_and_grad,
                x0=omega_0,                # initial guess (numpy array)
                method='L-BFGS-B',
                args=(omega_names, mle_input_dict,  sims_theta, sims_c, n_blocks),
                jac=True,                  # our wrapper returns both objval and gradient
                options={
                'disp': True,
                'maxiter':10000, # 'maxiter': 10000,
                'gtol': gtol,     # Stop when all elements of gradient is smaller than this value
                'ftol': ftol    # Conclude convergence when marginal improvement/max(Objfn(t), Objfn(t+1)) is smaller than this
                }
            )
            
            print(mle_fit)
            print(mle_fit["x"])
            print(mle_fit["hess_inv"])
            
            conv_ind = 1 if mle_fit["status"] == 0 else 0
            final_objfn = mle_fit["fun"]

            # If the optimization failed, add small noise to omega_0 and refine omega_0 with Nelder-Mead
            epsilon = 0.5
            omega_0 = omega_0 + np.random.uniform(-epsilon, epsilon, size=omega_0.shape)
            n_attempt += 1
  
        # Save the output as an npz file
        # Edit this function, so it works for even blocking set up
        write_npz(OUT, mle_fit, omega_names, mle_input_dict, sims_c, sims_theta, n_blocks, last_yr)

        # Clean up tensors and memory
        del mle_input_dict, sims_c, sims_theta, omega_0, omega_names, mle_fit
        gc.collect()  # Collect garbage from CPU memory

        if device.type == "cuda":
            torch.cuda.empty_cache()       # Release cached memory (helps prevent fragmentation)
            torch.cuda.ipc_collect()       # Force inter-process communication cleanup (optional but useful in loops)
