import dadi
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime
import sys

# -------------------------------
# Step 1: Load the SFS (updated path)
sfs_file_path = "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/data/chamaea/sfs/dadi/North-South.sfs"
data = dadi.Spectrum.from_file(sfs_file_path)  # Ensure this file exists!
# data = data.fold()  # Uncomment if you need to fold the SFS

ns = data.sample_sizes

# Print number of SNPs in the dataset
num_snps = int(data.S())
print(f"Number of SNPs in projected dataset: {num_snps}")

# -------------------------------
# Step 2: Define the models



#######################################################
#vicariance scenarios
###################################################
#Model 1: vicariance no migration
def vic_no_mig(params, ns, pts):
    """
    Split into two populations, no migration. Populations are fractions of the reference 
    population (Nref = 1). Population 1 (North) is the island/founder population with size s * Nref.
    Population 2 (South) is the mainland population with size (1 - s) * Nref.
    
    Parameters:
        T: Time in the past of split (in units of 2*Na generations).
        s: Fraction of Nref that founds pop 1 (North, island). Pop 1 size = s. Pop 2 size = 1-s.
    """
    T, s = params
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T, nu1=s, nu2=1-s, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


#Model 2: vicariance asymetric migration
def vic_asym_mig(params, ns, pts):
    """
    Split into two populations with continuous asymmetric migration.
    
    Assumes:
        - North = Pop1 = founder (island) population, smaller size.
        - South = Pop2 = mainland population, larger size.
    
    Parameters:
        T: Time since split (in units of 2*Na generations).
        s: Fraction of Nref that founds the North (Pop1).
        m12: Migration rate from South (Pop2) to North (Pop1).
        m21: Migration rate from North (Pop1) to South (Pop2).
    """
    T, s, m12, m21 = params
    xx = dadi.Numerics.default_grid(pts)

    # Ancestral population -> split with asymmetric migration
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T, nu1=s, nu2=1 - s, m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


#Model 3: vicariance secondary contact with asymetric migration
def vic_sec_contact_asym_mig(params, ns, pts):
    """
    Initial isolation following a split, then secondary contact with asymmetric migration.
    
    Assumes:
        - North = Pop1 = founder (island) population, smaller size.
        - South = Pop2 = mainland population, larger size.
    
    Parameters:
        m12: Migration rate from South (Pop2) to North (Pop1) during secondary contact.
        m21: Migration rate from North (Pop1) to South (Pop2) during secondary contact.
        T1: Duration of isolation post-split (no migration).
        T2: Duration of secondary contact with asymmetric migration.
        s: Fraction of Nref that founds the North (Pop1).
    """
    m12, m21, T1, T2, s = params
    xx = dadi.Numerics.default_grid(pts)

    # Split: Isolation phase (no migration)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=s, nu2=1 - s, m12=0, m21=0)

    # Secondary contact: Asymmetric migration phase
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=s, nu2=1 - s, m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


def vic_anc_asym_mig(params, ns, pts):
    """
    Ancient asymmetric migration following the split, then isolation.
    
    Assumes:
        - North = Pop1 = founder (island) population, smaller size.
        - South = Pop2 = mainland population, larger size.
    
    Parameters:
        m12: Migration rate from South (Pop2) to North (Pop1) during ancient phase.
        m21: Migration rate from North (Pop1) to South (Pop2) during ancient phase.
        T1: Duration of migration phase after the split.
        T2: Duration of subsequent isolation (no migration).
        s: Fraction of Nref that founds the North (Pop1).
    """
    m12, m21, T1, T2, s = params
    xx = dadi.Numerics.default_grid(pts)

    # Ancestral population → split → ancient migration phase
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=s, nu2=1 - s, m12=m12, m21=m21)

    # Followed by isolation
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=s, nu2=1 - s, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs
##########################################################################################
#Founder event scenarios
###############################################################################

# Model 1: Founder event with no migration
def founder_nomig(params, ns, pts):
    """
    Split into two populations, with no migration. Populations are fractions of the reference 
    population (Nref = 1). Population 1 (North) is the island/founder population with initial 
    size s * Nref, undergoing exponential growth to nu1 * Nref. Population 2 (South) is the 
    mainland population with constant size (1 - s) * Nref.
    
    Parameters:
        nu1: Final size of pop 1 (North), after exponential growth.
        T: Time in the past of split (in units of 2*Na generations).
        s: Fraction of Nref that founds pop 1 (North, island). Pop 1 initial size = s. Pop 2 size = 1-s.
    """
    nu1, T, s = params
    xx = dadi.Numerics.default_grid(pts)
    
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: s * (nu1/s)**(t/T)  # North grows from s to nu1
    
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nu1_func, nu2=1-s, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


# Model 5: Founder event with continuous asymmetric migration
def founder_asym(params, ns, pts):
    """
    Split into two populations, with asymmetric migration. Populations are fractions of the 
    reference population (Nref = 1). Population 1 (North) is the island/founder population 
    with initial size s * Nref, undergoing exponential growth to nu1 * Nref. Population 2 
    (South) is the mainland population with constant size (1 - s) * Nref.
    
    Parameters:
        nu1: Final size of pop 1 (North), after exponential growth.
        m12: Migration rate from South to North (2*Na*m12).
        m21: Migration rate from North to South (2*Na*m21).
        T: Time in the past of split (in units of 2*Na generations).
        s: Fraction of Nref that founds pop 1 (North, island). Pop 1 initial size = s. Pop 2 size = 1-s.
    """
    nu1, m12, m21, T, s = params
    xx = dadi.Numerics.default_grid(pts)
    
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: s * (nu1/s)**(t/T)  # North grows from s to nu1
    
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nu1_func, nu2=1-s, m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

# Model 6: Founder event with continuous asymmetric migration and secondary contact.
def founder_sec_contact_asym(params, ns, pts):
    """
    Split into two populations, with initial isolation followed by asymmetric migration (secondary contact).
    Populations are fractions of the reference population (Nref = 1). Population 1 (North) is the 
    island/founder population with initial size s * Nref, undergoing exponential growth to nu1 * Nref.
    Population 2 (South) is the mainland population with constant size (1 - s) * Nref.
    
    Parameters:
        nu1: Final size of pop 1 (North), after exponential growth.
        m12: Migration rate from South to North (2*Na*m12) during secondary contact.
        m21: Migration rate from North to South (2*Na*m21) during secondary contact.
        T1: Time of isolation post-split (in units of 2*Na generations).
        T2: Time of secondary contact migration until present.
        s: Fraction of Nref that founds pop 1 (North, island). Pop 1 initial size = s. Pop 2 size = 1-s.
    """
    nu1, m12, m21, T1, T2, s = params
    xx = dadi.Numerics.default_grid(pts)
    
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: s * (nu1/s)**(t/(T1+T2))  # North grows from s to nu1 over T1+T2
    
    # Isolation phase (T1) with no migration
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nu1_func, nu2=1-s, m12=0, m21=0)
    
    # Secondary contact phase (T2) with asymmetric migration
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=nu1_func, nu2=1-s, m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs



def founder_anc_asym_mig(params, ns, pts):
    """
    Founder event with ancient asymmetric migration followed by isolation.
    
    Assumes:
        - Pop1 = North = founder (island), smaller population
        - Pop2 = South = mainland population
        - North population grows from s * Nref to nu1 * Nref
    
    Parameters:
        nu1: Final size of pop 1 (North) after growth.
        m12: Migration rate from South (Pop2) to North (Pop1) during ancient phase.
        m21: Migration rate from North (Pop1) to South (Pop2) during ancient phase.
        T1: Time duration of ancient asymmetric migration.
        T2: Time duration of isolation (no migration) until present.
        s: Initial fraction of Nref that founds pop 1 (North).
    """
    nu1, m12, m21, T1, T2, s = params
    xx = dadi.Numerics.default_grid(pts)

    # North grows exponentially from s to nu1 over T1+T2
    nu1_func = lambda t: s * (nu1/s)**(t / (T1 + T2))

    # Begin with ancestral population and split
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Ancient asymmetric migration phase
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nu1_func, nu2=1 - s, m12=m12, m21=m21)

    # Isolation phase
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=nu1_func, nu2=1 - s, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs
#######################################################################
#########################################################################
# -------------------------------
# Step 3: Choose which model to run by setting run_model variable
# Options: "asym_mig", "no_mig", "ancient_mig", "secondary_contact"
run_model = "founder_anc_asym_mig"  # Change this to test different models




#######################################################
#vicariance model settings
###################################################
###################################################################
#Model 1
if run_model == "vic_no_mig":
    upper_bound = [10, 0.5]  # [T, s]
    lower_bound = [1e-2, 1e-4]
    p0 = [0.01842637, 0.27172629]
    func_ex = dadi.Numerics.make_extrap_log_func(vic_no_mig)
    model_name = "vic_no_mig"


elif run_model == "vic_asym_mig":
    upper_bound = [10, 0.5, 10, 10]  # [T, s, m12, m21]
    lower_bound = [1e-2, 1e-4, 0, 0]
    p0 = [0.82506907, 0.24405045, 4.67763869, 9.97423897]
    func_ex = dadi.Numerics.make_extrap_log_func(vic_asym_mig)
    model_name = "vic_asym_mig"


elif run_model == "vic_sec_contact_asym_mig":
    upper_bound = [10, 10, 5, 5, 0.5]  # [m12, m21, T1, T2, s]
    lower_bound = [0, 0, 1e-2, 1e-2, 1e-4]
    p0 = [3.04764208, 9.13191145, 0.02839705, 1.62006253, 0.34687854]
    func_ex = dadi.Numerics.make_extrap_log_func(vic_sec_contact_asym_mig)
    model_name = "vic_sec_contact_asym_mig"


elif run_model == "vic_anc_asym_mig":
    upper_bound = [10, 10, 5, 5, 0.5]  # [m12, m21, T1, T2, s]
    lower_bound = [0, 0, 1e-2, 1e-2, 1e-4]
    p0 = [0.39100985, 0.13993596, 0.01034213, 0.01083159, 0.33229878]
    func_ex = dadi.Numerics.make_extrap_log_func(vic_anc_asym_mig)
    model_name = "vic_anc_asym_mig"
################################################################



#######################################################
#founder model settings
###################################################
###############################################################
#Model 1
elif run_model == "founder_nomig":
    upper_bound = [10, 10, 0.5]  # [nu1, T, s]
    lower_bound = [1e-2, 1e-2, 1e-4]
    p0 = [7.38029186, 0.01707085, 0.03648225]
    func_ex = dadi.Numerics.make_extrap_log_func(founder_nomig)
    model_name = "founder_nomig"


#Model 4
elif run_model == "founder_asym":
    upper_bound = [10, 20, 20, 5, 0.5]  # [nu1, m12, m21, T, s]
    lower_bound = [1e-2, 0, 0, 1e-2, 1e-4]
    p0 = [0.5, 1, 1, 0.1, 0.05]  # Initialize m12, m21 with founder_sym’s m
    func_ex = dadi.Numerics.make_extrap_log_func(founder_asym)
    model_name = "founder_asym"



#Model 5
elif run_model == "founder_sec_contact_asym":
    upper_bound = [10, 20, 20, 5, 5, 0.5]  # [nu1, m12, m21, T1, T2, s]
    lower_bound = [1e-2, 0, 0, 1e-2, 1e-2, 1e-4]
    p0 = [3.74363751e+00, 8.04906822e+00, 9.27213551e+00, 2.93843604e-01, 4.88056078e-01, 8.06552985e-03]  # Split T=0.46 into T1=T2
    func_ex = dadi.Numerics.make_extrap_log_func(founder_sec_contact_asym)
    model_name = "founder_sec_contact_asym"



elif run_model == "founder_anc_asym_mig":
    upper_bound = [10, 20, 20, 5, 5, 0.5]  # [nu1, m12, m21, T1, T2, s]
    lower_bound = [1e-2, 0, 0, 1e-2, 1e-2, 1e-4]
    p0 = [0.5, 1.0, 1.0, 0.3, 0.1, 0.05]
    func_ex = dadi.Numerics.make_extrap_log_func(founder_anc_asym_mig)
    model_name = "founder_anc_asym_mig"
################################################################





else:
    raise ValueError("Invalid value for run_model. Use 'asym_mig', 'no_mig', 'ancient_mig', or 'secondary_contact'.")

# Perturb initial parameters
p0 = dadi.Misc.perturb_params(p0, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)
print(f"Model: {model_name}, Perturbed initial parameters: {p0}")

# Define extrapolation grid points
pts_l = [170, 180, 190]  # Adjusted to reflect projection size

# -------------------------------
# Step 4: Prepare organized output directory structure
run_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")  # Unique folder for each run

# Base directory
base_output_dir = "/Users/michaeltofflemire/Library/CloudStorage/Dropbox/storage/storage/local/projects/chamaea_genomics/output/dadi-models/runs_2025_06_30_16-00-pub-version"

# Create model directory if it doesn’t exist
model_output_dir = os.path.join(base_output_dir, model_name)
os.makedirs(model_output_dir, exist_ok=True)

# Create unique run directory
run_output_dir = os.path.join(model_output_dir, f"run_{run_timestamp}")
os.makedirs(run_output_dir, exist_ok=True)

# File paths for saving outputs
params_path   = os.path.join(run_output_dir, "parameters.txt")
logs_path     = os.path.join(run_output_dir, "logs.txt")
model_sfs_path = os.path.join(run_output_dir, "model_sfs.fs")
metrics_path  = os.path.join(run_output_dir, "metrics.txt")
plot_path     = os.path.join(run_output_dir, "observed_vs_model.png")
run_info_path = os.path.join(run_output_dir, "run_info.txt")

# Master results file (Tracks all dadi runs)
#master_results_path = os.path.join(base_output_dir, "master_results.txt")
# New location (per-model master file)
master_results_path = os.path.join(model_output_dir, "master_results.txt")

# -------------------------------
# Step 5: Save run configuration details
with open(run_info_path, "w") as f:
    f.write("Run Information:\n")
    f.write(f"Model Name: {model_name}\n")
    f.write(f"Run Timestamp: {run_timestamp}\n")
    f.write(f"Initial Parameters: {p0}\n")
    f.write(f"Upper Bound: {upper_bound}\n")
    f.write(f"Lower Bound: {lower_bound}\n")
    f.write(f"Extrapolation Points: {pts_l}\n")
    f.write(f"Optimization Settings: maxiter=10\n")
    f.write(f"Number of SNPs: {num_snps}\n\n")

# -------------------------------
# Step 6: Redirect stdout to log file and run optimization
print('Beginning optimization ************************************************')
with open(logs_path, "w") as log_file:
    original_stdout = sys.stdout
    sys.stdout = log_file
    try:
        popt = dadi.Inference.optimize_log(
            p0, data, func_ex, pts_l, 
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            verbose=10, maxiter=10
        )
    finally:
        sys.stdout = original_stdout

print('Finished optimization ************************************************')
print("Optimized parameters:", popt)

# Save optimized parameters
with open(params_path, "w") as f:
    f.write(f"Initial parameters: {p0}\nOptimized parameters: {popt}")

# -------------------------------
# Step 7: Compare the Model SFS with the Observed SFS
model_sfs = func_ex(popt, ns, pts_l)

# Calculate log-likelihood and optimal theta
ll_model = dadi.Inference.ll_multinom(model_sfs, data)
theta = dadi.Inference.optimal_sfs_scaling(model_sfs, data)
print(f"Log-likelihood of the model: {ll_model}")
print(f"Optimal value of theta: {theta}")

# Save model SFS and metrics
model_sfs.to_file(model_sfs_path)
with open(metrics_path, "w") as f:
    f.write(f"Log-likelihood: {ll_model}\nOptimal theta: {theta}")

# -------------------------------
# Step 8: Append results to master results file
with open(master_results_path, "a") as summary_file:
    summary_file.write(f"Run Timestamp: {run_timestamp}\n")
    summary_file.write(f"Model: {model_name}\n")
    summary_file.write(f"Optimized Parameters: {popt}\n")
    summary_file.write(f"Log-likelihood: {ll_model}\n")
    summary_file.write(f"Optimal Theta: {theta}\n")
    summary_file.write(f"Number of SNPs: {num_snps}\n")
    summary_file.write(f"Projection Sizes: {pts_l}\n")
    summary_file.write("-" * 60 + "\n")

# -------------------------------
# Step 9: Plot and save the observed vs. model SFS
plt.figure(figsize=(8, 6))
dadi.Plotting.plot_2d_comp_multinom(model_sfs, data, vmin=1, resid_range=3, pop_ids=('North', 'South'))
plt.savefig(plot_path)
plt.show()
