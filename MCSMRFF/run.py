from functions_for_gradient import *
from run_mcsmrff import run as run_mcsmrff

# Generate a large unit_cell for a test system
test_system = get_test_system()
#test_system.name = "test_start"

# Run MCSMRFF with initial parameters
#old_parameters = read_params("test_o", exact=True)
#old_parameters = read_params("test")
old_parameters = read_params("reasonable.tersoff", exact=True)

run_mcsmrff("test_start", test_system, old_parameters, RUN="200000")

# Optimize the parameters
perturbate_these = [
					"Pb,Pb,Pb",
					"Pb,Pb,Cl",
					"Pb,Cl,Pb",
					"Pb,Cl,Cl",
					"Cl,Pb,Pb",
					"Cl,Pb,Cl",
					"Cl,Cl,Pb",
					"Cl,Cl,Cl"
				]

new_parameters = steepest_descent("test", alpha=0.2, maxiter=5000, perturbation=1.01, param_file="reasonable.tersoff", three_body=perturbate_these)

# See how the updates did
run_mcsmrff("test_end", test_system, new_parameters, RUN="200000")