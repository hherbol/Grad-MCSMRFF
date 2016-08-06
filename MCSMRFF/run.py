from mcsmrff_gradient import *
from mcsmrff_run import run as mcsmrff_run

run_name = "test"
run_leng = "300000"

# Generate a large unit_cell for a test system
test_system = get_test_system()
test_system.name = run_name

# Run MCSMRFF with initial parameters ONLY IF IT HASN'T BEEN ALREADY DONE!
old_parameters = read_params("parameters/reasonable.tersoff", exact=True)
if not os.path.isdir("lammps/reasonable"):
	mcsmrff_run("reasonable", test_system, old_parameters, RUN=run_leng)

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

new_parameters = steepest_descent(run_name, alpha=0.001, maxiter=10000, perturbation=1.01, param_file="parameters/reasonable.tersoff", three_body=perturbate_these)
#new_parameters = steepest_descent(run_name, alpha=0.05, maxiter=20, perturbation=1.01, param_file="parameters/reasonable.tersoff", three_body=None)

# See how the updates did
mcsmrff_run("%s_SD" % run_name, test_system, new_parameters, RUN=run_leng)
