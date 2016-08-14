from mcsmrff_gradient import *
from mcsmrff_run import run as mcsmrff_run
from mcsmrff_bfgs import *
from time import time

run_name = "debug"
parameters = "run_0_output.tersoff"

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

start = time()
new_parameters_1 = steepest_descent(run_name, alpha=0.1, maxiter=10, perturbation=1.01, param_file="parameters/%s" % parameters, three_body=perturbate_these, tersoff=None, lj_coul=None)
stop = time()
print("Time for SD = %.2f seconds" % (stop - start))
start = time()
new_parameters_2 = bfgs(run_name, step_size=0.1, maxiter=10, perturbation=1.01, param_file="parameters/%s" % parameters, three_body=perturbate_these, tersoff=None, lj_coul=None)
stop = time()
print("Time for bfgs = %.2f seconds" % (stop - start))
