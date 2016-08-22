'''
mcsmrff_opt
-----------

Optimization algorithm for parameterization of MCSMRFF.  Works as follows:

1. Generate M sets of N dimensional parameters using LHC.  That is, there is a
   set of parameters we need for MCSMRFF to run which is N dimensional, and we
   have M sets of these to search the energy landscape.

2. Submit M simulations in parallel on the cluster.  These will run either bfgs
   or steepest descent (SD) to minimize the parameter space.

3. As simulations begin converging, compare converged parameters P_i with those
   converged in other clusters.  In the case of
		error(P_i) > error(P_j)
   kill all jobs in neighbourhood of P_j and expand P_i neighbourhood.
'''

import numpy as np
import os

from mcsmrff_bfgs import bfgs as BFGS
from mcsmrff_gradient import steepest_descent as SD
from mcsmrff_lhs import create_lhs
from mcsmrff_queue_job import job
from mcsmrff_files import read_params, write_params

## Constants
NUMBER_OF_SETS = 100
EXPANSION_NUMBER = 5
SIMULATION_NAME = "Large_Run"

# Part 1, generation of parameter sets
PARAMETER_SETS = create_lhs(num_samples=NUMBER_OF_SETS)
P1,P2,_ = read_params("parameters/run_0_output.tersoff",exact=True)

# Part 2, submit jobs to queue
if not os.path.isdir(SIMULATION_NAME):
	os.mkdir(SIMULATION_NAME)
os.chdir(SIMULATION_NAME)
if not os.path.isdir("parameters"):
	os.mkdir("parameters")

for i,P3 in enumerate(PARAMETER_SETS):
	job_name = "%s_%d" % (SIMULATION_NAME, i)
	pname = "parameters/" + job_name
	write_params(P1, P2, P3, pname, append="")
	job(job_name, pname)

# Part 3, monitor jobs. Follow lowest error jobs. Whenever one finishes with
#         higher error, kill neighbourhood and resubmit EXPANSION_NUMBER more
