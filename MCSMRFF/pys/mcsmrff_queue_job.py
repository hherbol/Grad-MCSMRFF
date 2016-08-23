import os
from utils import pysub

def job(run_name, parameters, path=""):
	
	file_string = '''from mcsmrff_sd import steepest_descent as SD
from mcsmrff_bfgs import bfgs as BFGS
import mcsmrff_run
import mcsmrff_utils

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

new_parameters = BFGS("$RUN_NAME$", step_size=0.1, maxiter=1000, perturbation=1.01, param_file="$PARAMETERS$", three_body=perturbate_these, tersoff=None, lj_coul=None, opt="Force", reset_step_size=5, training_set_file_path="/fs/home/hch54/MCSMRFF/Grad-MCSMRFF/MCSMRFF/training_sets/training_set.pickle")

mcsmrff_run.get_glimpse($RUN_NAME$)

rms, _ = mcsmrff_utils.pdf_metric("glimpse_$RUN_NAME$", lammps_job=True, persist=False, start=0.0, stop=3.0, step=0.01, cutoff=1.5)

print("\\n\\nRMS for run '$RUN_NAME$' is %.2f\\n\\n" % rms)
'''

	for i in range(4):
		file_string = file_string.replace("$RUN_NAME$",run_name)
	if parameters.endswith(".tersoff"):
		file_string = file_string.replace("$PARAMETERS$",parameters)
	else:
		file_string = file_string.replace("$PARAMETERS$",parameters+".tersoff")

	fptr = open(run_name+".py",'w')
	fptr.write(file_string)
	fptr.close()

	path_to_py = os.getcwd()

	if path != "":
		if path.startswith("/"):
			path_to_py = path
		else:
			path_to_py = path_to_py + "/" + path

	if path_to_py.endswith("/"):
		path_to_py = path_to_py[:-1]


	pysub(run_name+".py",nprocs='1',queue="long", path=path_to_py)