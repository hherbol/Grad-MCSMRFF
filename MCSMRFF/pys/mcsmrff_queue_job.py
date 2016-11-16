import os
from job import pysub


def job(run_name, parameters,
        path="", new_pdf_props={}, new_opt_props={}, disregard=[]):

    training_set_pickle_path = (
        "/fs/home/hch54/"
        "Grad-MCSMRFF/MCSMRFF/training_sets/"
        "training_set.pickle")
    pdf_props = {"start": 0.0,
                 "stop": 5.0,
                 "step": 0.01,
                 "cutoff": 10.0,
                 "persist": "False",
                 "quanta": 0.001}
    opt_props = {"step_size": 0.1,
                 "maxiter": 1000,
                 "perturbation": 1.01,
                 "training_set_file_path": training_set_pickle_path}
    for prop in new_pdf_props:
        if prop in pdf_props:
            pdf_props[prop] = new_pdf_props[prop]
        else:
            raise Exception("The pdf parameter %s does not exist." % prop)
    for prop in new_opt_props:
        if prop in opt_props:
            opt_props[prop] = new_opt_props[prop]
        else:
            raise Exception("The optimizer parameter %s does not exist." %
                            prop)

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

new_parameters = BFGS("$RUN_NAME$",
                      step_size=$STEP_SIZE$,
                      maxiter=$MAXITER$,
                      perturbation=$PERTURBATION$,
                      param_file="$PARAMETERS$",
                      three_body=perturbate_these,
                      tersoff=None,
                      lj_coul=None,
                      opt="Force",
                      reset_step_size=5,
                      training_set_file_path=\"$TRAINING_SET_FILE_PATH$\")

mcsmrff_run.get_glimpse("$RUN_NAME$")

rms, _ = mcsmrff_utils.pdf_metric("glimpse_$RUN_NAME$",
                                  lammps_job=True,
                                  persist=$PERSIST$,
                                  start=$START$,
                                  stop=$STOP$,
                                  step=$STEP$,
                                  cutoff=$CUTOFF$,
                                  quanta=$QUANTA$,
                                  disregard=[$DISREGARD$])
print("\\n\\nRMS for run '$RUN_NAME$' is %.5f\\n\\n" % rms)
'''

    for i in range(4):
        file_string = file_string.replace("$RUN_NAME$", run_name)
    if parameters.endswith(".tersoff"):
        file_string = file_string.replace("$PARAMETERS$", parameters)
    else:
        file_string = file_string.replace("$PARAMETERS$",
                                          parameters + ".tersoff")
    for prop in pdf_props:
        file_string = file_string.replace("$%s$" %
                                          prop.upper(), str(pdf_props[prop]))
    for prop in opt_props:
        file_string = file_string.replace("$%s$" %
                                          prop.upper(), str(opt_props[prop]))
    s = ['\"%s\"' % elem for elem in disregard]
    s = ','.join(s)
    file_string = file_string.replace("$DISREGARD$", s)

    fptr = open(run_name + ".py", 'w')
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

    pysub(run_name + ".py", nprocs='1', queue="long", path=path_to_py)
