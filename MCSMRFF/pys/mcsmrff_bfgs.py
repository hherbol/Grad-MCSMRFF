import os
import sys
import numpy as np
import copy

import mcsmrff_files
import mcsmrff_utils
import mcsmrff_gradient

# A function for steepest descent optimization of parameters
def bfgs(run_name, step_size=0.05, step_size_adjustment=0.5, maxiter=1000, gtol=1E-3, perturbation=1.01,
	     param_file=None, three_body=None, tersoff=None, lj_coul=None, opt = "Force",
		 linesearch="armijo", armijio_line_search_factor=1E-4, reset_step_size=10, max_step_size=0.2, display=0, callback=None,
		 training_set_file_path=None):

	def rebuild_params(params, N):
		parameters[0] = params[:N].reshape((-1,2))
		parameters[2] = params[N:]
		return parameters

	def apply_restrictions(params, p_hold, dr, nLJ):
		# Prevent beta from being negative
		for i,p in enumerate(params[2]):
			if abs(p) < 1e-9: params[2][i] = 0
			if i%14 == 7 and p < 0:
				params[2][i] = 1e-8

		# If we improve by changing m (dr[i] < 0 st i%14 ==0) then flip it
		# Recall, dr is the change in error. Thus a negative one implies the new error is smaller than the old one
		for i,p in enumerate(p_hold):
			if i % 14 == 0:
				if dr[i+nLJ] < 0:
					params[2][i] = int(1 if np.round(p) == 3 else 3)
				else:
					params[2][i] = int(np.round(p))
		return params

	# It will look for an input file of type "input_runname.tersoff" by default
	if param_file is None:
		parameters = mcsmrff_files.read_params(run_name)
	else:
		parameters = mcsmrff_files.read_params(param_file, exact=True)
	parameters = list(parameters)

	# Remove old old, and move old to old old
	if os.path.isdir("parameters/%s" % run_name):
		if os.path.isdir("parameters/OLD_%s" % run_name):
			os.system("rm -rf parameters/OLD_%s" % run_name)
		os.rename("parameters/%s" % run_name, "parameters/OLD_%s" % run_name)

	atoms, systems_by_composition = mcsmrff_utils.get_training_set(run_name, use_pickle=True, pickle_file_name=training_set_file_path)	

	print("\n\nStep        Avg Energy        Avg rms_force        error_force (%)        error_energy (%)\
			 \n------------------------------------------------------------------------------------------")
	sys.stdout.flush()

	# Initialize inv Hess and Identity matrix
	current_parameters = copy.deepcopy(parameters)
	nLJ = len(np.array(current_parameters[0]).flatten())
	I = np.eye(nLJ+len(current_parameters[2]), dtype=np.float64)
	current_Hessian = I.copy()

	# Get gradient and store your old func_max
	mcsmrff_gradient.run_lammps(atoms, systems_by_composition, current_parameters[0], current_parameters[1], current_parameters[2], run_name)
	energies, rms_forces = mcsmrff_files.parse_lammps_output(run_name, len(atoms.atoms))
	E_avg = sum(energies)/len(energies)
	F_avg = sum(rms_forces)/len(rms_forces)

	# Hold original values
	ALPHA_CONST = step_size
	BETA_CONST = step_size_adjustment
	RESET_CONST = reset_step_size
	MIN_STEP = 1E-8
	BACKTRACK_EPS = 1E-3
	
	# Get function to describe linesearch
	if linesearch is 'armijo':
		if display > 1: print("armijo linesearch "),
		def check_backtrack(f1,f0,gk,pk,armijio_line_search_factor,step_size):
			return f1-f0 > armijio_line_search_factor*step_size*np.dot(gk,pk)
	else:
		if display > 1: print("default linesearch "),
		def check_backtrack(f1,f0,pk,gk,armijio_line_search_factor,step_size):
			return (f1-f0)/(abs(f1)+abs(f0)) > BACKTRACK_EPS

	loop_counter = 0
	while (loop_counter < maxiter):
		# Get gradient and error of the system with given parameters
		# Here we drop ljparams and tersoff_params into a 1D array and get the gradient
		working_parameters = np.append(np.array(current_parameters[0]).flatten().copy(), np.array(current_parameters[2]).copy())
		current_gradient = mcsmrff_gradient.get_gradient(list(current_parameters), atoms, systems_by_composition, run_name, perturbation=perturbation, three_body=three_body)

		# Scale the gradient by the largest value. This is because our "gradient" is relatively arbitrary and we want full control by the specified
		# step_size and step_size_adjustment variables
		current_gradient /= max(abs(current_gradient))
		
		error_force, error_energy = mcsmrff_gradient.calculate_error(run_name, len(atoms.atoms))
		error_force *= 100.0
		error_energy *= 100.0

		# We'll focus on minimizing Average RMS Force
		if opt == "Force":
			old_fval = error_force
		elif opt == "Energy":
			old_fval = error_energy
		else:
			amount_F = float(opt)
			old_fval = amount_F*error_force + (1.0-amount_F)*error_energy

		# Print output
		print("%d            %.2f            %.2f            %.2f            %.2f" % (loop_counter, sum(energies)/len(energies), sum(rms_forces)/len(rms_forces), error_force, error_energy))
		sys.stdout.flush()
			
		# Find step to take
		step_direction = -np.dot(current_Hessian, current_gradient)
		force_mags = current_gradient**2
		
		# Prevent division by 0 in situations that we have 0/0 because no parameter is changed
		for i,f in enumerate(force_mags):
			if abs(f) == abs(step_direction[i]) and abs(f) == 0:
				step_direction[i] = 1

		# Remove the magnitude imparted by the hessian through rescaling
		scalar = np.sqrt(force_mags / (step_direction**2))
		step_direction = (step_direction.T * scalar).T
		
		# If we are doing unreasonably small step sizes, quit
		step_lengths = np.sqrt(step_direction**2) * step_size
		if max(step_lengths) < MIN_STEP:
			if display > 0: print("Error - Step size unreasonable (%lg)" % abs(max(step_length))),
			warnflag = 2
			break

		# If we have too large of a step size, set to max
		# TODO: Maybe scale everything so we remain on the eigenvector
		indices_of_large_steps = [(i,s) for i,s in enumerate(step_lengths) if s > max_step_size]
		for i,s in indices_of_large_steps: step_direction[i] *= max_step_size/s
		max_step_flag = len(indices_of_large_steps) > 0
		step_direction = step_direction.flatten()

		# As we are changing values manually, this is no longer
		# the BFGS(Hess) algorithm so reset the Inverse Hessian
		# -> This is because we're no longer on the eigendirection
		if max_step_flag:
			if display > 0: print("Warning - Setting step to max step size"),
			current_Hessian = I.copy()

		# Hold new parameters
		# We do so because we may not want to keep them if it makes the next step bad.
		new_working_parameters = working_parameters + step_size * step_direction

		new_parameters = rebuild_params(new_working_parameters, nLJ)

		new_parameters = apply_restrictions(new_parameters, current_parameters[2], step_direction, nLJ)

		# Get the new gradient and check if max has increased
		mcsmrff_gradient.run_lammps(atoms, systems_by_composition, new_parameters[0], new_parameters[1], new_parameters[2], run_name)
		new_gradient = mcsmrff_gradient.get_gradient(list(new_parameters), atoms, systems_by_composition, run_name, perturbation=perturbation, three_body=three_body)
		# Scale the gradient by the largest value. This is because our "gradient" is relatively arbitrary and we want full control by the specified
		# step_size and step_size_adjustment variables
		new_gradient /= max(abs(new_gradient))
		new_energies, new_rms_forces = mcsmrff_files.parse_lammps_output(run_name, len(atoms.atoms))
		new_E_avg = sum(new_energies)/len(new_energies)
		new_F_avg = sum(new_rms_forces)/len(new_rms_forces)

		error_force, error_energy = mcsmrff_gradient.calculate_error(run_name, len(atoms.atoms))
		error_force *= 100.0
		error_energy *= 100.0

		# We'll focus on minimizing Average RMS Force
		if opt == "Force":
			fval = error_force
		elif opt == "Energy":
			fval = error_energy
		else:
			amount_F = float(opt)
			fval = amount_F*error_force + (1.0-amount_F)*error_energy

		if check_backtrack(fval, old_fval, new_gradient, step_direction, armijio_line_search_factor, step_size):
			# Step taken overstepped the minimum.  Lowering step size
			if display > 0:
				print("\tResetting System as %lg > %lg!" % (fval, old_fval))
				print("\talpha: %lg" % step_size),

			step_size *= np.float64(step_size_adjustment)

			if display > 0: print("-> %lg\n" % step_size)

			# Reset the Inverse Hessian if desired.
			# It is still up for debate if this is to be recommended or not.  As the 
			# inverse hessian corects itself, it might not be important to do this.
			current_Hessian = I
			reset_step_size = RESET_CONST
			continue

		# This allows for the edge case in which after decreasing step_size, a situation arises
		# in which larger alphas are acceptable again. Thus, we reset to the original step_size
		else:
			reset_step_size -= 1
			# If we want to reset_step_size and step_size has been decreased before, set to initial vals
			if reset_step_size < 0 and step_size < ALPHA_CONST:
				if display > 1: print("\tResetting Alpha, Beta, Reset and Inverse Hessian")
				step_size, step_size_adjustment, reset_step_size = ALPHA_CONST, BETA_CONST, RESET_CONST
				# Once again, debatable if we want this here.  When reseting step sizes we
				# might have a better H inverse than the Identity would be.
				current_Hessian = I
				continue
		
		# Recalculate change_in_parameters to maintain the secant condition
		change_in_parameters = new_working_parameters - working_parameters
		
		# Store new max value in old_max for future comparison
		old_fval = fval

		# Get difference in gradients for further calculations
		change_in_gradient = new_gradient - current_gradient

		try:  # this was handled in numeric, let it remaines for more safety
			rhok = 1.0 / (np.dot(change_in_gradient, change_in_parameters))
		except ZeroDivisionError:
			rhok = 1000.0
			if display > 1:
				print("Divide-by-zero encountered: rhok assumed large")
		if np.isinf(rhok):  # this is patch for np
			rhok = 1000.0
			if display > 1:
				print("Divide-by-zero encountered: rhok assumed large")


		# Run BFGS Update for the Inverse Hessian
		A1 = I - change_in_parameters[:, np.newaxis] * change_in_gradient[np.newaxis, :] * rhok
		A2 = I - change_in_gradient[:, np.newaxis] * change_in_parameters[np.newaxis, :] * rhok
		current_Hessian = np.dot(A1, np.dot(current_Hessian, A2)) + (rhok * change_in_parameters[:, np.newaxis] * change_in_parameters[np.newaxis, :])

		if display > 1: print("fval %lg" % (fval))
		# Increment the loop counter
		loop_counter += 1

		# Store new parameters, as it has passed the check
		current_parameters = new_parameters
		current_gradient = new_gradient

		if not os.path.isdir("parameters"): os.mkdir("parameters")
		if not os.path.isdir("parameters/%s" % run_name): os.mkdir("parameters/%s" % run_name)
		os.chdir("parameters/%s" % run_name)
		mcsmrff_files.write_params(parameters[0],parameters[1],parameters[2],run_name,append="_%d" % loop_counter)
		os.chdir("../../")

		# If callback is desired
		if callback is not None:
			callback(current_parameters)


		mcsmrff_gradient.run_lammps(atoms, systems_by_composition, current_parameters[0], current_parameters[1], current_parameters[2], run_name)

		# Get current energies and rms_force for the given parameters
		energies, rms_forces = mcsmrff_files.parse_lammps_output(run_name, len(atoms.atoms))
		# Save this iteration of parameters
		if not os.path.isdir("parameters"): os.mkdir("parameters")
		os.chdir("parameters")
		mcsmrff_files.write_params(parameters[0],parameters[1],parameters[2],run_name,append="_output")
		os.chdir("../")

	# Save the final parameters
	os.chdir("parameters")
	mcsmrff_files.write_params(parameters[0],parameters[1],parameters[2],run_name,append="_output")
	os.chdir("../")

	return parameters
