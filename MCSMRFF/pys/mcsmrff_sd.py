import os
import mcsmrff_files, mcsmrff_utils, mcsmrff_gradient

# A function for steepest descent optimization of parameters
def steepest_descent(run_name, alpha=0.05, maxiter=1000, gtol=1E-3, perturbation=1.01, param_file=None, three_body=None, tersoff=None, lj_coul=None): #better, but tends to push error up eventually, especially towards endpoints.
	if param_file is None:
		# It will look for an input file of type "input_runname.tersoff" by default
		parameters = read_params(run_name)
	else:
		parameters = read_params(param_file, exact=True)
	parameters = list(parameters)

	atoms, systems_by_composition = get_training_set(run_name, use_pickle=True)	

	print("\n\nStep        Avg Energy        Avg rms_force        error_force (%)        error_energy (%)\
             \n------------------------------------------------------------------------------------------")

	# Get current Energy and rms_force for the given parameters
	mcsmrff_gradient.run_lammps(atoms, systems_by_composition, parameters[0], parameters[1], parameters[2], run_name)
	energies, rms_forces = mcsmrff_files.parse_lammps_output(run_name, len(atoms.atoms))

	nLJ = len(np.array(parameters[0]).flatten())

	step = 0
	while (step < maxiter):
		# Get gradient and error of the system with given parameters
		gradient = mcsmrff_gradient.get_gradient(list(parameters), atoms, systems_by_composition, run_name, perturbation=perturbation, three_body=three_body)
		error_force, error_energy = mcsmrff_gradient.calculate_error(run_name, len(atoms.atoms))
		error_force *= 100.0
		error_energy *= 100.0

		# Print output
		print("%d            %.2f            %.2f            %.2f            %.2f" % (step, sum(energies)/len(energies), sum(rms_forces)/len(rms_forces), error_force, error_energy))

		# Store parameters used in new array for perterbation later
		f = np.append(np.array(parameters[0]).flatten().copy(), np.array(parameters[2]).copy())
		
		# Find step to take
		max_step_length = np.sqrt((np.square(gradient, dtype=np.float64)).max(),dtype=np.float64)

		# This is to deal with scenarios in which an overflow is encountered
		if np.isinf(max_step_length):
			max_step_length = 1000.0

		# Generate step to take
		if max_step_length > 1.0:
			dr = gradient * alpha / max_step_length
		else:
			dr = gradient * alpha

		# Calculate new parameters
		# Note, we subtract the step
		f -= np.array(dr).flatten()
		parameters[0] = f[:nLJ].reshape((-1,2))
		p_hold = np.array(parameters[2]).copy()
		parameters[2] = f[nLJ:]

		# Prevent beta from being negative
		for i,p in enumerate(parameters[2]):
			if abs(p) < 1e-9: parameters[2][i] = 0
			if i%14 == 7 and p < 0:
				parameters[2][i] = 1e-8

		# If we improve by changing m (dr[i] < 0 st i%14 ==0) then flip it
		# Recall, dr is the change in error. Thus a negative one implies the new error is smaller than the old one
		for i,p in enumerate(p_hold):
			if i % 14 == 0:
				if dr[i+nLJ] < 0:
					parameters[2][i] = int(1 if np.round(p) == 3 else 3)
				else:
					parameters[2][i] = int(np.round(p))

		mcsmrff_gradient.run_lammps(atoms, systems_by_composition, parameters[0], parameters[1], parameters[2], run_name)

		# Get current energies and rms_force for the given parameters
		energies, rms_forces = mcsmrff_files.parse_lammps_output(run_name, len(atoms.atoms))
		# Save this iteration of parameters
		if not os.path.isdir("parameters"): os.mkdir("parameters")
		os.chdir("parameters")
		mcsmrff_files.write_params(parameters[0],parameters[1],parameters[2],run_name,append="_output")
		os.chdir("../")

		step += 1

	# Save the final parameters
	os.chdir("parameters")
	mcsmrff_files.write_params(parameters[0],parameters[1],parameters[2],run_name,append="_output")
	os.chdir("../")

	return parameters