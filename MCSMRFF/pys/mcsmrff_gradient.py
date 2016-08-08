#This file contains a few functions for use when making the gradient optimization methods for MCSMRFF
#It is now fully set up to run one LAMMPS simulation

import os, sys
import cPickle as pickle
import numpy as np
import copy

from merlin import *
from re import findall
from mcsmrff_constants import *
from mcsmrff_utils import *
from mcsmrff_files import *

# This function will read in all orca simulations from a folder "training_sets" in the current working directory
# If use_pickle is True, then it will first check to see if a pickle file exists in the training_sets folder under the name "training_set.pickle"
#     If file does not exist - Read in all data from subfolders and generate the file
#     If file does exist - Read in the pickle file if it exists
# If pickle_file_name is a string, then it will read in the pickle file specified by the string path
def get_training_set(run_name, use_pickle=True, pickle_file_name=None):
	# Take care of pickle file I/O
	# Get file name
	if pickle_file_name is None:
		pfile = "training_sets/training_set.pickle"
	else:
		pfile = pickle_file_name

	system = None
	# If the pickle file does not exist, then make it
	# If use_pickle is False, then make the read in the data from the training_sets folder
	if not os.path.isfile(pfile) or not use_pickle:
		if pickle_file_name is not None:
			raise Exception("Requested file %s, but unable to read it in." % pickle_file_name)

		# Generate the pickle itself if it doesn't exist
		# Create the size of the box to be 1000 x 100 x 100 to hold your training sets
		system = utils.System(box_size=[1e3, 100.0, 100.0], name="training_set") 
		systems_by_composition = {}

		# For each folder in the training_sets folder lets get the cml file we want and write the energies and forces for that file
		for name in os.listdir("training_sets"):
			# We'll read in any training subset that succeeded and print a warning on those that failed
			try:
				result = orca.read("training_sets/%s/%s.out" % (name, name))
			except IOError:
				print("Warning - Training Subset %s not included as results not found..." % name)
				continue

			# Parse the force output and change units. In the case of no force found, do not use this set of data
			try:
				forces = orca.engrad_read("training_sets/%s/%s.orca.engrad" % (name, name), pos="Ang")[0]
				# Convert force from Ha/Bohr to kcal/mol-Ang
				convert = lambda x: units.convert_dist("Ang","Bohr",units.convert_energy("Ha","kcal",x))
				for a,b in zip(result.atoms, forces): a.fx, a.fy, a.fz = convert(b.fx), convert(b.fy), convert(b.fz)
			except (IndexError, IOError) as e:
				print("Warning - Training Subset %s not included as results not found..." % name)
				continue
				#for a in result.atoms: a.fx, a.fy, a.fz = 0.0, 0.0, 0.0

			# Get the bonding information
			with_bonds = utils.Molecule("training_sets/%s/system.cml" % name, extra_parameters=extra_Pb, test_charges=False)

			# Copy over the forces read in into the system that has the bonding information
			for a,b in zip(with_bonds.atoms, result.atoms):
				a.fx, a.fy, a.fz = b.fx, b.fy, b.fz
				if utils.dist(a,b)>1e-4: raise Exception('Atoms are different:', (a.x,a.y,a.z), (b.x,b.y,b.z)) # sanity check on atom positions

			# Rename some things
			with_bonds.energy = result.energy
			with_bonds.name = name

			# Now, we read in all the potential three-body interactions that our training set takes into account
			# This will be in a 1D array
			composition = ' '.join(sorted([a.element for a in result.atoms]))
			if composition not in systems_by_composition:
				systems_by_composition[composition] = []
			systems_by_composition[composition].append(with_bonds)

		# Generate (1) xyz file of various systems as different time steps and (2) system to simulate
		xyz_atoms = []
		for composition in systems_by_composition:
			# Sort so that the lowest energy training subset is first in the system
			systems_by_composition[composition].sort(key=lambda s:s.energy)
			baseline_energy = systems_by_composition[composition][0].energy
			# Offset the energies by the lowest energy, convert units of the energy
			for s in systems_by_composition[composition]:
				s.energy -= baseline_energy
				s.energy = units.convert_energy("Ha","kcal/mol",s.energy)
				# Don't use high-energy systems, because these will not likely be sampled in MD
				if s.energy > 500.0: continue
				# For testing purposes, output
				print "DEBUG:", s.name, s.energy
				xyz_atoms.append(s.atoms)
				system.add(s, len(system.molecules)*1000.0)

		# Make the box just a little bigger (100) so that we can fit all our systems
		system.xhi = len(system.molecules)*1000.0+100.0

		# Write all of the states we are using to training_sets.xyz
		if not os.path.isdir("training_sets"): os.mkdir("training_sets")
		os.chdir("training_sets")
		files.write_xyz(xyz_atoms, 'training_sets')
		os.chdir("../")
		# Generate our pickle file if desired
		if use_pickle:
			print("Saving pickle file %s..." % pfile)
			fptr = open(pfile, "wb")
			pickle.dump([system, systems_by_composition], fptr)
			fptr.close()

	# If use_pickle is true AND the pickle file exists, then we can just read it in
	if system is None and use_pickle:
		print("Reading pickle file %s..." % pfile)
		fptr = open(pfile, "rb")
		system, systems_by_composition = pickle.load(fptr)
		system.name = run_name
		fptr.close()
	elif system is None:
		raise Exception("Requested file %s, but unable to read it in." % pfile)

	# Now we have the data, save it to files for this simulation of "run_name" and return parameters
	if not os.path.isdir("lammps"): os.mkdir("lammps") #go to the lammps folder, or make one
	if not os.path.isdir("lammps/%s" % run_name): os.mkdir("lammps/%s" % run_name)
	os.chdir("lammps/%s" % run_name)
	write_system_and_training_data(run_name, system, systems_by_composition)
	os.chdir("../../")

	return system, systems_by_composition

# This function will run a LAMMPS single point calculation using the MCSMRFF force field
def run_lammps(system,systems_by_composition,lj_params,atom_list,tersoff_params,run_name):
	# Ensure cwd/lammps/run_name exists
	directory_to_work = '%s/lammps/%s'%(os.getcwd(),run_name)
	directory_to_work = directory_to_work.split("/")
	for i,d in enumerate(directory_to_work):
		if d.strip() == "": continue
		if not os.path.isdir("/".join(directory_to_work[:i+1])):
			os.mkdir("/".join(directory_to_work[:i+1]))
	directory_to_work = "/".join(directory_to_work)
	os.chdir(directory_to_work)
	system.name = run_name

	# Make a file that has all the names of the systems we are using here
	files.write_lammps_data(system)

	# Begin code for the LAMMPS input file
	#what units are you using - units real
	#atom_style full - type of atoms
	#force field
	#bonds
	#angles 
	#dihedrals
	#no lennard jones bewteen nearest and second nearest neighbors but multiplied by 0.5 by third nearest neighbors
	#
	#fixed bounaries
	#read the data file that we just wrote
	commands = ('''
	units real
	atom_style full
	pair_style hybrid/overlay lj/cut/coul/inout 0.2 3.5 15 tersoff
	bond_style harmonic
	angle_style harmonic
	dihedral_style opls
	special_bonds lj/coul 0.0 0.0 0.5

	boundary f f f
	read_data	'''+run_name+'''.data
	''').splitlines()

	# Indices of OPLS parameters
	tersoff_types = [t for t in system.atom_types if t.index in [Pb,Cl]]
	charges    = lj_params[0]
	lj_sigma   = lj_params[1]
	lj_epsilon = lj_params[2]

	# This is how we are reading in the 104 tersoff parameters (i.e. m,gamma,N,D...etc)
	tersoff_strings, i, j = [], 0, 0
	while i < (len(tersoff_params)-1): # Concatenate the atom names and params
		tmp = atom_list[j:j+3]
		for param in tersoff_params[i:i+14]:
			tmp.append(str(param))
		tersoff_strings.append(tmp)
		i += 14
		j += 3

	# R and D are the distances from the center of the first atom, so this is how we get the cutoff distance between those
	inner_cutoffs = {}
	for type_i in tersoff_types:
		for type_j in tersoff_types:
			for s in tersoff_strings:
				types = s[:3]
				R, D = float(s[13]), float(s[14])
				if types == (type_i.element_name, type_j.element_name, type_j.element_name):
					inner_cutoffs[ (type_i,type_j) ] = R+D

	index = 0
	for t in system.atom_types:
		if t in tersoff_types:
			t.charge = charges[index]
			t.vdw_e = lj_epsilon[index]
			t.vdw_r = lj_sigma[index]
			index += 1

	# Write to the lammps file the atom types and vdw radii
	for i in range(len(system.atom_types)):
		for j in range(i, len(system.atom_types)):
			type_i = system.atom_types[i]
			type_j = system.atom_types[j]
			commands.append('pair_coeff %d %d lj/cut/coul/inout %f %f %f' % (
				i+1, 
				j+1, 
				(type_i.vdw_e*type_j.vdw_e)**0.5, 
				(type_i.vdw_r*type_j.vdw_r)**0.5, 
				inner_cutoffs[ (i,j) ] if (i,j) in inner_cutoffs else 0.0) 
			)
		commands.append('set type %d charge %f' % (i+1, type_i.charge) )

	# Generate a lmp object to make the LAMMPS input file
	lmp = utils.Struct()
	lmp.file = open(run_name+'.in', 'w')
	def writeline(line):
		lmp.file.write(line+'\n')
	lmp.command = writeline
	for line in commands:
		lmp.command(line)

	#write the pair_coeff command
	write_params(lj_params,atom_list,tersoff_params,run_name,append="_input")
	write_system_and_training_data(run_name, system, systems_by_composition)
	lmp.command('pair_coeff * * tersoff '+run_name+'_input.tersoff '+(' '.join([ (t.element_name if t in tersoff_types else 'NULL') for t in system.atom_types])))

	for t in system.bond_types:
		lmp.command('bond_coeff %d	%f %f' % (t.lammps_type, t.e, t.r) )
	for t in system.angle_types:
		lmp.command('angle_coeff %d	%f %f' % (t.lammps_type, t.e, t.angle) )
	for t in system.dihedral_types:
		lmp.command('dihedral_coeff %d	%f %f %f %f' % ((t.lammps_type,)+t.e))

	commands = '''
	compute atom_pe all pe/atom
	compute sum_pe all reduce sum c_atom_pe
	neigh_modify once yes

	dump 1 all custom 1 '''+run_name+'''.dump id type x y z fx fy fz 
	run 0
	undump 1
	'''
	for line in commands.splitlines():
		lmp.command(line)

	lmp.file.close()

	# Run the simulation the change directory back to the parent one
	os.system('%s -in %s.in -log %s.log >> out.log' % (lammps_mcsmrff, run_name,run_name))
	os.chdir("../../")

# This function calculates the error between a MCSMRFF simulation of "run_name" and the stored training_forces and training_energies
def calculate_error(run_name, natoms):
	energy, rms_force = parse_lammps_output(run_name, natoms)

	# Now we know the energy and rms_force, just need to get difference from dft results and return
	f_training_forces = open("lammps/%s/%s_training_forces.txt" % (run_name, run_name) ).read().split("\n")
	f_training_energies = open("lammps/%s/%s_training_energies.txt" % (run_name, run_name) ).read().split("\n")
	training_energy = sum([float(x) for x in f_training_energies if x.strip() != ""])
	training_rms_force = np.array([float(x) for x in f_training_forces if x.strip() != ""])
	training_rms_force = np.sqrt((training_rms_force**2).sum() / len(training_rms_force))

	error = (training_rms_force - rms_force) / training_rms_force
  	return abs(error)

# This allows the user to specify which three-body interactions they want to optimize parameters for
def indices_of_desired_three_body(element_strings, three_body):
	element_list = np.array(element_strings).reshape((-1,3))
	if three_body is None:
		index = range(len(element_list))
	else:
		index = []
		for i,s in enumerate(element_list):
			if ",".join(s) in three_body:
				index.append(i)
		if index == []:
			print("Error - No three_body set specified in three_body variable.")
			sys.exit()
	return index

# A function to calculate the gradient of the MCSMRFF parameters given an initial guess and atomic positions
# This will perturb each parameter slightly and build up a gradient
def get_gradient(parameters, system, systems_by_composition, run_name, perturbation=1.01, three_body=None):
	p_lj, p_atoms, p_tersoff = [np.array(p) for p in parameters]
	n_trios = len(p_atoms.reshape((-1,3))) # Find the number of three-body terms
	
	index = indices_of_desired_three_body(parameters[1], three_body)

	atoms = copy.deepcopy(system)
	atoms.name = atoms.name + "_grad_calc"
	
	error_0 = calculate_error(run_name, len(atoms.atoms))
	gradient = np.empty(len(p_tersoff))

	for i,p in enumerate(p_tersoff):
		if int(i / 14) not in index: continue

		perturbed_parameters = p_tersoff.copy()

		if i%14 == 0:
			perturbed_parameters[i] = np.float64(3) if p == 1 else np.float64(1)
		else:
			perturbed_parameters[i] = p * perturbation

		run_lammps(atoms, systems_by_composition, parameters[0], parameters[1], perturbed_parameters, "%s_grad_calc" % run_name)
		
		gradient[i] = np.float64(calculate_error("%s_grad_calc" % run_name, len(atoms.atoms)) - error_0)

	return gradient

# A function for steepest descent optimization of parameters
def steepest_descent(run_name, alpha=0.05, maxiter=1000, gtol=1E-3, perturbation=1.01, param_file=None, three_body=None): #better, but tends to push error up eventually, especially towards endpoints.
	if param_file is None:
		# It will look for an input file of type "input_runname.tersoff" by default
		parameters = read_params(run_name)
	else:
		parameters = read_params(param_file, exact=True)
	parameters = list(parameters)

	atoms, systems_by_composition = get_training_set(run_name, use_pickle=True)	

	print("\n\nStep\tEnergy\t    rms_force    error\n-------------------------------------")
	
	# Get current Energy and rms_force for the given parameters
	run_lammps(atoms, systems_by_composition, parameters[0], parameters[1], parameters[2], run_name)
	energy, rms_gradient = parse_lammps_output(run_name, len(atoms.atoms))

	step = 0
	while (rms_gradient > gtol) and (step < maxiter):
		# Get gradient and error of the system with given parameters
		gradient = get_gradient(list(parameters), atoms, systems_by_composition, run_name, perturbation=perturbation, three_body=three_body)
		error = abs(calculate_error(run_name, len(atoms.atoms))*100.0)

		# Print output
		print("%d\t%.2f\t%.2f\t%.2f" % (step, energy, rms_gradient, error))

		# Store parameters used in new array for perterbation later
		f = np.array(parameters[2]).flatten().copy()

		# Find step to take
		max_step_length = np.sqrt((np.square(gradient, dtype=np.float64)).max(),dtype=np.float64)

		# This is to deal with scenarios in which an overflow is encountered
		if max_step_length == np.inf:
			max_step_length = 1000.0

		# Generate step to take
		if max_step_length > 1.0:
			dr = gradient * alpha / max_step_length
		else:
			dr = gradient * alpha

		# Only change parameters of three-body interactions we specified using "three_body"
		index = indices_of_desired_three_body(parameters[1], three_body)
		for i,p in enumerate(parameters[2]):
			if int(i/14) not in index:
				dr[i] = 0
				f[i] = parameters[2][i]

		# Calculate new parameters
		# Note, we subtract the step
		parameters[2] = f - np.array(dr).flatten()

		for i,p in enumerate(parameters[2]):
			if abs(p) < 1e-9: parameters[2][i] = 0
			if i%14 == 7 and p < 0:
				parameters[2][i] = 1e-8

		# If we improve by changing m (dr[i] < 0 st i%14 ==0) then flip it
		# Recall, dr is the change in error. Thus a negative one implies the new error is smaller than the old one
		for i,p in enumerate(parameters[2]):
			if i % 14 == 0:
				p = f[i]
				if dr[i] < 0:
					parameters[2][i] = int(1 if p == 3 else 3)
				else:
					parameters[2][i] = int(p)

		run_lammps(atoms, systems_by_composition, parameters[0], parameters[1], parameters[2], run_name)

		# Get current Energy and rms_force for the given parameters
		energy, rms_gradient = parse_lammps_output(run_name, len(atoms.atoms))
		# Save this iteration of parameters
		if not os.path.isdir("parameters"): os.mkdir("parameters")
		os.chdir("parameters")
		write_params(parameters[0],parameters[1],parameters[2],run_name,append="_output")
		os.chdir("../")

		step += 1

	# Save the final parameters
	os.chdir("parameters")
	write_params(parameters[0],parameters[1],parameters[2],run_name,append="_output")
	os.chdir("../")

	return parameters

# A function to get a 2x2 perovskite crystal for test simulations
def get_test_system():
	L, N = 6.0, 3
	dim = L*N+0.5
	test_system = utils.System(box_size=[dim, dim, dim], name="test_run")
	PbMACl3 = utils.Molecule('systems/unit_cell', extra_parameters=extra_Pb, test_charges=False)

	count = 0
	for xi in range(N):
		for yi in range(N):
			for zi in range(N):
				count += 1
				x, y, z = (xi-0.5)*L, (yi-0.5)*L, (zi-0.5)*L
				test_system.add(PbMACl3, x, y, z)

	return test_system
