# Imports from system
import os
import copy
import shutil
import itertools
import numpy as np
import cPickle as pickle
# Imports from Clancelot
import orca
import files
import structures
import geometry
import units
# Imports from MCSMRFF
import mcsmrff_files
from mcsmrff_constants import *

# A function to read in the seed folder
def read_seed(path="./seed"):
	"""
	Read in all cml files from the seed directory.

	**Parameters**

		path: *str, optional*
			A path to the seed directory.

	**Returns**

		molecules_A: *list, list, molecules*
			A list of molecules from the seed directory
		molecules_B: *list, molecules*
			A list of molecules from the seed directory.  In this case, we merge child molecules into one.
	"""
	if path.endswith("/"): path = path[:-1]
	if not os.path.exists(path):
		raise Exception("Unable to find seed directory")

	molecules_A = []
	for fptr in os.listdir(path):
		if not fptr.endswith(".cml"): continue
		molecules_A.append(files.read_cml(path+"/"+fptr, return_molecules=True, allow_errors=True, test_charges=False))

	if molecules_A == []:
		raise Exception("Seed directory is empty")

	molecules_B = []
	for seed in molecules_A:
		atoms, bonds, angles, dihedrals = [], [], [], []
		for mol in seed:
			atoms += mol.atoms
			bonds += mol.bonds
			angles += mol.angles
			dihedrals += mol.dihedrals
		molecules_B.append(structures.Molecule(atoms, bonds, angles, dihedrals))

	return molecules_A, molecules_B

def compile_training_set(path="./training_sets"):
	if path.endswith("/"): path = path[:-1]
	if not os.path.exists(path):
		raise Exception("Unable to find training set directory")

	frames = []
	i = 0
	while os.path.exists(path+"/"+str(i)+".cml"):
		atoms = files.read_cml(path+"/"+str(i)+".cml", return_molecules=False, allow_errors=True, test_charges=False)[0]
		frames.append(atoms)
		i += 1
	files.write_xyz(frames, "full_training_set")

def generate_training_set(N_perterbations_per_seed = 5,
							perterbation_dx = 0.3,
							perterbation_dr = 10,
							expansion_perterbation_dx = 0.1,
							expansion_perterbation_dr = 0.5,
							expansion_step = 0.5,
							N_expansion = 10,
							N_perterbations_per_expansion = 5):
	"""
	From a seed directory of cml files, generate a more robust training set.
	"""
	# Ensure output directory is available
	if os.path.exists("./training_sets"):
		raise Exception("Warning - Would overwrite training set folder already in existance.")
	os.mkdir("training_sets")
	training_set_counter = 0

	# Generate seed list
	seeds_A, seeds_B = read_seed()

	# Add seed list to training set list
	for mol in seeds_B:
		files.write_cml(mol, name="training_sets/%d" % training_set_counter)
		training_set_counter += 1

	# Now, for each seed, add N perterbations
	for mol in seeds_B:
		for i in range(N_perterbations_per_seed):
			m = copy.deepcopy(mol)
			m.perterbate(dx=perterbation_dx, dr=perterbation_dr)
			files.write_cml(m, name="training_sets/%d" % training_set_counter)
			training_set_counter += 1

	# Now, for each seed of multiple molecules, add varying distances
	for mols in seeds_A:
		if len(mols) == 1: continue
		to_run = []

		# Get every pair combination
		for m in range(0, len(mols)+1):
			for subset in itertools.combinations(mols, m):
				if len(subset) == 2:
					to_run.append(copy.deepcopy(subset))
		if to_run == []:
			raise Exception("Unable to get combinations")

		# Loop through combos and expand
		for m1,m2 in to_run:
			axis = np.array(m2.get_center_of_mass()) - np.array(m1.get_center_of_mass())
			axis /= np.linalg.norm(axis)
			step = axis * expansion_step
			for i in range(N_expansion):
				m2.translate(step)
				for j in range(N_perterbations_per_expansion):
					m3 = copy.deepcopy(m1)
					m4 = copy.deepcopy(m2)
					m3.perterbate(dx=expansion_perterbation_dx, dr=expansion_perterbation_dr)
					m4.perterbate(dx=expansion_perterbation_dx, dr=expansion_perterbation_dr)
					files.write_cml([m3,m4], name="training_sets/%d" % training_set_counter)
					training_set_counter += 1

def run_low_level():
	if not os.path.exists("training_sets"):
		raise Exception("No training set folder to run.")

	frange = [int(a.split('.cml')[0]) for a in os.listdir("training_sets") if a.endswith(".cml")]
	frange.sort()
	if len(frange) == 0:
		raise Exception("No viable files in training sets folder.")

	route = "! B97-D3 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} Grid7"
	extra_section = ""

	running_jobs = []
	for i in frange:
		atoms = files.read_cml("training_sets/%d.cml" % i, allow_errors=True, test_charges=False, return_molecules=False)[0]
		charge = sum([a.type.charge for a in atoms])
		running_jobs.append( orca.job("ts_%d" % i, route, atoms=atoms, extra_section=extra_section, charge=charge, grad=True, queue="batch", procs=2, sandbox=True) )
	return running_jobs

def run_high_level():
	if not os.path.exists("training_sets"):
		raise Exception("No training set folder to run.")

	frange = [int(a.split('.cml')[0]) for a in os.listdir("training_sets") if a.endswith(".cml")]
	frange.sort()
	if len(frange) == 0:
		raise Exception("No viable files in training sets folder.")

	route = "! PW6B95 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} Grid7"
	extra_section = ""

	running_jobs = []
	previous_failed = []
	for i in frange:
		atoms = files.read_cml("training_sets/%d.cml" % i, allow_errors=True, test_charges=False, return_molecules=False)[0]
		charge = sum([a.type.charge for a in atoms])
		prev_converged = orca.read("ts_%d" % i).converged
		if prev_converged:
			running_jobs.append( orca.job("ts_%d_high" % i, route, atoms=[], extra_section=extra_section, charge=charge, grad=True, queue="batch", procs=2, previous="ts_%d" % i, sandbox=True) )
		else:
			previous_failed.append(i)
	return running_jobs, previous_failed

# This function will read in all orca simulations from a folder "training_sets" in the current working directory
# If use_pickle is True, then it will first check to see if a pickle file exists in the training_sets folder under the name "training_set.pickle"
#     If file does not exist - Read in all data from subfolders and generate the file
#     If file does exist - Read in the pickle file if it exists
# If pickle_file_name is a string, then it will read in the pickle file specified by the string path
def pickle_training_set(run_name,
						training_sets_folder="training_sets",
						pickle_file_name="training_set",
						high_energy_cutoff=500.0,
						system_x_offset=1000.0,
						verbose=False):
	"""
	A function to picle together the training set in a manner that is readable
	for MCSMRFF.  This is a single LAMMPs data file with each training set
	offset alongst the x-axis by system_x_offset.

	**Parameters**

		run_name: *str*
		training_sets_folder: *str, optional*
			Path to the folder where all the training set data is.
		pickle_file_name: *str, optional*
			A name for the pickle file and training set system.
		high_energy_cutoff: *float, optional*
			A cutoff for systems that are too large in energy, as MD is likely
			never to sample them.
		system_x_offset: *float, optional*
			The x offset for the systems to be added by.
		verbose: *bool, optional*
			Whether to have additional stdout or not.

	**Returns**

		Stuff
	"""
	# Take care of pickle file I/O
	if training_sets_folder.endswith("/"):
		training_sets_folder = training_sets_folder[:-1]
	if pickle_file_name is not None and pickle_file_name.endswith(".pickle"):
		pickle_file_name = pickle_file_name.split(".pickle")[0]
	pfile = training_sets_folder + "/" + pickle_file_name + ".pickle"
	sys_name = pickle_file_name
	if os.path.isfile(pfile):
		raise Exception("Pickled training set already exists!")

	# Generate empty system for your training set
	system = None
	system = structures.System(box_size=[1e3, 100.0, 100.0], name=sys_name) 
	systems_by_composition = {}

	# For each folder in the training_sets folder lets get the cml file we
	# want and write the energies and forces for that file
	for name in os.listdir(training_sets_folder):
		# We'll read in any training subset that succeeded and print a warning
		# on those that failed
		try:
			result = orca.read("%s/%s/%s.out" %
								(training_sets_folder, name, name)
							  )
		except IOError:
			print("Warning - Training Subset %s not included as \
out file not found..." % name)
			continue

		# Check for convergence
		if not result.converged:
			print("Warning - Results for %s have not converged." % name)
			continue

		# Parse the force output and change units. In the case of no force
		# found, do not use this set of data
		try:
			forces = orca.engrad_read("%s/%s/%s.orca.engrad" %
										(training_sets_folder, name, name),
										pos="Ang"
									 )[0]
			# Convert force from Ha/Bohr to kcal/mol-Ang
			convert = lambda x: units.convert_dist("Ang","Bohr",
										units.convert_energy("Ha","kcal",x))
			for a,b in zip(result.atoms, forces):
				a.fx, a.fy, a.fz = convert(b.fx), convert(b.fy), convert(b.fz)
		except (IndexError, IOError) as e:
			print("Warning - Training Subset %s not included as \
results not found..." % name)
			continue

		# Get the bonding information
		with_bonds = structures.Molecule("%s/%s/%s.cml" % 
										(training_sets_folder, name, name),
										extra_parameters=extra_Pb,
										allow_errors=True,
										test_charges=False
									)

		# Copy over the forces read in into the system that has the bonding
		# information
		for a,b in zip(with_bonds.atoms, result.atoms):
			a.fx, a.fy, a.fz = b.fx, b.fy, b.fz
			# sanity check on atom positions
			if geometry.dist(a,b)>1e-4:
				raise Exception('Atoms are different:', (a.x,a.y,a.z),
													    (b.x,b.y,b.z)
								) 

		# Rename and save energy
		with_bonds.energy = result.energy
		with_bonds.name = name

		# Now, we read in all the potential three-body interactions that our
		# training set takes into account.  This will be in a 1D array
		composition = ' '.join(sorted([a.element for a in result.atoms]))
		if composition not in systems_by_composition:
			systems_by_composition[composition] = []
		systems_by_composition[composition].append(with_bonds)

	# Generate:
	#  (1) xyz file of various systems as different time steps
	#  (2) system to simulate
	xyz_atoms = []
	to_delete = []
	for i,composition in enumerate(systems_by_composition):
		# Sort so that the lowest energy training subset is first
		# in the system
		systems_by_composition[composition].sort(key=lambda s:s.energy)
		baseline_energy = systems_by_composition[composition][0].energy
		# Offset the energies by the lowest energy, and convert energy units
		for j,s in enumerate(systems_by_composition[composition]):
			s.energy -= baseline_energy
			s.energy = units.convert_energy("Ha","kcal/mol",s.energy)
			# Don't use high-energy systems, because these will not likely
			# be sampled in MD
			if s.energy > high_energy_cutoff:
				to_delete.append([composition,j])
				continue
			# For testing purposes, output
			if verbose:
				print "Using:", s.name, s.energy
			xyz_atoms.append(s.atoms)
			system.add(s, len(system.molecules)*system_x_offset)
	
	# Delete the system_names that we aren't actually using due to energy
	# being too high
	to_delete = sorted(to_delete, key=lambda x: x[1])[::-1]
	for d1,d2 in to_delete:
		if verbose:
			print "Warning - Training Subset %s not included as energy \
is too high..." % systems_by_composition[d1][d2].name
		del systems_by_composition[d1][d2]

	# Make the box just a little bigger (100) so that we can fit all our
	# systems
	system.xhi = len(system.molecules)*system_x_offset+100.0

	# Write all of the states we are using to training_sets.xyz
	files.write_xyz(xyz_atoms, training_sets_folder+'/' + pickle_file_name)
	# Generate our pickle file
	print("Saving pickle file %s..." % pfile)
	fptr = open(pfile, "wb")
	pickle.dump([system, systems_by_composition], fptr)
	fptr.close()

	# Now we have the data, save it to files for this simulation of
	# "run_name" and return parameters
	if not os.path.isdir("lammps"):
		os.mkdir("lammps")
	if not os.path.isdir("lammps/%s" % run_name):
		os.mkdir("lammps/%s" % run_name)
	os.chdir("lammps/%s" % run_name)
	mcsmrff_files.write_system_and_training_data(run_name,
												 system,
												 systems_by_composition
												)
	os.chdir("../../")

	return system, systems_by_composition

def minimize_seeds():
	seeds = []
	seed_names = []
	route_low = "! OPT B97-D3 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} Grid7"
	extra_section_low = ""
	route_high = "! OPT PW6B95 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} Grid7"
	extra_section_high = ""
	for seed in os.listdir("seed"):
		seeds.append(files.read_cml("seed/%s" % seed, allow_errors=True, test_charges=False)[0])
		seed_names.append(seed)
	jobs = []
	for i,seed in enumerate(seeds):
		charge = sum([a.type.charge for a in seed])
		jobs.append( orca.job("seed_%d_low" % i, route_low, atoms=seed, extra_section=extra_section_low, charge=charge, grad=True, queue="batch", procs=2, sandbox=True) )
	for j in jobs: j.wait()
	for i,seed in enumerate(seeds):
		charge = sum([a.type.charge for a in seed])
		jobs.append( orca.job("seed_%d_high" % i, route_high, atoms=[], extra_section=extra_section_high, charge=charge, grad=True, queue="batch", procs=2, previous="seed_%d_low" % i, sandbox=True) )
	for j in jobs: j.wait()
	for i,seed in enumerate(seeds):
		new_pos = orca.read("seed_%d_high" % i).atoms
		if not raw_data.converged:
			print("Failed to optimize %s" % seed_names[i])
			continue
		cml_file = files.read_cml("seed/%s" % seed_names[i], allow_errors=True, test_charges=False, return_molecules=True)
		j=0
		for mol in cml_file:
			for k,a in enumerate(mol.atoms):
				b = new_pos[j]
				a.x, a.y, a.z = b.x, b.y, b.z
				j += 1
		files.write_cml("seed/%s_opt" % seed_names[i], cml_file)

def create_training_set(min_seeds=True):
	if min_seeds:
		minimize_seeds()
	generate_training_set()
	compile_training_set()
	jobs = run_low_level()
	for j in jobs:
		j.wait()
	jobs2, failed_sims = run_high_level()

	print("\nThe following simulations failed to converge at low level:")
	for i in failed_sims:
		print(" %d," % i),

	for j in jobs2:
		j.wait()
		# Copy over cml file
		i = int(j.name.split("_")[1])
		shutil.copyfile("training_sets/%d.cml" % i,
						"orca/%s/%s.cml" % (j.name, j.name)
						)

	if os.path.exists("ts_hybrid"):
		raise Exception("Folder for final training set already exists!")
	os.mkdir("ts_hybrid")

	for fptr in os.listdir("orca"):
		if not os.path.isdir("orca/"+fptr): continue
		if fptr.endswith("_high"):
			os.system("cp -rf orca/%s ts_hybrid/" % fptr)

	pickle_training_set("set1",training_sets_folder="ts_hybrid")

if __name__ == "__main__":
	#create_training_set()
	minimize_seeds()
