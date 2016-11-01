# Imports from system
import os
import copy
import itertools
import numpy as np
# Imports from Clancelot
import orca
import files
import structures

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

	route = "! B97-D3 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} Grid7 COSMO"
	extra_section = '''%cosmo
	SMD true
	epsilon 36.7
	end
%geom
	MaxIter 500
	end
'''

	for i in frange:
		atoms = files.read_cml("training_sets/%d.cml" % i, allow_errors=True, test_charges=False)
		orca.job("ts_%d" % i, route, atoms=atoms, extra_section=extra_section, queue="batch", procs=2)

generate_training_set()
compile_training_set()
run_low_level()
