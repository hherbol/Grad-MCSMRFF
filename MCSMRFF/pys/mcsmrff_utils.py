import os, sys
import numpy as np
import cPickle as pickle
import copy

import files
import debyer
import structures
import geometry
import mcsmrff_files
from mcsmrff_constants import *

def is_float(x):
	try:
		float(x)
		return True
	except ValueError:
		return False

def read_dump_worker(fptr, data=["element","x","y","z"]):
	ids = None
	atoms_flag = False

	frames, frame = [], []
	for line in open(fptr,'r'):
		if line.startswith("ITEM: ATOMS"):
			atoms_flag = True
			frame = []
			if ids is None:
				ids = {name:index for index,name in enumerate(line.split()[2:])}
				continue
		elif line.startswith("ITEM:"):
			atoms_flag = False
			if frame != []: frames.append(frame)
			frame = []
			continue
		elif atoms_flag:
			# Read in the data
			atom = []
			values = line.split()
			for d in data:
				try:
					atom.append( np.float64(values[ids[d]]) )
				except:
					atom.append( values[ids[d]] )
					print values
			frame.append(atom)

	return frames

def frames_to_xyz(filename, frames):
	fptr = open(filename,'w')
	for frame in frames:
		fptr.write("%d\nAtoms\n" % len(frame))
		for atom in frame:
			fptr.write("%s\t%lg\t%lg\t%lg\n" % tuple(atom))
	fptr.close()

def read_dump(fptr, persist=False):
	frames = read_dump_worker(fptr, data=["type","xu","yu","zu"])

	swap = {1:"Pb",2:"Cl",3:"N",4:"C",5:"H",6:"H"}

	for i,frame in enumerate(frames):
		for j,atom in enumerate(frame):
			frames[i][j][0] = swap[atom[0]]

	frames_to_xyz("tmp_read_dump.xyz", frames)
	frames = files.read_xyz("tmp_read_dump.xyz")

	if not persist:
		os.system("rm tmp_read_dump.xyz")

	return frames

# Given frames (xyz file format read in using files.read_xyz), get metric
def pdf_metric(A, ref=None, persist=False, lammps_job=False, start=0.0, stop=10.0, step=0.1, cutoff=3.0, quanta=0.001, disregard=[]):
	# If we are checking a lammps job, grab the xyz file from lammps/run_name/run_name.xyz
	if lammps_job==True:
		if A.endswith(".xyz"): A.split(".xyz")[0]
		if A.endswith(".dump"): A.split(".dump")[0]
		if A.endswith(".data"): A.split(".data")[0]
		A = read_dump("lammps/%s/%s2.dump" % (A,A))

	# If we passed a string, then read in the file
	if type(A) is str:
		if not A.endswith(".xyz"): A += ".xyz"
		A = files.read_xyz(A)

	# Assume reference is the first frame
	if ref is None:
		ref = A[0]
	else:
		raise Exception("This is not coded yet.")

	B = copy.deepcopy([A[0],A[-1]])
	# Remove anything in disregard
	for i,frame in enumerate(B):
		to_kill = []
		for j,atom in enumerate(frame):
			if atom.element in disregard:
				to_kill.append(j)
		to_kill = sorted(to_kill)[::-1]
		for k in to_kill:
			del B[i][k]

	pdf_ref = debyer.get_pdf(B[0], persist=persist, output="tmp_pdf_ref", start=start, stop=stop, step=step, cutoff=cutoff, quanta=quanta)
	pdf_final = debyer.get_pdf(B[-1], persist=persist, output="tmp_pdf_final", start=start, stop=stop, step=step, cutoff=cutoff, quanta=quanta)
	files.write_xyz([B[0],B[-1]], "pdf_metric_debug")

	# Split lists
	pdf_ref = zip(*pdf_ref)
	pdf_final = zip(*pdf_final)

	difference = [(a-b)**2 for a,b in zip(pdf_ref[1], pdf_final[1])]
	rms = (np.asarray(difference).sum() / float(len(difference)))**0.5

	if not persist:
		os.system("rm pdf_metric_debug.xyz")

	return rms, [pdf_ref, pdf_final]

# A function to get a 2x2 perovskite crystal for test simulations
def get_test_system(length_in_ang=6.0, number_per_side=3, path_to_unit_cell='/fs/home/hch54/MCSMRFF/Grad-MCSMRFF/MCSMRFF/systems/unit_cell'):
	L = length_in_ang
	N = number_per_side
	dim = L*N+0.5
	test_system = structures.System(box_size=[dim, dim, dim], name="test_run")
	PbMACl3 = structures.Molecule(path_to_unit_cell, extra_parameters=extra_Pb, test_charges=False)

	count = 0
	for xi in range(N):
		for yi in range(N):
			for zi in range(N):
				count += 1
				x, y, z = (xi-0.5)*L, (yi-0.5)*L, (zi-0.5)*L
				test_system.add(PbMACl3, x, y, z)

	return test_system

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
		system = structures.System(box_size=[1e3, 100.0, 100.0], name="training_set") 
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
			with_bonds = structures.Molecule("training_sets/%s/system.cml" % name, extra_parameters=extra_Pb, test_charges=False)

			# Copy over the forces read in into the system that has the bonding information
			for a,b in zip(with_bonds.atoms, result.atoms):
				a.fx, a.fy, a.fz = b.fx, b.fy, b.fz
				if geometry.dist(a,b)>1e-4: raise Exception('Atoms are different:', (a.x,a.y,a.z), (b.x,b.y,b.z)) # sanity check on atom positions

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
		to_delete = []
		for i,composition in enumerate(systems_by_composition):
			# Sort so that the lowest energy training subset is first in the system
			systems_by_composition[composition].sort(key=lambda s:s.energy)
			baseline_energy = systems_by_composition[composition][0].energy
			# Offset the energies by the lowest energy, convert units of the energy
			for j,s in enumerate(systems_by_composition[composition]):
				s.energy -= baseline_energy
				s.energy = units.convert_energy("Ha","kcal/mol",s.energy)
				# Don't use high-energy systems, because these will not likely be sampled in MD
				if s.energy > 500.0:
					to_delete.append([composition,j])
					continue
				# For testing purposes, output
				print "DEBUG:", s.name, s.energy
				xyz_atoms.append(s.atoms)
				system.add(s, len(system.molecules)*1000.0)
		
		# Delete the system_names that we aren't actually using due to energy being too high
		to_delete = sorted(to_delete, key=lambda x: x[1])[::-1]
		for d1,d2 in to_delete:
			print "Warning - Training Subset %s not included as energy is too high..." % systems_by_composition[d1][d2].name
			del systems_by_composition[d1][d2]

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
	mcsmrff_files.write_system_and_training_data(run_name, system, systems_by_composition)
	os.chdir("../../")

	return system, systems_by_composition