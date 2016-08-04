import os, sys
import cPickle as pickle
import numpy as np
import copy

from merlin import *
from re import findall

#This file contains a few functions for use when making the gradient optimization methods for MCSMRFF
#It is now fully set up to run one LAMMPS simulation

def read_params(run_name): #read in a file and output a list of all the TERSOFF parameters for a file of the format shown below in the order shown below
	filename = "input_%s.tersoff"%run_name 
	#####################################################################################################
	# Generated by min_style params, run = 'june15'. LAMMPS units = real.
	# Error metric=1.00602
	# Charges:           0.4       -0.2 
	# LJ-sigma:            3          2 
	# LJ-epsilon:        0.1        0.1 
	# i, j, k,       m,        gamma,    lambda3,  c,        d,        costheta0,
	#                n,        beta,     lambda2,  B,        R,        D,        lambda1,  A
	#Pb Pb Pb         3         1   1.33856   4240.71    1.0451  0.248252
        #         1 1.4621e-07   2.77144   14721.8   3.00116 0.0721795    4.8965      2142

	#Pb Pb Cl         3         1   2.35497   991.032  0.783726  0.415971
	#                 1 1.4621e-07   2.77144   14721.8   2.81598  0.137382    4.8965      2142

	#Pb Cl Pb         3         1  0.142447   624.254  0.757654 0.000200989
	#                 1 1.4621e-07   2.55336   38004.7  0.483911 0.0368801   3.44825   8876.94

	#Pb Cl Cl         3         1 1.44822e-05   1743.67  0.398607 -0.439513
	#                 1 0.00019596   1.43378   7254.15         3   0.18053    3.2259     50000

	#Cl Pb Pb         3         1 6.41602e-05   49806.5  10 -0.740097
	#                 1      0.01       1.4    3916.4    3.3541 0.0104918   2.49099      5000

	#Cl Pb Cl         3         1 3.87838e-06   764.695        300         0
	#                 1 0.00841335   2.55336   38004.7   3.61493 0.0004688   3.44825   8876.94

	#Cl Cl Pb         3         1         5     800        300  0.389188
	#                 1 0.00841335   2.33527   98110.6   3.18962 0.00270539         2     36788

	#Cl Cl Cl         3         1   1.91467    100000   12.1841       0.5
	#                 1 0.00841335   2.33527   98110.6   4.30146  0.187761         1.9     36788
	#####################################################################################################
	lj_paramlist = []
	for line in open(filename): #this gives you the lj_sigma and epsilon you dont need this but its here just in case
		if line.startswith('# Charges:'): charges = [float(x) for x in line.split()[2:]]
		if line.startswith('# LJ-sigma:'): lj_sigma = [float(x) for x in line.split()[2:]]
		if line.startswith('# LJ-epsilon:'): lj_epsilon = [float(x) for x in line.split()[2:]]
	lj_paramlist.append(charges)
	lj_paramlist.append(lj_sigma)
	lj_paramlist.append(lj_epsilon)

	
	tersoff_paramlist = []
	atom_list = []
	#this is how we are reading in the 112 tersoff parameters INCLUDING m...DONT MINIMIZE m (i.e. m,gamma,N,D...etc) only minimize the actual 104 parameters
	tersoff_strings = findall('\n'+('(\S+) +'*9)[:-2]+' *\n +'+('(\S+) +'*8)[:-2], open(filename).read()) #re is regex library which is a way of pseudo-code to parse strings
	for atom_triplet in tersoff_strings:
		for param in atom_triplet[3:]:
			tersoff_paramlist.append(float(param))
		for atom in atom_triplet[:3]:
			atom_list.append(atom) #list of all atoms (Pb Pb Cl, Pb Cl Cl, etc)
	
	return lj_paramlist, atom_list, tersoff_paramlist #returns a list of all the parameters for both LJ and Tersoff as well as the atoms used

def write_params(lj_params,atom_list,tersoff_params,filename,append="_o"): #read in a list of LJ and tersoff in order as shown above in read_params and write this to a file
	ofile = open(filename+append+".tersoff","w") #write a file run_name_o.tersoff as the output
	ofile.write('''# Generated by min_style params, run = 'None'. LAMMPS units = real.
# Error metric=0.0
# Charges:           '''+str(lj_params[0][0])+'''       '''+str(lj_params[0][1])+''' 
# LJ-sigma:          '''+str(lj_params[1][0])+'''       '''+str(lj_params[1][1])+'''
# LJ-epsilon:        '''+str(lj_params[2][0])+'''       '''+str(lj_params[2][1])+''' 
# i, j, k,       m,        gamma,    lambda3,  c,        d,        costheta0,
#                n,        beta,     lambda2,  B,        R,        D,        lambda1,  A\n''')
	i = 0
	j = 0
	while i < (len(tersoff_params)-1): #write the output of the tersoff params, atoms in the list, etc.
		ofile.write(write_atom_line(atom_list,j)+write_param_line(tersoff_params,i))
		i += 14
		j += 3
	return

def write_param_line(tersoff,i): #write a single line of a tersoff input file (given the tersoff list and the index of the line. primarily written in conjunction with write_params
	return ("%f   %f   %f   %f   %f   %f\n\t\t %f   %f   %f   %f   %f   %f   %f   %f\n\n" % (tersoff[i],tersoff[i+1],tersoff[i+2],tersoff[i+3],tersoff[i+4],tersoff[i+5],tersoff[i+6],tersoff[i+7],tersoff[i+8],tersoff[i+9],tersoff[i+10],tersoff[i+11],tersoff[i+12],tersoff[i+13]))

def write_atom_line(atom,i): #write a single line of atoms (given the atom list and the index of the line. primarily written in conjunction with write_params
	return ("%s %s %s\t" % (atom[i],atom[i+1],atom[i+2]))

def get_training_set(run_name, use_pickle=True, pickle_file_name=None): #initialize the box and get the training set data from /training_set/ and return a System object with this information
	if not os.path.isdir("lammps"): os.mkdir("lammps") #go to the lammps folder, or make one
	if not os.path.isdir("lammps/%s" % run_name): os.mkdir("lammps/%s" % run_name)
	os.chdir("lammps/%s" % run_name)

	if pickle_file_name is not None:
		if ".pickle" not in pickle_file_name:
			pickle_file_name += ".pickle"
	else:
		pickle_file_name = run_name + ".pickle"

	if use_pickle and os.path.isfile(pickle_file_name):
		print("\nReading pickle file %s" % pickle_file_name)
		fptr = open(pickle_file_name, "rb")
		system, systems_by_composition = pickle.load(fptr)
		system.name = run_name
		fptr.close()

		files.write_lammps_data(system) #make a file that has all the names of the systems we are using here

		f = open(system.name+'_training_data.txt', 'w')
		for composition in systems_by_composition:
			for s in systems_by_composition[composition]:
				f.write(s.name+'\n')
		f.close()

		# write forces to a file
		f = open(system.name+'_training_forces.txt', 'w') #DFT forces
		for a in system.atoms:
			f.write("%e\n%e\n%e\n" % (a.fx, a.fy, a.fz) )
		f.close()

		# write energies to a file 
		f = open(system.name+'_training_energies.txt', 'w') #DFT energies these are what we are comparing against and calculating energy error from
		for m in system.molecules:
			f.write("%e\n" % (m.energy) )
		f.close()
		os.chdir("../../")

		return system
	os.chdir("../../")

	I_ = 66 
	Cl_ = 21
	H_ = 54
	N_ = 53
	Pb_ = 111

	Pb = 907
	I = 838
	Cl = 344
	HN = 233

	extra = {
		Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.4, vdw_e=10.1, vdw_r=3.0),
	}

	#create the size of the box to be 1000 x 100 x 100 to hold your training sets
	system = utils.System(box_size=[1e3, 100.0, 100.0], name=run_name) 
	systems_by_composition = {}

	#for each folder in the training_sets folder lets get the cml file we want and write the energies and forces for that file
	for outer in [os.getcwd()]:
		directories = next(os.walk(outer+'/training_sets'))[1] #go into the training_sets dir
		for name in directories:
			try:
				result = orca.read(outer+'/training_sets/'+name+'/'+name+'.out')
			except IOError:
				continue
			try:
				forces = orca.engrad_read(outer+'/training_sets/'+name+'/'+name+'.orca.engrad', pos='Ang')[0] #read the engrad file
				convert = lambda x: units.convert_dist("Ang","Bohr",units.convert_energy("Ha","kcal",x)) #convert the forces
				for a,b in zip(result.atoms, forces):
					a.fx, a.fy, a.fz = convert(b.fx), convert(b.fy), convert(b.fz)
			except (IndexError, IOError) as e: # no forces available, so use blanks
				print "WARNING RUNNING WITH FORCES = 0.0" #FIXME Henry fix this so it can find the forces from the MP2 perturbation calculations
				print name
				for a in result.atoms:
					a.fx, a.fy, a.fz = 0.0, 0.0, 0.0
		
			with_bonds = utils.Molecule(outer+'/training_sets/'+name+'/system.cml', extra_parameters=extra, test_charges=False)
		
			for a,b in zip(with_bonds.atoms, result.atoms):
				a.fx, a.fy, a.fz = b.fx, b.fy, b.fz # copy forces
				if utils.dist(a,b)>1e-4: raise Exception('Atoms are different:', (a.x,a.y,a.z), (b.x,b.y,b.z)) # sanity check on atom positions
		
			with_bonds.energy = result.energy
			with_bonds.name = name
		
			composition = ' '.join(sorted([a.element for a in result.atoms]))
			if composition not in systems_by_composition:
				systems_by_composition[composition] = []
			systems_by_composition[composition].append(with_bonds) #adding all the atom types to the system. systems_by_composition is list of all systems you want by training set

	xyz_atoms = []

	for composition in systems_by_composition: #within each type of system, lowest energy must be first and equal to 0.0
		systems_by_composition[composition].sort(key=lambda s:s.energy) #sort so that the lowest energy is first
		baseline_energy = systems_by_composition[composition][0].energy
		print composition
		for s in systems_by_composition[composition]:
			s.energy -= baseline_energy
			s.energy = units.convert_energy("Ha","kcal/mol",s.energy)
			if s.energy > 500.0: continue #don't use high-energy systems, because these will not likely be sampled in MD
			print s.name, s.energy #for testing purposes
			xyz_atoms.append(s.atoms) #for testing purposes
			system.add(s, len(system.molecules)*1000.0)

	system.xhi = len(system.molecules)*1000.0+100.0 #make the box just a little bigger (100) so that we can fit all our systems

	files.write_xyz(xyz_atoms, 'states') #write all of the states we are using to states.xyz in this dir
	
	os.chdir("lammps/%s" % run_name)
	
	#write trainingset names to a file _data.txt
	files.write_lammps_data(system) #make a file that has all the names of the systems we are using here
	f = open(system.name+'_training_data.txt', 'w')
	for composition in systems_by_composition:
		for s in systems_by_composition[composition]:
			f.write(s.name+'\n')
	f.close()
	# write forces to a file
	f = open(system.name+'_training_forces.txt', 'w') #DFT forces
	for a in system.atoms:
		f.write("%e\n%e\n%e\n" % (a.fx, a.fy, a.fz) )
	f.close()
	# write energies to a file 
	f = open(system.name+'_training_energies.txt', 'w') #DFT energies these are what we are comparing against and calculating energy error from
	for m in system.molecules:
		f.write("%e\n" % (m.energy) )
	f.close()

	if use_pickle and not os.path.isfile(pickle_file_name):
		print("\nSaving pickle file %s" % pickle_file_name)
		fptr = open(pickle_file_name, "wb")
		pickle.dump([system, systems_by_composition], fptr)
		fptr.close()

	os.chdir("../../")
	return system

def run_lammps(system,lj_params,atom_list,tersoff_params,run_name): #read in the training set (system) LJ params [], tersoff params [],and run name
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
	files.write_lammps_data(system) #make a file that has all the names of the systems we are using here

	I_ = 66 
	Cl_ = 21
	H_ = 54
	N_ = 53
	Pb_ = 111

	Pb = 907
	I = 838
	Cl = 344
	HN = 233
	commands = ('''units real
atom_style full
pair_style hybrid/overlay lj/cut/coul/inout 0.2 3.5 15 tersoff
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5

boundary f f f
read_data	'''+run_name+'''.data
	''').splitlines() 
	###header for the lammps input script as shown above###
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


	tersoff_types = [t for t in system.atom_types if t.index in [Pb,Cl]]  #indices of OPLS parameters (used to have a HN in this list)
	charges    = lj_params[0]
	lj_sigma   = lj_params[1]
	lj_epsilon = lj_params[2]

	#this is how we are reading in the 104 tersoff parameters (i.e. m,gamma,N,D...etc)
	tersoff_strings = []
	i = 0
	j = 0
	while i < (len(tersoff_params)-1): #concatenate the atom names and params
		tmp = atom_list[j:j+3]
		for param in tersoff_params[i:i+14]:
			tmp.append(str(param))
		tersoff_strings.append(tmp)
		i += 14
		j += 3

	#R and D are the distances from the center of the first atom, so this is how we get the cutoff distance between those
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
	#write to the lammps file the atoms types and the vdw radii for each
	for i in range(len(system.atom_types)):
		for j in range(i, len(system.atom_types)):
			type_i = system.atom_types[i]
			type_j = system.atom_types[j]
			commands.append('pair_coeff %d %d lj/cut/coul/inout %f %f %f' % (i+1, j+1, (type_i.vdw_e*type_j.vdw_e)**0.5, (type_i.vdw_r*type_j.vdw_r)**0.5, inner_cutoffs[ (i,j) ] if (i,j) in inner_cutoffs else 0.0) )
		commands.append('set type %d charge %f' % (i+1, type_i.charge) )

	lmp = utils.Struct()
	lmp.file = open(run_name+'.in', 'w')
	def writeline(line):
		lmp.file.write(line+'\n')
	lmp.command = writeline
	for line in commands:
		lmp.command(line)

	#write the pair_coeff command
	write_params(lj_params,atom_list,tersoff_params,run_name,append="_input")
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
	os.system('/fs/home/afh72/lammps/lammps-7Dec15/src/lmp_serial -in %s.in -log %s.log >> out.log' % (run_name,run_name))
	os.chdir("../../")

def parse_lammps_output(run_name, natoms):
	#grad_run_name = run_name + "_grad_calc"
	grad_run_name = run_name
	# Parse output
	energy = lammps_job.read("%s" % grad_run_name, trj_file=None)[0].data[0][-2]

	## Read in forces
	force_file = open("lammps/%s/%s.dump" % (grad_run_name, grad_run_name)).read().split("\n")
	i = 0
	while i < len(force_file) and "ITEM: ATOMS" not in force_file[i]:
		i += 1

	ids = {name:index for index,name in enumerate(force_file[i].split()[2:])}
	# Now i points to the first line of data from the dump file. Start reading in
	i += 1
	j = 0

	forces = np.empty(natoms*3)
	while i < len(force_file):
		line = force_file[i].strip()
		if line == "":
			i += 1
			continue
		else:
			values = line.split()
			forces[j] = values[ids["fx"]]
			forces[j+1] = values[ids["fy"]]
			forces[j+2] = values[ids["fz"]]
			j += 3
			i += 1
	rms_force = np.sqrt((forces**2).sum() / len(forces))
	
	return energy, rms_force

def calculate_error(run_name, natoms):
	energy, rms_force = parse_lammps_output(run_name, natoms)

	# Now we know the energy and rms_force, just need to get difference from dft results and return
	f_training_forces = open("lammps/%s/%s_training_forces.txt" % (run_name, run_name) ).read().split("\n")
	f_training_energies = open("lammps/%s/%s_training_energies.txt" % (run_name, run_name) ).read().split("\n")
	training_energy = sum([float(x) for x in f_training_energies if x.strip() != ""])
	training_rms_force = np.array([float(x) for x in f_training_forces if x.strip() != ""])
	training_rms_force = np.sqrt((training_rms_force**2).sum() / len(training_rms_force))

	error = (training_rms_force - rms_force) / training_rms_force

  	return error

# A function to calculate the gradient of the MCSMRFF parameters given an initial guess and atomic positions
# This will perturb each parameter slightly and build up a gradient
def get_gradient(parameters, system, run_name, perturbation=1.01):
	p_lj, p_atoms, p_tersoff = [np.array(p) for p in parameters]
	n_trios = len(p_atoms.reshape((-1,3))) # Find the number of three-body terms
	
	atoms = copy.deepcopy(system)
	atoms.name = atoms.name + "_grad_calc"
	
	error_0 = calculate_error(run_name, len(atoms.atoms))
	gradient = np.empty(len(p_tersoff))
	for i,p in enumerate(p_tersoff):
		if i%14 == 0:
			gradient[i] = 0
			continue
		perturbed_parameters = p_tersoff.copy()
		perturbed_parameters[i] = p * perturbation

		run_lammps(atoms, parameters[0], parameters[1], perturbed_parameters, "%s_grad_calc" % run_name)
		
		gradient[i] = calculate_error(run_name, len(atoms.atoms)) - error_0

	return gradient

# A function for steepest descent optimization of parameters
def steepest_descent(run_name, alpha=0.05, maxiter=1000, gtol=1E-3): #better, but tends to push error up eventually, especially towards endpoints.
	parameters = read_params(run_name) #it will look for an input file of type "input_runname.tersoff"
	parameters = list(parameters)
	atoms = get_training_set(run_name, use_pickle=True, pickle_file_name=run_name)	

	step = 0
	print("Step\tEnergy\t    rms_force    error\n-----------------------------------")
	# Get current Energy and rms_force for the given parameters
	run_lammps(atoms, parameters[0], parameters[1], parameters[2], run_name)
	energy, rms_gradient = parse_lammps_output(run_name, len(atoms.atoms))
	while (rms_gradient > gtol) and (step < maxiter):
		# Get gradient of the system
		gradient = get_gradient(list(parameters), atoms, run_name)
		error = calculate_error(run_name, len(atoms.atoms))*100.0
		if error < 0: error *= -1
		print("%d\t%.2f\t%.2f\t%.2f" % (step, energy, rms_gradient, error))

		# Calculate new parameters
		f = np.array(parameters[2]).flatten().copy()

		max_step_length = np.sqrt(((f)**2).sum(axis=0).max())
		if max_step_length > 1.0:
			dr = f * alpha / max_step_length
		else:
			dr = f * alpha

		parameters[2] = f + np.array(dr).flatten()
		for i,p in enumerate(parameters[2]):
			if i % 14 == 0: parameters[2][i] = int(p)

		step += 1
		run_lammps(atoms, parameters[0], parameters[1], parameters[2], run_name)
		# Get current Energy and rms_force for the given parameters
		energy, rms_gradient = parse_lammps_output(run_name, len(atoms.atoms))
	

	write_params(parameters[0],parameters[1],parameters[2],run_name,append="_o")


steepest_descent("test", alpha=1.0, maxiter=5)

#THIS CODE WILL NOW RUN 1 LAMMPS SIMULATION
#parameters = read_params("test") #it will look for an input file of type "input_runname.tersoff"
#atoms = get_training_set("test", use_pickle=True, pickle_file_name="test")
#
#grad = get_gradient(parameters, atoms, "test", perturbation=1.01)
#f = open("Gradient1","w")
#f.write(str(grad))
#f.close()
#
#grad = get_gradient(parameters, atoms, "test", perturbation=1.05)
#f = open("Gradient2","w")
#f.write(str(grad))
#f.close()
#
#grad = get_gradient(parameters, atoms, "test", perturbation=1.1)
#f = open("Gradient3","w")
#f.write(str(grad))
#f.close()
#
#grad = get_gradient(parameters, atoms, "test", perturbation=1.2)
#f = open("Gradient4","w")
#f.write(str(grad))
#f.close()

#run_lammps(atoms,parameters[0],parameters[1],parameters[2],"test")
