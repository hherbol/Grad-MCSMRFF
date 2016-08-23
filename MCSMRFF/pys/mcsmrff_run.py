# A function to run a MCSMRFF simulation using the given (1) run_name, (2) system, and (3) parameter file
# Currently it is highly tuned to Pb, Cl, I systems
# Note, parameters must be a list of [LJ, atom_strs, tersoff]
from merlin import *
import shutil

import files
from sysconst import lammps_mcsmrff as LAMMPS_DIR
from hashlib import md5
from re import findall

import mcsmrff_files, mcsmrff_utils
from mcsmrff_constants import *

def run_mcsmrff(run_name, system, parameters, seed=None):

	print("\n\n\n")
	for i in range(5):
		for j in range(50): print("#"),
		print("")
	print("\n\n\n")

	if seed is None:
		seed = str(int(md5(run_name).hexdigest(), 16)%(2**16))
	else:
		seed = str(seed)
	system.name = run_name

	# Change to the correct working directory
	if not os.path.isdir("lammps"): os.mkdir("lammps")
	if not os.path.isdir("lammps/%s" % run_name): os.mkdir("lammps/%s" % run_name)
	os.chdir("lammps/%s" % run_name)

	# Generate files
	files.write_xyz(system.atoms)
	files.write_lammps_data(system)
	mcsmrff_files.write_params(parameters[0], parameters[1], parameters[2], run_name, append="_input")

	# Begin generating string to hold LAMMPS input
	commands = ('''units real
	atom_style full
	pair_style hybrid/overlay lj/cut/coul/inout 0.2 0.0 12 tersoff
	bond_style harmonic
	angle_style harmonic
	dihedral_style opls
	special_bonds lj/coul 0.0 0.0 0.5

	boundary p p p
	read_data	'''+run_name+'''.data
	''').splitlines()

	# Removed HN for no H3 input types
	tersoff_types = [t for t in system.atom_types if t.index in [Pb,Cl]]

	# Grab the parameters. Could in theory just use "parameters" variable, but too lazy to change code
	print os.getcwd()
	for line in open('%s_input.tersoff' % run_name):
		if line.startswith('# Charges:'): charges = [float(x) for x in line.split()[2:]]
		if line.startswith('# LJ-sigma:'): lj_sigma = [float(x) for x in line.split()[2:]]
		if line.startswith('# LJ-epsilon:'): lj_epsilon = [float(x) for x in line.split()[2:]]

	tersoff_strings = findall('\n'+('(\S+) +'*9)[:-2]+' *\n +'+('(\S+) +'*8)[:-2], open(run_name+'_input.tersoff').read())
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

	# Write out parameters to LAMMPS input string
	for i in range(len(system.atom_types)):
		for j in range(i, len(system.atom_types)):
			type_i = system.atom_types[i]
			type_j = system.atom_types[j]
			commands.append('pair_coeff %d %d lj/cut/coul/inout %f %f %f' % (i+1, j+1, (type_i.vdw_e*type_j.vdw_e)**0.5, (type_i.vdw_r*type_j.vdw_r)**0.5, inner_cutoffs[ (i,j) ] if (i,j) in inner_cutoffs else 0.0) )
		commands.append('set type %d charge %f' % (i+1, type_i.charge) )

	# Generate a LAMMPS structure for an input file
	lmp = utils.Struct()
	lmp.file = open(run_name+'.in', 'w')
	def writeline(line):
		lmp.file.write(line+'\n')
	lmp.command = writeline
	for line in commands:
		lmp.command(line)

	lmp.command('pair_coeff * * tersoff '+run_name+'_input.tersoff '+(' '.join([ (t.element_name if t in tersoff_types else 'NULL') for t in system.atom_types])))

	for t in system.bond_types:
		lmp.command('bond_coeff %d	%f %f' % (t.lammps_type, t.e, t.r) )
	for t in system.angle_types:
		lmp.command('angle_coeff %d	%f %f' % (t.lammps_type, t.e, t.angle) )
	for t in system.dihedral_types:
		lmp.command('dihedral_coeff %d	%f %f %f %f' % ((t.lammps_type,)+t.e))

	commands = '''
	neigh_modify every 1 check yes delay 0
	dump	1 all xyz 100 '''+run_name+'''.xyz
	dump    2 all custom 100 '''+run_name+'''2.dump type xu yu zu
	thermo_style custom step temp press ke pe epair emol vol
	thermo 1000

	group mobile id > 0

	velocity all create 10.0 '''+seed+''' rot yes dist gaussian
	timestep 0.01

	fix press all npt temp 10.0 10.0 10.0 iso 0.0 0.0 100.0
	run 50000
	#run 100000
	unfix press

	fix motion all nvt temp 10.0 300.0 100.0
	run 200000
	#run 300000
	'''
	for line in commands.splitlines():
		lmp.command(line)

	lmp.file.close()

	os.system('%s -in %s.in -log %s.log' % (LAMMPS_DIR,run_name,run_name))
	os.chdir("../../")

	print("\n\n\n")
	for i in range(5):
		for j in range(50): print("#"),
		print("")
	print("\n\n\n")

def get_glimpse(s):
	print("Running glimpse for s = %s\n\n" % s)

	run_name = "glimpse_%s" % s
	parameters = "%s_output.tersoff" % s

	# Generate a large unit_cell for a test system
	test_system = mcsmrff_utils.get_test_system()
	test_system.name = run_name

	# Run MCSMRFF with initial parameters ONLY IF IT HASN'T BEEN ALREADY DONE!
	running_parameters = mcsmrff_files.read_params("parameters/%s" % parameters, exact=True)
	run_mcsmrff("%s" % run_name, test_system, running_parameters)