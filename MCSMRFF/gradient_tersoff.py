from merlin import *
import shutil, hashlib, re
import os, sys
#format: % python gradient_tersoff.py name_of_input_file
#get the run name or call this run test
run_name = sys.argv[1] if len(sys.argv)>1 else 'test'

#indices of OPLS parameters
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

if not os.path.isdir("lammps"): os.mkdir("lammps") #go to the lammps folder, or make one

shutil.copy('%s/input_%s.tersoff'% (os.getcwd(),run_name), '%s/lammps/%s_input.tersoff'%(os.getcwd(),system.name)) #NOTE: This is reading the input file from the parent directory named input_run_name.tersoff and copying it to the lammps dir

os.chdir('lammps')
files.write_lammps_data(system) #make a file that has all the names of the systems we are using here
f = open(system.name+'_data.txt', 'w')
for composition in systems_by_composition:
	for s in systems_by_composition[composition]:
		f.write(s.name+'\n')
f.close()



# write forces to a file
f = open(system.name+'_forces.txt', 'w') #DFT forces
for a in system.atoms:
	f.write("%e\n%e\n%e\n" % (a.fx, a.fy, a.fz) )
f.close()
# write energies to a file 
f = open(system.name+'_energies.txt', 'w') #DFT energies these are what we are comparing against and calculating energy error from
for m in system.molecules:
	f.write("%e\n" % (m.energy) )
f.close()

commands = ('''units real
atom_style full
pair_style hybrid/overlay lj/cut/coul/inout 0.2 3.5 15 tersoff
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5

boundary f f f
read_data	'''+system.name+'''.data
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

for line in open(system.name+'_input.tersoff'): #reading in the parameters from the input file
	if line.startswith('# Charges:'): charges = [float(x) for x in line.split()[2:]]
	if line.startswith('# LJ-sigma:'): lj_sigma = [float(x) for x in line.split()[2:]]
	if line.startswith('# LJ-epsilon:'): lj_epsilon = [float(x) for x in line.split()[2:]]

#this is how we are reading in the 104 tersoff parameters (i.e. m,gamma,N,D...etc)
tersoff_strings = re.findall('\n'+('(\S+) +'*9)[:-2]+' *\n +'+('(\S+) +'*8)[:-2], open(system.name+'_input.tersoff').read()) #re is regex library which is a way of pseudo-code to parse strings

#R and D are the distances from the center of the first atom, so this is how we get the cutoff distance between those
inner_cutoffs = {}
for type_i in tersoff_types:
	for type_j in tersoff_types:
		for s in tersoff_strings:
			print s
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
lmp.file = open(system.name+'.in', 'w')
def writeline(line):
	lmp.file.write(line+'\n')
lmp.command = writeline
for line in commands:
	lmp.command(line)

#write the pair_coeff command
lmp.command('pair_coeff * * tersoff '+system.name+'_input.tersoff '+(' '.join([ (t.element_name if t in tersoff_types else 'NULL') for t in system.atom_types])))

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

dump 1 all custom 1 '''+system.name+'''.dump id type x y z fx fy fz 
run 0
undump 1
'''
for line in commands.splitlines():
	lmp.command(line)

lmp.file.close()
os.system('/fs/home/afh72/lammps/lammps-7Dec15/src/lmp_serial -in %s.in -log %s.log' % (system.name,system.name))

