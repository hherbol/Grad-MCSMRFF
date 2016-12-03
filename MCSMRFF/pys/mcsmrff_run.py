# A function to run a MCSMRFF simulation using the given
#      (1) run_name, (2) system, and (3) parameter file
# Currently it is highly tuned to Pb, Cl, I systems
# Note, parameters must be a list of [LJ, atom_strs, tersoff]

import files
import sysconst
import subprocess
import jobs
# import lammps_job
import structures
import os
import sys
import shutil
# from sysconst import lammps_mcsmrff as LAMMPS_DIR
from hashlib import md5
from re import findall

import mcsmrff_files
import mcsmrff_utils

import units


# A function to run an LAMMPS Simulation. Requires a run name and a
# string of lammps code (run_name and input_script)
def job(run_name, input_script, system, parameters, queue=None, procs=1,
        email=None, pair_coeffs_included=True, hybrid_pair=False,
        hybrid_angle=False, TIP4P=False):
    """
    Wrapper to submitting a LAMMPs simulation.

    **Parameters**

        run_name: *str*
            Name of the simulation to be run.
        input_script: *str*
            Input script for LAMMPs simulation.
        system: :class:`structures.System`
            System object for our simulation.
        parameters:
            PARAMETERS FOR THIS MCSMRFF SIMULATION.
        tersoff_atoms:
            TERSOFF ATOMS IN THIS SIMULATION.
        queue: *str, optional*
            What queue to run the simulation on (queueing system dependent).
        procs: *int, optional*
            How many processors to run the simulation on.
        email: *str, optional*
            An email address for sending job information to.
        pair_coeffs_included: *bool, optional*
            Whether we have included the pair coefficients to be written
            to our lammps data file.
        hybrid_pair: *bool, optional*
            Whether to detect different treatments of pairing interactions
            amongs different atom types(True), or not (False).
        hybrid_angle: *bool, optional*
            Whether to detect different treatments of angles amongst different
            atom types (True), or not (False).
        TIP4P: *bool, optional*
            Whether to identify TIP4P settings within the lammps data file and
            update the input file (True), or not (False).

    **Returns**

        job: :class:`jobs.Job`
            If running locally, return the process handle, else return the
            job container.
    """
    if len(run_name) > 31 and queue is not None:
        raise Exception("Job name too long (%d) for NBS. Max character \
length is 31." % len(run_name))

    # Change to correct directory
    os.system('mkdir -p lammps/%s' % run_name)
    os.chdir('lammps/%s' % run_name)

    # Generate the lammps data file
    files.write_lammps_data(system,
                            pair_coeffs_included=False,
                            hybrid_pair=hybrid_pair,
                            hybrid_angle=hybrid_angle)
    mcsmrff_files.write_params(parameters[0], parameters[1], parameters[2],
                               run_name, append="")

    # Write the lammps input script. Expects lines of lammps code
    f = open(run_name + '.in', 'w')
    f.write(input_script)
    f.close()

    # Run the simulation
    if queue is None:
        cmd_to_run = sysconst.lmp_env_vars
        cmd_to_run = cmd_to_run + "\n%s -in %s.in -log %s.log"\
            % (sysconst.lammps_mcsmrff, run_name, run_name)

        process_handle = subprocess.Popen(cmd_to_run, shell=True)
        job_handle = jobs.Job(run_name, process_handle)
    else:
        cmd_to_run = sysconst.lammps_mcsmrff +\
            " -in " +\
            (os.getcwd() + '/' + run_name) +\
            ".in -echo log -log " +\
            (os.getcwd() + '/' + run_name) +\
            ".log"

        job_handle = jobs.submit_job(run_name,
                                     cmd_to_run,
                                     procs=procs,
                                     queue=queue,
                                     additional_env_vars=sysconst.lmp_env_vars,
                                     email=email,
                                     preface="mpi")

    # Copy run script
    fname = sys.argv[0]
    if '/' in fname:
        fname = fname.split('/')[-1]
    try:
        shutil.copyfile('../../%s' % fname, fname)
    except IOError:
        # Submitted a job oddly enough that sys.argv[0]
        # is not the original python file name, so don't do this
        pass

    # Return to the appropriate directory
    os.chdir('../..')

    return job_handle


def run_mcsmrff(run_name, system, parameters, tersoff_atoms,
                queue=None, procs=1, email=None,
                pair_coeffs_included=True, hybrid_pair=False,
                hybrid_angle=False, TIP4P=False, seed=None):
    if seed is None:
        seed = str(int(md5(run_name).hexdigest(), 16) % (2**16))
    else:
        seed = str(seed)
    system.name = run_name

    # Begin generating string to hold LAMMPS input
    commands = ('''units real
atom_style full
pair_style hybrid/overlay lj/cut/coul/inout 0.2 3.5 15 tersoff
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5

boundary p p p
read_data   ''' + run_name + '''.data

''').splitlines()

    # Removed HN for no H3 input types
    tersoff_types = [t for t in system.atom_types
                     if t.index in tersoff_atoms]

    elems_by_index = [(t.lammps_type, t.element) for t in system.atom_types]
    elems_by_index = sorted(elems_by_index, key=lambda x: [0])
    elems_by_index = ' '.join([units.elem_i2s(e[1]) for e in elems_by_index])

    # Grab the parameters. Could in theory just use "parameters" variable,
    # but too lazy to change code
    for line in open('%s.tersoff' % run_name):
        if line.startswith('# Charges:'):
            charges = [float(x) for x in line.split()[2:]]
        if line.startswith('# LJ-sigma:'):
            lj_sigma = [float(x) for x in line.split()[2:]]
        if line.startswith('# LJ-epsilon:'):
            lj_epsilon = [float(x) for x in line.split()[2:]]

    tersoff_strings = findall('\n' +
                              ('(\S+) +' * 9)[:-2] +
                              ' *\n +' +
                              ('(\S+) +' * 8)[:-2],
                              open(run_name + '.tersoff').read())
    inner_cutoffs = {}
    for type_i in tersoff_types:
        for type_j in tersoff_types:
            for s in tersoff_strings:
                types = s[:3]
                R, D = float(s[13]), float(s[14])
                if types == (type_i.element_name,
                             type_j.element_name,
                             type_j.element_name):
                    inner_cutoffs[(type_i, type_j)] = R + D

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
            commands.append('pair_coeff %d %d lj/cut/coul/inout %f %f %f' %
                            (i + 1, j + 1,
                             (type_i.vdw_e * type_j.vdw_e)**0.5,
                             (type_i.vdw_r * type_j.vdw_r)**0.5,
                             inner_cutoffs[(i, j)]
                             if (i, j) in inner_cutoffs else 0.0))
        commands.append('set type %d charge %f' % (i + 1, type_i.charge))

    # Generate a LAMMPS structure for an input file

    lmp = structures.Struct()
    lmp.file = open(run_name + '.in', 'w')

    def writeline(line):
        lmp.file.write(line + '\n')
    lmp.command = writeline

    for line in commands:
        lmp.command(line)

    commands.append('pair_coeff * * tersoff ' + run_name + '.tersoff ' +
                    (' '.join([(t.element_name
                              if t in tersoff_types else 'NULL')
                     for t in system.atom_types])))

    for t in system.bond_types:
        commands.append('bond_coeff %d  %f %f' % (t.lammps_type, t.e, t.r))
    for t in system.angle_types:
        commands.append('angle_coeff %d %f %f' % (t.lammps_type, t.e, t.angle))
    for t in system.dihedral_types:
        commands.append('dihedral_coeff %d  %f %f %f %f' %
                        ((t.lammps_type,) + t.e))

    commands.append('''
neigh_modify every 1 check yes delay 0
dump    1 all xyz 100 ''' + run_name + '''.xyz
dump_modify 1 element ''' + elems_by_index + '''
dump    2 all custom 100 ''' + run_name + '''2.dump type xu yu zu
thermo_style custom step temp press ke pe epair emol vol
thermo 1000

group mobile id > 0

velocity all create 10.0 ''' + seed + ''' rot yes dist gaussian
timestep 0.01

fix press all npt temp 10.0 10.0 10.0 iso 0.0 0.0 100.0
run 50000
unfix press

fix motion all nvt temp 10.0 300.0 100.0
run 200000
''')

    # job(run_name, commands, system, queue=None, hybrid_angle=False)
    J = job(run_name, "\n".join(commands), system, parameters, queue=queue,
            procs=procs, email=email,
            pair_coeffs_included=pair_coeffs_included,
            hybrid_pair=hybrid_pair, hybrid_angle=hybrid_angle, TIP4P=TIP4P)
    return J

# os.system('%s -in %s.in -log %s.log' % (LAMMPS_DIR, run_name, run_name))


def get_glimpse(s, tersoff_atoms):
    print("Running glimpse for s = %s\n\n" % s)

    run_name = "glimpse_%s" % s
    parameters = "%s_output.tersoff" % s

    # Generate a large unit_cell for a test system
    test_system = mcsmrff_utils.get_test_system()
    test_system.name = run_name

    # Run MCSMRFF with initial parameters ONLY IF IT HASN'T BEEN ALREADY DONE!
    running_parameters = mcsmrff_files.read_params("parameters/%s" %
                                                   parameters, exact=True)
    run_mcsmrff("%s" % run_name,
                test_system,
                running_parameters,
                tersoff_atoms)
