# This file contains a few functions for use when making the gradient
# optimization methods for MCSMRFF. It is now fully set up to run one
# LAMMPS simulation

import os
import numpy as np
import copy

import files
import structures

import mcsmrff_files
import mcsmrff_constants


def run_lammps(system, systems_by_composition, tersoff_atoms,
               lj_params, atom_list, tersoff_params, run_name):
    """
    Run a single step calculation using the given parameters.

    **Parameters**

        system:
        systems_by_composition:
        tersoff_atoms:
        lj_params:
        atom_list:
        tersoff_params:
        run_name:

    **Returns**

        Stuff.
    """

    # Ensure cwd/lammps/run_name exists
    directory_to_work = '%s/lammps/%s' % (os.getcwd(), run_name)
    directory_to_work = directory_to_work.split("/")
    for i, d in enumerate(directory_to_work):
        if d.strip() == "":
            continue
        if not os.path.isdir("/".join(directory_to_work[:i + 1])):
            os.mkdir("/".join(directory_to_work[:i + 1]))
    directory_to_work = "/".join(directory_to_work)
    os.chdir(directory_to_work)
    system.name = run_name

    # Make a file that has all the names of the systems we are using here
    files.write_lammps_data(system)

    # Begin code for the LAMMPS input file
    commands = ('''
    units real
    atom_style full
    pair_style hybrid/overlay lj/cut/coul/inout 0.2 3.5 15 tersoff
    bond_style harmonic
    angle_style harmonic
    dihedral_style opls
    special_bonds lj/coul 0.0 0.0 0.5

    boundary f f f
    read_data   ''' + run_name + '''.data
    ''').splitlines()

    # Indices of OPLS parameters
    tersoff_types = [t for t in system.atom_types
                     if t.index in tersoff_atoms]

    charges = lj_params[0]
    lj_sigma = lj_params[1]
    lj_epsilon = lj_params[2]

    # This is how we are reading in the 104 tersoff parameters
    # (i.e. m,gamma,N,D...etc)
    tersoff_strings, i, j = [], 0, 0
    # Concatenate the atom names and params
    while i < (len(tersoff_params) - 1):
        tmp = atom_list[j].split(",")
        for param in tersoff_params[i:i + 14]:
            tmp.append(str(param))
        tersoff_strings.append(tmp)
        i += 14
        j += 1

    # R and D are the distances from the center of the first atom, so this is
    # how we get the cutoff distance between those
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

    # Write to the lammps file the atom types and vdw radii
    for i in range(len(system.atom_types)):
        for j in range(i, len(system.atom_types)):
            type_i = system.atom_types[i]
            type_j = system.atom_types[j]
            commands.append('pair_coeff %d %d lj/cut/coul/inout %f %f %f' % (
                i + 1,
                j + 1,
                (type_i.vdw_e * type_j.vdw_e)**0.5,
                (type_i.vdw_r * type_j.vdw_r)**0.5,
                inner_cutoffs[(i, j)] if (i, j) in inner_cutoffs else 0.0)
            )
        commands.append('set type %d charge %f' % (i + 1, type_i.charge))

    # Generate a lmp object to make the LAMMPS input file
    lmp = structures.Struct()
    lmp.file = open(run_name + '.in', 'w')

    def writeline(line):
        lmp.file.write(line + '\n')

    lmp.command = writeline
    for line in commands:
        lmp.command(line)

    # Write the pair_coeff command
    mcsmrff_files.write_params(lj_params, atom_list,
                               tersoff_params, run_name, append="_input")
    mcsmrff_files.write_system_and_training_data(run_name, system,
                                                 systems_by_composition)
    lmp.command('pair_coeff * * tersoff ' + run_name + '_input.tersoff ' +
                (' '.join([(t.element_name if t in tersoff_types else 'NULL')
                          for t in system.atom_types])))

    for t in system.bond_types:
        lmp.command('bond_coeff %d  %f %f' % (t.lammps_type, t.e, t.r))
    for t in system.angle_types:
        lmp.command('angle_coeff %d %f %f' % (t.lammps_type, t.e, t.angle))
    for t in system.dihedral_types:
        lmp.command('dihedral_coeff %d  %f %f %f %f' %
                    ((t.lammps_type,) + t.e))

    commands = '''
    compute atom_pe all pe/atom
    compute sum_pe all reduce sum c_atom_pe
    neigh_modify once yes

    dump 2 all xyz 1 ''' + run_name + '''.xyz
    dump 1 all custom 1 ''' + run_name + '''.dump id type x y z fx fy fz c_atom_pe

    fix temp all nvt temp 10.0 10.0 100.0

    run 0
    undump 1
    '''
    for line in commands.splitlines():
        lmp.command(line)

    lmp.file.close()

    # Run the simulation the change directory back to the parent one
    os.system('%s -in %s.in -log %s.log >> out.log' %
              (mcsmrff_constants.lammps_mcsmrff, run_name, run_name))
    os.chdir("../../")


def calculate_error(run_name, natoms):
    """
    This function calculates the error between a MCSMRFF simulation of
    "run_name" and the stored training_forces and training_energies.

    **Parameters**

        run_name:
        natoms:

    **Returns**

        stuff.
    """
    energies, rms_forces = mcsmrff_files.parse_lammps_output(run_name, natoms)

    # Offset energies
    e_offset = energies[0]
    energies = [e - e_offset for e in energies]

    # Now we know the energy and rms_force, just need to get difference from
    # dft results and return
    f_training_forces = open("lammps/%s/%s_training_forces.txt" %
                             (run_name, run_name)).read().split("\n")
    f_training_energies = open("lammps/%s/%s_training_energies.txt" %
                               (run_name, run_name)).read().split("\n")
    training_energies = [float(x) for x in f_training_energies
                         if x.strip() != ""]
    training_rms_forces = np.array([float(x) for x in f_training_forces
                                    if x.strip() != ""])

    error_force = [((a - b) / a)**2
                   for a, b in zip(training_rms_forces, rms_forces)]
    # Issue with offset of first energy being 0. So divide by a+kT instead
    # kT = units.convert_energy("J","kcal/mol",10.0*constants.K_b)
    kT = 0.6
    error_energy = [((a - b) / (a + kT))**2
                    for a, b in zip(training_energies, energies)]

    N = len(rms_forces)
    error_force = np.sqrt(sum(error_force) / N)
    error_energy = np.sqrt(sum(error_energy) / N)

    return error_force, error_energy


def indices_of_desired_three_body(element_strings, three_body):
    """
    This allows the user to specify which three-body interactions they want to
    optimize parameters for.

    **Parameters**

        element_strings:
        three_body:

    **Returns**

        stuff:
    """
    element_list = np.array(element_strings).flatten()
    if three_body is None:
        index = range(len(element_list))
    else:
        concat_three_body = [", ".join(s) for s in three_body]
        index = []
        for i, s in enumerate(element_list):
            if s in concat_three_body:
                index.append(i)
        if index == []:
            raise Exception("No three_body set in three_body variable.")
    return index


def get_gradient(parameters, system, systems_by_composition, run_name,
                 perturbation=1.01, three_body=None,
                 tersoff=None, lj_coul=None, constant_charge=True,
                 tersoff_atoms=mcsmrff_constants.tersoff_atoms):
    """
    A function to calculate the gradient of the MCSMRFF parameters given an
    initial guess and atomic positions. This will perturb each parameter
    slightly and build up a gradient. There are three possible specifiers on
    what to perturb:
      three_body - Which three-body interactions you want to parameterize
      tersoff - Which tersoff parameters you want to parameterize
      lj_coul - Which of the charges, LJ-sigma and LJ-epsilon to perturb
    """

    if perturbation <= 1.0:
        raise Exception("Perturbation must be greater than 1.0.")

    p_lj, p_atoms, p_tersoff = [np.array(p) for p in parameters]

    indices_three_body = indices_of_desired_three_body(parameters[1],
                                                       three_body)

    N_tersoff = len(p_lj[0])
    if tersoff is None:
        tersoff = range(14)
    if lj_coul is None:
        lj_coul = range(N_tersoff * 3)
    if constant_charge:
        lj_coul = np.array([x for x in lj_coul if x not in range(N_tersoff)])

    atoms = copy.deepcopy(system)
    atoms.name = atoms.name + "_grad_calc"

    error_force, error_energy = calculate_error(run_name, len(atoms.atoms))
    error_0 = error_force

    nLJ, ntersoff = len(p_lj.flatten()), len(p_tersoff)
    gradient = np.zeros(ntersoff + nLJ)

    full_parameter_list = np.append(p_lj.flatten().copy(), p_tersoff.copy())
    for i, p in enumerate(full_parameter_list):
        # Only perturb systems we want to.  That is, leave gradient = 0 for
        # systems we do not want to perturb
        if i >= nLJ and int((i - nLJ) / 14) not in indices_three_body:
            continue
        if i >= nLJ and int((i - nLJ) % 14) not in tersoff:
            continue
        if i < nLJ and i not in lj_coul:
            continue

        # Hold a 1D array with all parameters we might perturb
        perturbed_parameters = full_parameter_list.copy()

        # In the case of m, we perturb it between 1 and 3. Everything else is
        # multiplied by perturbation
        if i >= nLJ and (i - nLJ) % 14 == 0:
            perturbed_parameters[i] = np.int64(3) if p == 1 else np.int64(1)
        else:
            perturbed_parameters[i] = p * perturbation

        p1 = perturbed_parameters[:nLJ].reshape((-1, len(p_lj)))
        p2 = perturbed_parameters[nLJ:]

        run_lammps(atoms, systems_by_composition, tersoff_atoms,
                   p1, parameters[1], p2, "%s_grad_calc" % run_name)

        # Use force
        error_1 = calculate_error("%s_grad_calc" % run_name,
                                  len(atoms.atoms))[0]
        gradient[i] = np.float64(error_1 - error_0)

    return gradient
