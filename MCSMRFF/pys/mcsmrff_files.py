##############################################################################
# This file holds any helper functions for file I/O
#
# parse_tersoff_file(filename) - A function that, given a file name, will
# parse the tersoff data into a list of strings
#
# read_params(run_name, exact=False) - A function to read in MCSMRFF
#    parameters from a file with the name input_run_name.tersoff
#    Using exact=True, it'll instead read in run_name.tersoff
#    This function returns three lists:
#        lj_paramlist, atom_list, tersoff_paramlist
#
# lj_paramlist - A list of lists.  The first is the charges, second is the
#                LJ_sigma value and third is the LJ_epsilon values.
# atom_list - A 1D list of strings holding all the atomic three-body
#             interactions
# tersoff_paramlist - A 1D list of floats holding the tersoff parameters for
#                     all three-body interactions
#
# write_params(lj_params,atom_list,tersoff_params,filename,append="_o") -
#     A function to write out parameters to a file
#
# write_param_line(tersoff,i) - Write a single line of a tersoff input file
#     (given the tersoff list and the index of the line. primarily written
#     in conjunction with write_params
#
# write_atom_line(atom,i) - Write a single line of atoms (given the atom list
#                           and the index of the line. primarily written in
#                           conjunction with write_params
#
# write_system_and_training_data(run_name, system, systems_by_composition) -
#     A function to write a system and training data to files
#
# parse_lammps_output(run_name, natoms) - A function to parse a LAMMPS output
#     and return the single point energy and RMS force
##############################################################################

import numpy as np

import files
import mcsmrff_utils


# A function to parse the tersoff parameters in the MCSMRFF file
def parse_tersoff_file(filename):
    tersoff_file = open(filename).read().split("\n")
    i = 0
    tersoff_strings, line = [], []
    while i < len(tersoff_file):
        if tersoff_file[i].strip() == "" or tersoff_file[i].strip()[0] == "#":
            i += 1
            continue
        tline = tersoff_file[i].strip().split()
        first = tline[0]
        # If the first 'word' is not a number,
        # then we are starting a new line of tersoff strings
        if not mcsmrff_utils.is_float(first):
            if line == []:
                line = tline
                i += 1
                continue
            else:
                tersoff_strings.append(line)
                line = tline
        else:
            line += tline
        i += 1
    tersoff_strings.append(line)
    return tersoff_strings


# A function to read in a MCSMRFF parameter file
def read_params(run_name, exact=False):
    """
    Read in the MCSMRFF style parameter file.

    **Parameters**

        run_name: *str*
            The name of the tersoff file.
        exact: *bool, optional*
            Whether the path supplied is an exact path or not.

    **Return**

        lj_paramlist: *list, list, float*
            Three lists, holding (1) charges, (2) LJ sigma, and (3) LJ epsilon
        atom_list: *list, str*
            1D array holding the different 3-body interactions within the file
        tersoff_paramlist: *list, float*
            1D array holding the different tersoff parameters
    """
    if ".tersoff" in run_name:
        run_name = run_name.split(".tersoff")[0]
    if not exact:
        filename = "input_%s.tersoff" % run_name
    else:
        filename = run_name + ".tersoff"
    print("Reading in parameter file %s..." % filename)

    lj_paramlist = []
    # This gives you the lj_sigma and epsilon
    # you dont need this but its here just in case
    for line in open(filename):
        if line.startswith('# Charges:'):
            charges = [float(x) for x in line.split()[2:]]
        if line.startswith('# LJ-sigma:'):
            lj_sigma = [float(x) for x in line.split()[2:]]
        if line.startswith('# LJ-epsilon:'):
            lj_epsilon = [float(x) for x in line.split()[2:]]
    lj_paramlist.append(charges)
    lj_paramlist.append(lj_sigma)
    lj_paramlist.append(lj_epsilon)

    tersoff_paramlist = []
    atom_list = []

    tersoff_strings = parse_tersoff_file(filename)

    for atom_triplet in tersoff_strings:
        for param in atom_triplet[3:]:
            tersoff_paramlist.append(np.float64(param))
        for atom in atom_triplet[:3]:
            # list of all atoms (Pb Pb Cl, Pb Cl Cl, etc)
            atom_list.append(atom)

    # Returns a list of all the parameters for both LJ and Tersoff
    # as well as the atoms used
    return lj_paramlist, atom_list, tersoff_paramlist


# Read in a list of LJ and tersoff in order as shown above in
# read_params and write this to a file
def write_params(lj_params, atom_list, tersoff_params, filename, append="_o"):
    """
    Write out parameters to a parameter file.

    **Parameters**

        lj_paramlist: *list, list, float*
            Three lists, holding (1) charges, (2) LJ sigma, and (3) LJ epsilon
        atom_list: *list, str*
            1D array holding the different 3-body interactions within the file
        tersoff_paramlist: *list, float*
            1D array holding the different tersoff parameters
        filename: *str*
            The name of the output file
        append: *str, optional*
            A string to append to the end of the filename

    **Return**

        None
    """
    ofile = open(filename + append + ".tersoff", "w")
    ofile.write('''# Generated by min_style params, run = 'None'. LAMMPS units = real.
# Error metric=0.0
# Charges:           ''' + '''    '''.join([str(x) for x in lj_params[0]]) + '''
# LJ-sigma:          ''' + '''    '''.join([str(x) for x in lj_params[1]]) + '''
# LJ-epsilon:        ''' + '''    '''.join([str(x) for x in lj_params[2]]) + '''
# i, j, k,       m,        gamma,    lambda3,  c,        d,        costheta0,
#                n,        beta,     lambda2,  B,        R,        D,''' +
                '''        lambda1,  A\n''')

    i, j = 0, 0
    # Write the output of the tersoff params, atoms in the list, etc.
    while i < (len(tersoff_params) - 1):
        ofile.write(_write_atom_line(atom_list, j) +
                    _write_param_line(tersoff_params, i))
        i += 14
        j += 3
    return


# Write a single line of a tersoff input file (given the tersoff list
# and the index of the line. primarily written in conjunction with write_params
def _write_param_line(tersoff, i):
    return (("%f   %f   %f   %f   %f   %f" +
            "\n\t\t %f   %f   %f   %f   %f   %f   %f   %f\n\n") %
            (tersoff[i],
             tersoff[i + 1],
             tersoff[i + 2],
             tersoff[i + 3],
             tersoff[i + 4],
             tersoff[i + 5],
             tersoff[i + 6],
             tersoff[i + 7],
             tersoff[i + 8],
             tersoff[i + 9],
             tersoff[i + 10],
             tersoff[i + 11],
             tersoff[i + 12],
             tersoff[i + 13]))


# Write a single line of atoms (given the atom list and the index of the line.
# primarily written in conjunction with write_params
def _write_atom_line(atom, i):
    return ("%s %s %s\t" % (atom[i], atom[i + 1], atom[i + 2]))


# NOTE - Systems are in order of "left" to "right" on the x axis
def write_system_and_training_data(run_name, system, systems_by_composition):
    system.name = run_name

    # Make a file that has all the names of the systems we are using here
    files.write_lammps_data(system)

    f = open(system.name + '_training_data.txt', 'w')
    for composition in systems_by_composition:
        for s in systems_by_composition[composition]:
            f.write(s.name + '\n')
    f.close()

    # write rms_forces to a file
    # DFT forces
    f = open(system.name + '_training_forces.txt', 'w')
    for m in system.molecules:
        f.write("%e\n" %
                np.sqrt(
                    np.sum(
                        [a.fx**2 + a.fy**2 + a.fz**2 for a in m.atoms]
                    ) / np.float64(len(m.atoms))))
    f.close()

    # Write energies to a file
    # DFT energies these are what we are comparing against and
    # calculating energy error from
    f = open(system.name + '_training_energies.txt', 'w')
    for m in system.molecules:
        f.write("%e\n" % (m.energy))
    f.close()


# A function to parse a LAMMPS output and return the single point
# energy and RMS force
def parse_lammps_output(run_name, natoms, buffer_len=100.0):
    # Read in forces
    force_file = open("lammps/%s/%s.dump" %
                      (run_name, run_name)).read().split("\n")
    i = 0
    while i < len(force_file) and "ITEM: ATOMS" not in force_file[i]:
        i += 1

    ids = {name: index for index, name in enumerate(force_file[i].split()[2:])}

    # Now i points to the first line of data from the dump file.
    # Start reading in
    i += 1
    energies, forces = {}, {}
    while i < len(force_file):
        line = force_file[i].strip()
        if line == "":
            i += 1
            continue
        else:
            values = line.split()
            x = np.float64(values[ids["x"]])
            fx = np.float64(values[ids["fx"]])
            fy = np.float64(values[ids["fy"]])
            fz = np.float64(values[ids["fz"]])
            pe = np.float64(values[ids["c_atom_pe"]])

            index = int((x + 100.0) / 1000)
            if index not in energies:
                energies[index] = [1, pe]
                forces[index] = [1, fx**2 + fy**2 + fz**2]
            else:
                energies[index][0] += 1
                energies[index][1] += pe
                forces[index][0] += 1
                forces[index][1] += fx**2 + fy**2 + fz**2
            i += 1

    energy = sorted([(j, v) for j, v in energies.items()], key=lambda x: x[0])
    energies = np.array([e[1][1] for e in energy])

    rms_force = sorted([(j, v) for j, v in forces.items()], key=lambda x: x[0])
    rms_forces = np.array([np.sqrt(e[1][1] / e[1][0]) for e in rms_force])

    return energies, rms_forces


# Example parameter file
#####################################################################################################
# Generated by min_style params, run = 'june15'. LAMMPS units = real.
# Error metric=1.00602
# Charges:           0.4       -0.2 
# LJ-sigma:            3          2 
# LJ-epsilon:        0.1        0.1 
# i, j, k,       m,        gamma,    lambda3,  c,        d,        costheta0,
#                n,        beta,     lambda2,  B,        R,        D,        lambda1,  A
#Pb Pb Pb         3         1   1.33856   4240.71    1.0451  0.248252
#                 1 1.4621e-07   2.77144   14721.8   3.00116 0.0721795    4.8965      2142

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