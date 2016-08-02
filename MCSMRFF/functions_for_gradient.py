from merlin import *
import shutil, hashlib, re
import os, sys
#This file contains a few functions for use when making the gradient optimization methods for MCSMRFF

def read_params(filename): #read in a file and output a list of all the TERSOFF parameters for a file of the format shown below in the order shown below
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
	#this is how we are reading in the 112 tersoff parameters INCLUDING m...DONT MINIMIZE m (i.e. m,gamma,N,D...etc) only minimize the actual 104 parameters
	tersoff_strings = re.findall('\n'+('(\S+) +'*9)[:-2]+' *\n +'+('(\S+) +'*8)[:-2], open(filename).read()) #re is regex library which is a way of pseudo-code to parse strings
	for atom_triplet in tersoff_strings:
		for param in atom_triplet[3:]:
			tersoff_paramlist.append(float(param))
	
	return lj_paramlist,tersoff_paramlist #returns a list of all the parameters for both LJ and Tersoff

def write_params(lj_params,tersoff_params,filename): #read in a list of LJ and tersoff in order as shown above in read_params and write this to a file
	ofile = open(filename,"w")
	ofile.write('''# Generated by min_style params, run = 'None'. LAMMPS units = real.
# Error metric=0.0
# Charges:           '''+str(lj_params[0][0])+'''       '''+str(lj_params[0][1])+''' 
# LJ-sigma:          '''+str(lj_params[1][0])+'''       '''+str(lj_params[1][1])+'''
# LJ-epsilon:        '''+str(lj_params[2][0])+'''       '''+str(lj_params[2][1])+''' 
# i, j, k,       m,        gamma,    lambda3,  c,        d,        costheta0,
#                n,        beta,     lambda2,  B,        R,        D,        lambda1,  A\n''')
	ofile.write('''
Pb Pb Pb         '''+write_line(tersoff_params,0)+'''
Pb Pb Cl         '''+write_line(tersoff_params,14)+'''
Pb Cl Pb         '''+write_line(tersoff_params,28)+'''
Pb Cl Cl         '''+write_line(tersoff_params,42)+'''
Cl Pb Pb         '''+write_line(tersoff_params,56)+'''
Cl Pb Cl         '''+write_line(tersoff_params,70)+'''
Cl Cl Pb         '''+write_line(tersoff_params,84)+'''
Cl Cl Cl         '''+write_line(tersoff_params,98))
	return

def write_line(tersoff,i): #write a single line of a tersoff input file (given the tersoff list and the index of the line. primarily written in conjunction with write_params
	return ("%f   %f   %f   %f   %f   %f\n\t\t %f   %f   %f   %f   %f   %f   %f   %f\n\n" % (tersoff[i],tersoff[i+1],tersoff[i+2],tersoff[i+3],tersoff[i+4],tersoff[i+5],tersoff[i+6],tersoff[i+7],tersoff[i+8],tersoff[i+9],tersoff[i+10],tersoff[i+11],tersoff[i+12],tersoff[i+13]))

lj,tersoff = read_params("input_test.tersoff")
write_params(lj,tersoff,"test.tersoff")






