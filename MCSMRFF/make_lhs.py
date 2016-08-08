####################################################################################################################################################################################
## This file generates a latin hypercube sample for the tersoff parameters only
## Be sure that you have pulled from clancelot to obtain the doe_lhs script
## create_lhs(sample_max,sample_min,num_samples=1) - A function to generate tersoff parameters for the number of samples denoted by num_samples
#                                      sample_max : This is a list of 13 numbers for the maximum of each of the 13 tersoff parameters 
#							(gamma,lambda3,c,d,costheta0,n,beta,lambda2,B,R,D,lambda1,A) in that order (default is in place)
#				       sample_min : This is a list of 13 numbers for the minimum of each of the 13 tersoff parameters 
#							(gamma,lambda3,c,d,costheta0,n,beta,lambda2,B,R,D,lambda1,A) in that order (default is in place)
#                                      num_samples : the total number of samples you want to generate (default is 1)
#
####################################################################################################################################################################################

import random
import doe_lhs as lhs

def create_lhs(sample_max=[10,10,1000,1000,1,10,10,10,100000,5,2,15,15],sample_min=[0,0,0,0,-1,0,0,0,0,2,0,0,0],num_samples=1): 

	lhs_points = lhs.lhs(351,num_samples) ##since there are 27 total entires and 13 parameters per entry we must sample 351 numbers
	lhs_list = []
	OldRange = 1  
	OldMax = 1
	OldMin = 0
	for sample in lhs_points:
		sample_list = []
		i = 0
		for OldValue in sample:
			if (i == 0) : #gamma
				sample_list.append(random.choice([1,3])) #This lets us append a 1 or 3 to the front for the m variable without changing LHS applicability 
				NewValue = (((OldValue - OldMin) * (sample_max[i]-sample_min[i])) / OldRange) + sample_min[i]
				i += 1
			elif (i == 9): #R
				NewValue = (((OldValue - OldMin) * (sample_max[i]-sample_min[i])) / OldRange) + sample_min[i]
				i += 1
				R        = NewValue
			elif (i == 10): #D
				NewValue = (((OldValue - OldMin) * (sample_max[i]-sample_min[i])) / OldRange) + sample_min[i]
				i += 1
				while NewValue > R:
					NewValue = random.uniform(NewMin,NewMax)
			elif (i == 12): #A
				NewValue = (((OldValue - OldMin) * (sample_max[i]-sample_min[i])) / OldRange) + sample_min[i]
				i = 0
			else : 
				NewValue = (((OldValue - OldMin) * (sample_max[i]-sample_min[i])) / OldRange) + sample_min[i]
				i += 1
			sample_list.append(NewValue)
		lhs_list.append(sample_list)
	return lhs_list

x = create_lhs()
print x
