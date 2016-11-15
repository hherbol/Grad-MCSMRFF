import random
import pyDOE.doe_lhs as lhs


def create_lhs(
        sample_max=[10, 10, 1000, 1000, 1, 10, 10, 10, 100000, 5, 2, 15, 15],
        sample_min=[0, 0, 0, 0, -1, 0, 0, 0, 0, 2, 0, 0, 0],
        num_samples=1):
    """
    Generate a latin hypercube sample for the tersoff parameters.
    The order in which parameters are stored are as follows:
        gamma, lambda3, c, d, cos(theta_0), n,
        beta, lambda2, B, R, D, lambda1, A

    **Parameters**

        sample_max: *list, float, optional*
            The maximum for each of the 13 tersoff parameters.
        sample_min: *list, float, optional*
            The minimum for each of the 13 tersoff parameters.
        num_samples: *int, optional*
            The number of samples you wish to generate.

    **Returns**

        Stuff.
    """

    # Since there are 27 total entires and 13 parameters per entry we must
    # sample 351 numbers
    lhs_points = lhs.lhs(351, num_samples)
    lhs_list = []
    OldRange = 1
    OldMin = 0

    for sample in lhs_points:
        sample_list = []
        i = 0
        for OldValue in sample:
            # Gamma
            if (i == 0):
                # This lets us append a 1 or 3 to the front for the m variable
                # without changing LHS applicability
                sample_list.append(random.choice([1, 3]))
                NewValue = (
                    (
                        (OldValue - OldMin) * (sample_max[i] - sample_min[i])
                    ) / OldRange
                ) + sample_min[i]

                i += 1
            # R
            elif (i == 9):
                NewValue = (
                    (
                        (OldValue - OldMin) * (sample_max[i] - sample_min[i])
                    ) / OldRange
                ) + sample_min[i]

                i += 1
                R = NewValue
            # D
            elif (i == 10):
                NewValue = (
                    (
                        (OldValue - OldMin) * (sample_max[i] - sample_min[i])
                    ) / OldRange
                ) + sample_min[i]

                i += 1

                if NewValue > R:
                    raise Exception("Generated D>R, which is unrealistic.")
            # A
            elif (i == 12):
                NewValue = (
                    (
                        (OldValue - OldMin) * (sample_max[i] - sample_min[i])
                    ) / OldRange
                ) + sample_min[i]

                i = 0
            else:
                NewValue = (
                    (
                        (OldValue - OldMin) * (sample_max[i] - sample_min[i])
                    ) / OldRange
                ) + sample_min[i]

                i += 1
            sample_list.append(NewValue)

        lhs_list.append(sample_list)

    return lhs_list


if __name__ == "__main__":
    x = create_lhs()
    print x

##############################################################################
# This file generates a latin hypercube sample for the tersoff parameters only
# Be sure that you have pulled from clancelot to obtain the doe_lhs script
# create_lhs(sample_max,sample_min,num_samples=1) - A function to generate
# tersoff parameters for the number of samples denoted by num_samples
#    sample_max : This is a list of 13
#    numbers for the maximum of each of the 13 tersoff parameters
#    (gamma,lambda3,c,d,costheta0,n,beta,lambda2,B,R,D,lambda1,A) in that order
#        (default is in place)
#    sample_min : This is a list of 13 numbers for the minimum of each of the
#        13 tersoff parameters
#    (gamma,lambda3,c,d,costheta0,n,beta,lambda2,B,R,D,lambda1,A) in that order
#        (default is in place)
#    num_samples : the total number of samples you want to generate
#        (default is 1)
##############################################################################
