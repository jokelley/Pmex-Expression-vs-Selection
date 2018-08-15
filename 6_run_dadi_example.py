################################################################################################
#
#           6. Estimating demographic models with dadi                       
#
#
################################################################################################

# This is an example script of how to run dadi
# The demographic model is defined in demographic_model_splitaysmmig.py
# This script is for actually using empirical data and the defined demographic model in dadi

import numpy
from numpy import array

import dadi

# In demographic_model_splitaysmmig.py, we've defined a custom model for this problem
import demographic_model_splitaysmmig

# Load the data
data = dadi.Spectrum.from_file('puy_shared_sfs.txt', mask_corners=True)

# Mask singletons
data.mask[1,0] = True
data.mask[0,1] = True
data.mask[19,0] = True
data.mask[0,19] = True
ns = data.sample_sizes

# These are the grid point settings will use for extrapolation.
pts_l = [20,30,40]

# Define our custom demographic model
func = demographic_model_splitaysmmig.split_asym_mig
# Define the starting array for parameter optimization (this should be randomized for subsequent optimizations)
params = array([1,1,1,1,1])
# The upper_bound array is for use in optimization
upper_bound = [100, 100, 50,50,50]
lower_bound = [1e-3, 1e-3, 0, 0, 0]

# Makde the extrapolating version of our demographic model function
func_ex = dadi.Numerics.make_extrap_log_func(func)
# Calculate the model AFS
model = func_ex(params, ns, pts_l)
# Likelihood of the data given the model AFS
ll_model = dadi.Inference.ll_multinom(model, data)
print 'Model log-likelihood:', ll_model
# The optimal value of theta given the model
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print 'Model theta:', theta

# Perturb our parameter array before optimization. This does so by taking each
# parameter a up to a factor of two up or down
p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound)
# Do the optimization
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, 
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(params), maxiter=100)

# Print Optimized parameters, log-likelihood, and theta to separate files
print 'Optimized parameters', repr(popt)
model = func_ex(popt, ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(model, data)
print 'Optimized log-likelihood:', ll_opt
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print 'Optimal theta:', theta
f = open('ll_opt1.txt', 'w')
f.write('Optimized log-likelihood1 = ' + repr(ll_opt) + '\n')
f.close()
f = open('theta1.txt', 'w')
f.write('Optimized theta1 = ' + repr(theta) + '\n')
f.close()
f = open('opt_param1.txt', 'w')
f.write('Optimized parameters1 = ' + repr(popt) + '\n')
f.close()