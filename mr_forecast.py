import numpy as np
from scipy.stats import norm
from scipy.stats import truncnorm 

## constant
mearth2mjup = 317.828
mearth2msun = 333060.4
rearth2rjup = 11.21
rearth2rsun = 109.2

## hyper file
hyper_file = 'h4_thin_hyper.out'

### fix the number of different populations
n_pop = 4

### indicate which M belongs to population i given transition parameter
def indicate(M, trans, i):
	ts = np.insert(np.insert(trans, n_pop-1, np.inf), 0, -np.inf)
	ind = (M>=ts[i]) & (M<ts[i+1])
	return ind

### split hyper and derive c
def split_hyper_linear(hyper):
	c0, slope,sigma, trans = \
	hyper[0], hyper[1:1+n_pop], hyper[1+n_pop:1+2*n_pop], hyper[1+2*n_pop:]

	c = np.zeros_like(slope)
	c[0] = c0
	for i in range(1,n_pop):
		# trans[0] * slope[0] + c[0] = trans[0] * slope[1] + c[1]
		# c[1] = c[0] + trans[0] * (slope[0]-slope[1])
		c[i] = c[i-1] + trans[i-1]*(slope[i-1]-slope[i])

	return c, slope, sigma, trans

### model: straight line
def piece_linear(hyper, M, prob_R):
	c, slope, sigma, trans = split_hyper_linear(hyper)
	R = np.zeros_like(M)
	for i in range(4):
		ind = indicate(M, trans, i)
		mu = c[i] + M[ind]*slope[i]
		R[ind] = norm.ppf(prob_R[ind], mu, sigma[i])

	return R

### given mass distribution, yield radius distribution
def Mpost2R(mass, unit='earth', size = -1):

	## check input
	# mass type
	if type(mass) != np.ndarray or np.ndim(mass) != 1:
		print 'Error: input mass must be 1D numpy array. '
		return None

	# unit
	if unit == 'earth':
		pass
	elif unit == 'jupiter':
		mass = mass * mearth2mjup
	elif unit == 'sun':
		mass = mass * mearth2msun
	else:
		print 'Error: input unit must be earth/jupiter/sun. '
		return None

	# mass range
	if np.min(mass) < 1e-4 or np.max(mass) > 1e6:
		print 'Error: mass range out of model expectation. '
		return None

	# size
	if size == -1:
		sample_size = len(mass)
		mass_sample = mass
	elif size > 1:
		sample_size = int(size)
		mass_sample = np.random.choice(mass, size=sample_size, replace=True)
	else:
		print 'Error: size should be a positive integer. Default is the size of the input array.'
		return None


	## convert to radius
	logm = np.log10(mass_sample)
	prob = np.random.random(sample_size)
	logr = np.ones_like(logm)
	
	all_hyper = np.loadtxt(hyper_file)
	hyper_ind = np.random.randint(low = 0, high = np.shape(all_hyper)[0], size = sample_size)	
	hyper = all_hyper[hyper_ind,:]

	for i in range(sample_size):
		logr[i] = piece_linear(hyper[i], logm[i], prob[i])

	radius_sample = 10.** logr

	## convert to right unit
	if unit == 'earth':
		radius = radius_sample
	elif unit == 'jupiter':
		radius = radius_sample / rearth2rjup
	elif unit == 'sun':
		radius = radius_sample / rearth2rsun
	else:
		print 'Error: input unit must be earth/jupiter/sun. '
		return None
	
	return radius

### given mass statistics, yield radius distribution
# make it a gaussian if symmetric,
# or two different gaussian joint at the median if asymmetric
def Mstat2R(mean, std, down_std=-1, unit='earth', sample_size=100):
	print 'Assuming normal distribution truncated at the mass range limit of the model.'
	mass = truncnorm.rvs(1e-4, 1e6, loc=mean, scale=std, size=sample_size)		
	radius = Mpost2R(mass, unit='earth', size = -1)
	return radius


### p(radii|M)
def ProbRGivenM(radii, M, hyper):

	c, slope, sigma, trans = split_hyper_linear(hyper)
	prob = np.zeros_like(M)
	
	for i in range(4):
		ind = indicate(M, trans, i)
		mu = c[i] + M[ind]*slope[i]
		sig = sigma[i]
		prob[ind] = norm.pdf(radii, mu, sig)

	prob = prob/np.sum(prob)

	return prob

### given radius posterior, yield mass
def Rpost2M(radius, unit='earth', size=-1, grid_size = 1e3):
	
	## check input
	# radius type
	if type(radius) != np.ndarray or np.ndim(radius) != 1:
		print 'Error: input radius must be 1D numpy array. '
		return None

	# unit
	if unit == 'earth':
		pass
	elif unit == 'jupiter':
		radius = radius * rearth2rjup
	elif unit == 'sun':
		radius = radius * rearth2rsun
	else:
		print 'Error: input unit must be earth/jupiter/sun. '
		return None

	# mass range
	if np.min(radius) <= 0.:
		print 'Error: radius range out of model expectation. '
		return None

	# size
	if size == -1:
		sample_size = len(radius)
		radius_sample = radius
	elif size > 1:
		sample_size = int(size)
		radius_sample = np.random.choice(radius, size=sample_size, replace=True)
	else:
		print 'Error: size should be a positive integer. Default is the size of the input array.'
		return None

	# sample_grid
	if grid_size < 5:
		print 'Error: the sample grid is too sparse. Suggest using at least 100 sample grid.'
		return None

	## convert to mass
	logr = np.log10(radius_sample)
	logm = np.ones_like(logr)

	all_hyper = np.loadtxt(hyper_file)
	hyper_ind = np.random.randint(low = 0, high = np.shape(all_hyper)[0], size = sample_size)	
	hyper = all_hyper[hyper_ind,:]

	logm_grid = np.linspace(-4., 5.5, grid_size)

	for i in range(sample_size):
		prob = ProbRGivenM(logr[i], logm_grid, hyper[i,:])
		logm[i] = np.random.choice(logm_grid, size=1, p = prob)

	mass_sample = 10.** logm

	## convert to right unit
	if unit == 'earth':
		mass = mass_sample
	elif unit == 'jupiter':
		mass = mass_sample / mearth2mjup
	elif unit == 'sun':
		mass = mass_sample / mearth2msun
	else:
		print 'Error: input unit must be earth/jupiter/sun. '
		return None
	
	return mass

### given R statistics, yield mass distribution
def Rstat2M(mean, std, unit='earth', sample_size=100):
	print 'Assuming normal distribution truncated from zero to infinity'
	radius = truncnorm.rvs(0., np.inf, loc=mean, scale=std, size=sample_size)		
	mass = Rpost2M(radius, unit='earth', size = -1)
	return mass

