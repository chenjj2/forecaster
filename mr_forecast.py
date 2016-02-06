import numpy as np
from scipy.stats import norm
from scipy.stats import truncnorm 

## constant
mearth2mjup = 317.828
mearth2msun = 333060.4
rearth2rjup = 11.21
rearth2rsun = 109.2

mlower = 3e-4
mupper = 3e5

## hyper file
hyper_file = 'h4_thin_hyper.dat'

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
def Mpost2R(mass, unit='Earth'):
	# unit
	if unit == 'Earth':
		pass
	elif unit == 'Jupiter':
		mass = mass * mearth2mjup
	else:
		print "Warning: input unit must be 'Earth' or 'Jupiter'. Using 'Earth' as default."

	# mass range
	if np.min(mass) < 3e-4 or np.max(mass) > 3e5:
		print 'Error: mass range out of model expectation. Returning None.'
		return None

	## convert to radius
	sample_size = len(mass)
	logm = np.log10(mass)
	prob = np.random.random(sample_size)
	logr = np.ones_like(logm)
	
	all_hyper = np.loadtxt(hyper_file)
	hyper_ind = np.random.randint(low = 0, high = np.shape(all_hyper)[0], size = sample_size)	
	hyper = all_hyper[hyper_ind,:]

	for i in range(sample_size):
		logr[i] = piece_linear(hyper[i], logm[i], prob[i])

	radius_sample = 10.** logr

	## convert to right unit
	if unit == 'Jupiter':
		radius = radius_sample / rearth2rjup
	else:
		radius = radius_sample 

	return radius


### given mass statistics, yield radius stats, assuming Normal Distribution
def Mstat2R(mean, std, unit='Earth', sample_size=100):	
	# unit
	if unit == 'Earth':
		pass
	elif unit == 'Jupiter':
		mean = mean * mearth2mjup
		std = std * mearth2mjup
	else:
		print "Warning: input unit must be 'Earth' or 'Jupiter'. Using 'Earth' as default."

	# draw samples
	print 'Assuming normal distribution truncated at the mass range limit of the model.'
	mass = truncnorm.rvs( (mlower-mean)/std, (mupper-mean)/std, loc=mean, scale=std, size=sample_size)		
	radius = Mpost2R(mass, unit='Earth')

	if unit == 'Jupiter':
		radius = radius / rearth2rjup

	return np.mean(radius), np.std(radius)


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
def Rpost2M(radius, unit='Earth', grid_size = 1e3):
	# unit
	if unit == 'Earth':
		pass
	elif unit == 'Jupiter':
		radius = radius * rearth2rjup
	else:
		print "Warning: input unit must be 'Earth' or 'Jupiter'. Using 'Earth' as default."

	# mass range
	if np.min(radius) <= 0.:
		print 'Error: there cannot be negative radius. Returning None.'
		return None

	# sample_grid
	if grid_size < 10:
		print 'Warning: the sample grid is too sparse. Using 10 sample grid instead.'
		grid_size = 10

	## convert to mass
	sample_size = len(radius)
	logr = np.log10(radius)
	logm = np.ones_like(logr)

	all_hyper = np.loadtxt(hyper_file)
	hyper_ind = np.random.randint(low = 0, high = np.shape(all_hyper)[0], size = sample_size)	
	hyper = all_hyper[hyper_ind,:]

	logm_grid = np.linspace(-3.522, 5.477, grid_size)

	for i in range(sample_size):
		prob = ProbRGivenM(logr[i], logm_grid, hyper[i,:])
		logm[i] = np.random.choice(logm_grid, size=1, p = prob)

	mass_sample = 10.** logm

	## convert to right unit
	if unit == 'Jupiter':
		mass = mass_sample / mearth2mjup
	else:
		mass = mass_sample
	
	return mass

### given R statistics, yield mass stat
def Rstat2M(mean, std, unit='Earth', sample_size=100):	
	# unit
	if unit == 'Earth':
		pass
	elif unit == 'Jupiter':
		mean = mean * rearth2rjup
		std = std * rearth2rjup
	else:
		print "Warning: input unit must be 'Earth' or 'Jupiter'. Using 'Earth' as default."

	# draw samples
	print 'Assuming normal distribution truncated from zero on.'
	radius = truncnorm.rvs( (0.-mean)/std, np.inf, loc=mean, scale=std, size=sample_size)		
	mass = Rpost2M(radius, unit='Earth')

	if mass is None:
		return None

	if unit=='Jupiter':
		mass = mass / mearth2mjup

		return np.mean(mass), np.std(mass)

