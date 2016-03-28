import numpy as np
from scipy.stats import norm
from scipy.stats import truncnorm 
import h5py 

## constant
mearth2mjup = 317.828
mearth2msun = 333060.4
rearth2rjup = 11.21
rearth2rsun = 109.2

## boundary
mlower = 3e-4
mupper = 3e5

## number of category
n_pop = 4

## read parameter file
hyper_file = 'fitting_parameters.h5'
h5 = h5py.File(hyper_file, 'r')
all_hyper = h5['hyper_posterior'][:]
h5.close()

## function
from func import piece_linear, ProbRGivenM, classification

##############################################

def Mpost2R(mass, unit='Earth', classify='No'):
	"""
	Forecast the Radius distribution given the mass distribution.

	Parameters
	---------------
	mass: one dimensional array
		The mass distribution.
	unit: string (optional)
		Unit of the mass. 
		Options are 'Earth' and 'Jupiter'. Default is 'Earth'.
	classify: string (optional)
		If you want the object to be classifed. 
		Options are 'Yes' and 'No'. Default is 'No'.
		Result will be printed, not returned.

	Returns
	---------------
	radius: one dimensional array
		Predicted radius distribution in the input unit.
	"""

	# mass input
	mass = np.array(mass)
	assert len(mass.shape) == 1, "Input mass must be 1-D."

	# unit input
	if unit == 'Earth':
		pass
	elif unit == 'Jupiter':
		mass = mass * mearth2mjup
	else:
		print "Input unit must be 'Earth' or 'Jupiter'. Using 'Earth' as default."

	# mass range
	if np.min(mass) < 3e-4 or np.max(mass) > 3e5:
		print 'Mass range out of model expectation. Returning None.'
		return None

	## convert to radius
	sample_size = len(mass)
	logm = np.log10(mass)
	prob = np.random.random(sample_size)
	logr = np.ones_like(logm)

	hyper_ind = np.random.randint(low = 0, high = np.shape(all_hyper)[0], size = sample_size)	
	hyper = all_hyper[hyper_ind,:]

	if classify == 'Yes':
		classification(logm, hyper[:,-3:])
		

	for i in range(sample_size):
		logr[i] = piece_linear(hyper[i], logm[i], prob[i])

	radius_sample = 10.** logr

	## convert to right unit
	if unit == 'Jupiter':
		radius = radius_sample / rearth2rjup
	else:
		radius = radius_sample 

	return radius



def Mstat2R(mean, std, unit='Earth', sample_size=1000, classify = 'No'):	
	"""
	Forecast the mean and standard deviation of radius given the mena and standard deviation of the mass.
	Assuming normal distribution with the mean and standard deviation truncated at the mass range limit of the model.

	Parameters
	---------------
	mean: float
		Mean (average) of mass.
	std: float
		Standard deviation of mass.
	unit: string (optional)
		Unit of the mass. Options are 'Earth' and 'Jupiter'.
	sample_size: int (optional)
		Number of mass samples to draw with the mean and std provided.
	Returns
	---------------
	mean: float
		Predicted mean of radius in the input unit.
	std: float
		Predicted standard deviation of radius.
	"""

	# unit
	if unit == 'Earth':
		pass
	elif unit == 'Jupiter':
		mean = mean * mearth2mjup
		std = std * mearth2mjup
	else:
		print "Input unit must be 'Earth' or 'Jupiter'. Using 'Earth' as default."

	# draw samples
	mass = truncnorm.rvs( (mlower-mean)/std, (mupper-mean)/std, loc=mean, scale=std, size=sample_size)	
	if classify == 'Yes':	
		radius = Mpost2R(mass, unit='Earth', classify='Yes')
	else:
		radius = Mpost2R(mass, unit='Earth')

	if unit == 'Jupiter':
		radius = radius / rearth2rjup

	return np.mean(radius), np.std(radius)



def Rpost2M(radius, unit='Earth', grid_size = 1e3, classify = 'No'):
	"""
	Forecast the mass distribution given the radius distribution.

	Parameters
	---------------
	radius: one dimensional array
		The radius distribution.
	unit: string (optional)
		Unit of the mass. Options are 'Earth' and 'Jupiter'.
	grid_size: int (optional)
		Number of grid in the mass axis when sampling mass from radius.
		The more the better results, but slower process.
	classify: string (optional)
		If you want the object to be classifed. 
		Options are 'Yes' and 'No'. Default is 'No'.
		Result will be printed, not returned.

	Returns
	---------------
	mass: one dimensional array
		Predicted mass distribution in the input unit.
	"""
	
	# unit
	if unit == 'Earth':
		pass
	elif unit == 'Jupiter':
		radius = radius * rearth2rjup
	else:
		print "Input unit must be 'Earth' or 'Jupiter'. Using 'Earth' as default."

	# mass range
	if np.min(radius) <= 0.:
		print 'There cannot be negative radius. Returning None.'
		return None

	# sample_grid
	if grid_size < 10:
		print 'The sample grid is too sparse. Using 10 sample grid instead.'
		grid_size = 10

	## convert to mass
	sample_size = len(radius)
	logr = np.log10(radius)
	logm = np.ones_like(logr)

	hyper_ind = np.random.randint(low = 0, high = np.shape(all_hyper)[0], size = sample_size)	
	hyper = all_hyper[hyper_ind,:]

	logm_grid = np.linspace(-3.522, 5.477, grid_size)

	for i in range(sample_size):
		prob = ProbRGivenM(logr[i], logm_grid, hyper[i,:])
		logm[i] = np.random.choice(logm_grid, size=1, p = prob)

	mass_sample = 10.** logm

	if classify == 'Yes':
		classification(logm, hyper[:,-3:])

	## convert to right unit
	if unit == 'Jupiter':
		mass = mass_sample / mearth2mjup
	else:
		mass = mass_sample
	
	return mass



def Rstat2M(mean, std, unit='Earth', sample_size=1e3, grid_size=1e3, classify = 'No'):	
	"""
	Forecast the mean and standard deviation of mass given the mean and standard deviation of the radius.

	Parameters
	---------------
	mean: float
		Mean (average) of radius.
	std: float
		Standard deviation of radius.
	unit: string (optional)
		Unit of the radius. Options are 'Earth' and 'Jupiter'.
	sample_size: int (optional)
		Number of radius samples to draw with the mean and std provided.
	grid_size: int (optional)
		Number of grid in the mass axis when sampling mass from radius.
		The more the better results, but slower process.
	Returns
	---------------
	mean: float
		Predicted mean of mass in the input unit.
	std: float
		Predicted standard deviation of mass.
	"""
	# unit
	if unit == 'Earth':
		pass
	elif unit == 'Jupiter':
		mean = mean * rearth2rjup
		std = std * rearth2rjup
	else:
		print "Input unit must be 'Earth' or 'Jupiter'. Using 'Earth' as default."

	# draw samples
	radius = truncnorm.rvs( (0.-mean)/std, np.inf, loc=mean, scale=std, size=sample_size)	
	if classify == 'Yes':
		mass = Rpost2M(radius, 'Earth', grid_size, classify='Yes')
	else:
		mass = Rpost2M(radius, 'Earth', grid_size)

	if mass is None:
		return None

	if unit=='Jupiter':
		mass = mass / mearth2mjup

	return np.mean(mass), np.std(mass)


	