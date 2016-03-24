import numpy as np
from scipy.stats import norm, truncnorm

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