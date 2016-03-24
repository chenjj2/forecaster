forecaster
==========

Predict mass (radius) given radius (mass) measurements.

See xxx paper for details. 

If you use it, please cite xxx paper (link).



Usage
-----

Check demo.ipynb for more details.


A simple example:

	import numpy as np
	import mr_forecast as mr
	
	# predict the mean and std of radius given mass measurements

	Rmean, Rstd = mr.Mstat2R(mean=1.0, std=0.1, unit='Earth', sample_size=100)





