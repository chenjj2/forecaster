forecaster
==========

Forecaster uses a probabilistic mass-radius relation as the underlying model.

It can forecast mass (radius) given radius (mass) measurements.

The conversion includes three sources of uncertainties, from measurement, model fitting (MCMC process), and intrinsic dispersion in radius.

See arXiv:1603.08614 for details. 

If you use it, please cite it.



Usage
-----

Check demo.ipynb for more details.


A simple example:

	import numpy as np
	import mr_forecast as mr
	
	# predict the mean and std of radius given mass measurements

	Rmedian, Rplus, Rminus = mr.Mstat2R(mean=1.0, std=0.1, unit='Earth', sample_size=100)

A simple interactive example:
    
    print '=== Forecaster ==='
    print ' '
    print 'Example: Radius-to-Mass Conversion (without posteriors)'
    print ' '
    print 'Radius = A +/- B [Earth units]'
    mean = float(raw_input("Enter A: "))
    std = float(raw_input("Enter B: "))

    # predict the mean and std of radius given mass measurements
    Mmedian, Mplus, Mminus = mr.Rstat2M(mean, std, unit='Earth', sample_size=1e3, grid_size=1e3)
    print ' '
    print 'Mass = ',Mmedian,'+',Mplus,'-',Mminus,' M_earth'



