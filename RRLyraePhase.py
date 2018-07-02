from astropy.time import Time
import os, sys, re
import numpy as np

#time has to be in the iso format: YYYY-MM-DD HH:MM:SS
#period in days
#epoch in MJD
#outputs a numpy array so it can be easily used to use numpy array inputs
def calculate_phase(time, period, eph_epoch):
	time= [str(time)]
#creates an astropy time object
	nt=Time(time, format='iso',scale='utc')
#converts are ISO UTC time to a MJD
	mjd=nt.mjd

	phase = np.mod((mjd - eph_epoch)/period, 1)
	return(phase)
	
if __name__ == '__main__':
#just a test case
	calculate_phase('2018-01-01', 0.5, 245000)
