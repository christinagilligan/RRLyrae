from astropy.time import Time
import os, sys, re
import numpy as np
import matplotlib.pyplot as plt

#used to figure out when to observe an RR Lyrae when going on an observing run
#calculates maximum and minimum light for 30 cycles starting now

def plot_lightcurve(epoch, period):
    nt=Time.now()

    mjdnow=nt.mjd

    ut=Time(nt, scale='utc')

    epoch=float(epoch)
    period=float(period)


    #figure out how many cycles have occurred since the epoch until now
    numberofperiods=np.round((mjdnow-epoch)/period,0)


    #which ones to actually plot when maximum occured
    whichmaximums=np.arange(numberofperiods-1,numberofperiods+30,1)


    maximums=epoch+period*whichmaximums

    #do the same for minimums which is halfway between maximums

    minimums=np.zeros(len(maximums))


    for i in range(len(maximums)-1):
        minimums[i]=(maximums[i]+maximums[i+1])/2

    mjdtimemax=Time(maximums, format='mjd')
    mjdtimemin=Time(minimums, format='mjd')
    utctimemax=mjdtimemax.iso
    utctimemin=mjdtimemin.iso
    print('maximums in UTC')
    #print(maximums, utctimemax)
    print(utctimemax)
    print('minimums in UTC')
    #print(minimums, utctimemin)
    print(utctimemin)

#optional graph of max and min (not really useful I find)
#    ones=np.ones(len(maximums))

#    plt.figure()
#    plt.scatter(maximums,ones)
#    plt.scatter(minimums,-ones)
#    plt.axis((mjdnow-1,mjdnow+7,-1.1,1.1))
#    plt.ticklabel_format(useOffset=False)
#    plt.show()

if __name__ == '__main__':
#just a test case
    plot_lightcurve(33181.404,0.55302814)
