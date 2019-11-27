import numpy as np
import glob
import sys, os
from VelocityandSmoothing import *
from scipy.interpolate import interp1d
from multiprocessing import Pool
from functools import partial

#want to find the smoothing and velocity values for our particular star
#i will then use these values on our models instead of running MOOG for all the various ones it could be

modelpath = '/dartfs-hpc/rc/home/g/f002bwg/rrlyraespectra/ModelGrids/synspec/*K/t*all'

script, starname, spectrafile, spectratype = sys.argv

velfwhm=np.loadtxt('/dartfs-hpc/rc/home/g/f002bwg/rrlyraespectra/velandfwhm.txt',dtype=str,delimiter=',')

star=velfwhm[:,0]

starIndex=np.flatnonzero(star==str(starname))

vel=np.float64(velfwhm[starIndex,1])
fwhm=np.float64(velfwhm[starIndex,2])


OCvalues=[]
chisquare=[]

#make a dictionary to be able to keep track of the model parameters
modelDict= {'modelFile': 'test', 'parameters':{'temp': None, 'logg': None, 'feh': None, 'micro': None, 'alpha': None, 'OC': None}}

modelList = glob.glob(modelpath)
print(str(spectrafile))
print(str(spectratype))
waveStar, specStar = readfile('/dartfs-hpc/rc/home/g/f002bwg/rrlyraespectra/'+str(spectrafile),str(spectratype),0)

def processData(model, waveStar, specStar, vel, fwhm):
    
    try:
        waveSynth, specSynth = readfile(model,'Synth',2)

        SynthInterp=interp1d(waveSynth,specSynth)
        StarInterp=interp1d(waveStar,specStar)

        shiftedWave=velocityShiftStar(waveSynth, specSynth, vel)
        ShiftedSynthInterp=interp1d(shiftedWave,specSynth)

        normalz=smoothModel(specSynth,fwhm)
        smoothed=interp1d(shiftedWave, normalz)

        if spectratype == 'SALT':
            wavelengths=np.linspace(4910,4960,100)

        if spectratype == 'Carnegie':
            wavelengths=np.linspace(4742,5148,2000)
        else:
            wavelengths=np.linspace(4742,5148,2000)

        OCvalues.append(computeOC(wavelengths, smoothed, StarInterp))
        OCvalue = computeOC(wavelengths, smoothed, StarInterp)

        newstr = ''.join((ch if ch in '0123456789.' else ' ') for ch in model)

        listOfNumbers = [float(i) for i in newstr.split()]
        chisquare.append(computechi(wavelengths, smoothed, StarInterp))
        chisquarevalue = computechi(wavelengths, smoothed, StarInterp)
        modelDict[model]={'temp':listOfNumbers[0],'logg':listOfNumbers[1],'feh':listOfNumbers[2],'micro':listOfNumbers[3],'alpha':listOfNumbers[4],'OC':computeOC(wavelengths, smoothed, StarInterp)}

        # Write results to file
        return "\n {0} {1} {2}".format(str(model), str(OCvalue), str(chisquarevalue))

    except Exception as e:
        return "\n {0} failed".format(str(e))


pool = Pool(processes=16)
modelFunc = partial(processData, waveStar=waveStar, specStar=specStar, vel=vel, fwhm=fwhm)
smallMod=[]
smallMod=modelList[0:100]
f = open('monteCarloResults'+str(starname)+'Chris.txt','a')
for result in pool.map(modelFunc, smallMod):
    f.write(result)
