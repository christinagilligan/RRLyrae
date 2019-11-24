import numpy as np
import glob
import sys, os
from VelocityandSmoothing import *
#import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
#want to find the smoothing and velocity values for our particular star
#i will then use these values on our models instead of running MOOG for all the various ones it could be

script, starname, spectrafile, spectratype, inputmodel = sys.argv
inputmodel=np.int(inputmodel)
velfwhm=np.loadtxt('/dartfs-hpc/rc/home/g/f002bwg/rrlyraespectra/velandfwhm.txt',dtype=str,delimiter=',')
star=velfwhm[:,0]

#starIndex=np.flatnonzero(star=='SXFor')
starIndex=np.flatnonzero(star==str(starname))

vel=np.float64(velfwhm[starIndex,1])
fwhm=np.float64(velfwhm[starIndex,2])

modelpath='/dartfs-hpc/rc/home/g/f002bwg/rrlyraespectra/ModelGrids/synspec/*K/t*all'
modelList=glob.glob(modelpath)

#waveStar, specStar = readfile('/Users/christinagilligan/RRLyraeProject/RRLyraeSpectra/SXFor_H_uwm_shifted_cont_boxcar.txt_new','SALT')
#waveStar, specStar = readfile('/Users/christinagilligan/RRLyraeProject/pymoogi-master/example/363sxfortot.txt','Carnegie')
waveStar, specStar = readfile('/dartfs-hpc/rc/home/g/f002bwg/rrlyraespectra/'+str(spectrafile),str(spectratype),0)
#waveStar, specStar = readfile('/Users/christinagilligan/RRLyraeProject/RRLyraeSpectra/curve.txt','SALT')

#to delete later

model=modelList[inputmodel]

OCvalues=[]
chisquare=[]

#make a dictionary to be able to keep track of the model parameters

modelDict= {'modelFile': 'test', 'parameters':{'temp': None, 'logg': None, 'feh': None, 'micro': None, 'alpha': None, 'OC': None}, 'chisquare': None}

#for model in smallMod:
#    f = open('temp.txt','a')
#    f.write('\n'+'Running model '+str(i)+' of '+str(len(modelList)))
#    f.close()
#    print('Running model '+str(i)+' of '+str(len(modelList)))
waveSynth, specSynth = readfile(model,'Synth',2)
#waveStar, specStar = readfile('/Users/christinagilligan/RRLyraeProject/RRLyraeSpectra/SXFor_H_uwm_shifted_cont_boxcar.txt_new','SALT')
waveStar, specStar = readfile('/dartfs-hpc/rc/home/g/f002bwg/rrlyraespectra/'+str(spectrafile),str(spectratype),0)
#waveStar, specStar = readfile('/Users/christinagilligan/RRLyraeProject/RRLyraeSpectra/curve.txt','SALT')
SynthInterp=interp1d(waveSynth,specSynth)
StarInterp=interp1d(waveStar,specStar)
shiftedWave=velocityShiftStar(waveSynth, specSynth, vel)
ShiftedSynthInterp=interp1d(shiftedWave,specSynth)
normalz=smoothModel(specSynth,fwhm)
smoothed=interp1d(shiftedWave, normalz)
if spectratype == 'SALT':
    #wavelengths=np.linspace(4750,5140,2000)
    wavelengths=np.linspace(4910,4960,100)
if spectratype == 'Carnegie':
    wavelengths=np.linspace(4742,5148,2000)
else:
    wavelengths=np.linspace(4742,5148,2000)
#wavelengths=np.linspace(4742,5148,2000)
#wavelengths=np.linspace(4444,4650,600)
#wavelengths=np.linspace(5149,5451,600)
OCvalues.append(computeOC(wavelengths, smoothed, StarInterp))
newstr = ''.join((ch if ch in '0123456789.' else ' ') for ch in model)
listOfNumbers = [float(i) for i in newstr.split()]
chisquare=computechi(wavelengths, smoothed, StarInterp)

modelArray=[]
modelDict[model]={'temp':listOfNumbers[0],'logg':listOfNumbers[1],'feh':listOfNumbers[2],'micro':listOfNumbers[3],'alpha':listOfNumbers[4],'OC':computeOC(wavelengths, smoothed, StarInterp), 'chisquare':chisquare}
modelArray.append(listOfNumbers)
#    plt.figure()
#    plt.plot(shiftedWave,smoothed(shiftedWave))
#    plt.plot(waveStar, specStar)
#    plt.xlim(4910,4960)
#    plt.ylim(0.2,1)
#    plt.show()
#    if i % 100 ==0 or i==1:
            #plt.figure()
            #plt.plot(shiftedWave,smoothed(shiftedWave))
            #plt.plot(waveStar, specStar)
            #plt.xlim(4575,4600)
            #plt.ylim(0.4,1.2)
            #plt.show()
            #plt.title('Temp: '+str(listOfNumbers[0])+ ' Log g: '+str(listOfNumbers[1])+' Fe/H: -'+str(listOfNumbers[2])+' Microturb: '+str(listOfNumbers[3]))
            #plt.savefig(str(i)+'model.jpg')
modelArray=np.asarray(modelArray,dtype=int)
modelArray=np.append(modelArray,computeOC(wavelengths, smoothed, StarInterp))
modelArray=np.append(modelArray,chisquare)
f = open('monteCarlo'+str(starname)+'.'+str(inputmodel),'a')
f.write('\n '+str(inputmodel)+' '+str(model)+' '+str(computeOC(wavelengths, smoothed, StarInterp))+' '+str(chisquare))
#np.savetxt(f,+modelArray,delimiter=',', newline=" ")
f.close()
print('\n '+str(inputmodel)+' '+str(model)+' '+str(computeOC(wavelengths, smoothed, StarInterp))+' '+str(chisquare))
#now time to plot the best fitting model
#minIndex=np.flatnonzero(chisquare==np.min(chisquare))[0]

#model=modelList[minIndex]

#waveSynth, specSynth = readfile(model,'Synth',2)
#waveStar, specStar = readfile('/Users/christinagilligan/RRLyraeProject/RRLyraeSpectra/SXFor_H_uwm_shifted_cont_boxcar.txt_new','SALT')
#waveStar, specStar = readfile('/Users/christinagilligan/RRLyraeProject/pymoogi-master/example/363sxfortot.txt','Carnegie')
#waveStar, specStar = readfile('/dartfs-hpc/rc/home/g/f002bwg/rrlyraespectra/'+str(spectrafile),str(spectratype),0)
#waveStar, specStar = readfile('/Users/christinagilligan/RRLyraeProject/RRLyraeSpectra/curve.txt','SALT')
#SynthInterp=interp1d(waveSynth,specSynth)
#StarInterp=interp1d(waveStar,specStar)
#shiftedWave=velocityShiftStar(waveSynth, specSynth, vel)
#ShiftedSynthInterp=interp1d(shiftedWave,specSynth)
#normalz=smoothModel(specSynth,fwhm)
#smoothed=interp1d(shiftedWave, normalz)
#wavelengths=np.linspace(4444,4650,600)
#wavelengths=np.linspace(5149,5451,600)
#plt.figure()
#plt.plot(shiftedWave,smoothed(shiftedWave))
#plt.plot(waveStar, specStar)
#plt.xlim(4444,4650)
#plt.ylim(0,1.2)
#plt.xlim(4575,4600)
#plt.ylim(0.4,1.2)
#plt.title(modelDict[model])
#plt.show()
#print('Temperature: '+str(modelDict[model]['temp']))
#print('log g: '+str(modelDict[model]['logg']))
#print('FeH: '+str(modelDict[model]['feh']))
#print('Micro: '+str(modelDict[model]['micro']))
#print('Alpha: '+str(modelDict[model]['alpha']))
#print('OC: '+str(modelDict[model]['OC']))
#print(model)


#minIndex=np.flatnonzero(chisquare<np.min(chisquare)+0.05)

#f = open('monteCarloResults.txt','a')
#f.write('\n'+str(starname))
#for i in range(0,len(modelList)):
#    f.write('\n'+str(model[i])+' '+str(OCvalues[i])+' '+str(chisquare[i]))
#f.close()
