import numpy as np
import matplotlib.pyplot as plt
import sys, os
from scipy.interpolate import interp1d
import glob
import matplotlib.animation as animation
import re

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx



#performs the gaussian smoothing
def gaussSmooth(fwhm):
    step=0.020
    sigma=fwhm/2
    aa=0.6932/sigma**2
    power=1.0
    p=[]
    for i in range(0,100):
        value=np.exp(-aa*(step*i)**2)
        if value > 0.01:
            p.append(value)
            power = power+2*value
    return p,power


def readfile(filepath,type,stepAdd):
    try:
        if type=='Synth':
            files=open(filepath,'r').read().split('ALL')
            names=['file'+str(num)for num in range(len(files))]
            for num,file in enumerate(files):
                open(names[num],'w').write(file)
                #really only want file 3 out for this program since it's the best model
            model=open('file1','r')
            #only get the values for the synth
            #lines=model.readlines()[4:]
            lines=model.readlines()[3:]
            model.close()
            #get the wavelength information here
            model=open('file1','r')
            wavelengthInfo=model.readlines()[2]
            #wavelengthInfo=model.readlines()[3]
            model.close()
            wavelengthInfo=wavelengthInfo.split()

            #get the wavelength info
            begWave=np.float(wavelengthInfo[0])
            endWave=np.float(wavelengthInfo[1])
            step=np.float(wavelengthInfo[2])
            num=(endWave-begWave)/step+stepAdd #sometimes this is 1?
            wave=np.linspace(begWave,endWave,num)


            tempArray=[]
            for i in lines:
                i=i.replace('-', ' ')
                test=i.split()
                tempArray.append(test)

            tempArray=np.asarray(tempArray)
            tempArray=np.hstack(tempArray)
            #get the spectra values
            spec=np.asarray(tempArray,dtype=float)
            spec=spec*-1+1
            return wave,spec
        if type=='SALT':
            star=np.loadtxt(filepath)
            wave=star[:,0]
            spec=star[:,1]

            return wave,spec
        if type=='Carnegie':
            star=np.loadtxt(filepath)
            wave=star[:,0]
            spec=star[:,1]
            
            return wave,spec
    except:
        print('Not a valid type! Try again.')

#smooth the synthesis model
def smoothModel(spec,fwhm):
    p,power=gaussSmooth(fwhm)
    z=[]
    for value in spec:
        z.append(value)

    for i in range(len(p),len(spec)-len(p)):
        for j in range(0,len(p)):
            z[i]=z[i]+p[j]*(spec[i-j]+spec[i+j])
    normalz=z/power

    return normalz

def velocityShiftStar(wave,spec,velocity):
    deltaLambda=velocity/299792*wave
    shiftedWave=wave+deltaLambda

    return shiftedWave

def computeOC(wavelengths, specSynthInterp, specStarInterp):
    specSynthInterpValue=[]
    specStarInterpValue=[]
    for wave in wavelengths:
        specSynthInterpValue.append(specSynthInterp(wave))
        specStarInterpValue.append(specStarInterp(wave))
    specSynthInterpValue=np.asarray(specSynthInterpValue,dtype=float)
    specStarInterpValue=np.asarray(specStarInterpValue,dtype=float)
    OCs=specStarInterpValue-specSynthInterpValue
    OCabs=np.sum(np.abs(OCs))
    return OCabs

def computechi(wavelengths, specSynthInterp, specStarInterp):
    specSynthInterpValue=[]
    specStarInterpValue=[]
    for wave in wavelengths:
        specSynthInterpValue.append(specSynthInterp(wave))
        specStarInterpValue.append(specStarInterp(wave))
    specSynthInterpValue=np.asarray(specSynthInterpValue,dtype=float)
    specStarInterpValue=np.asarray(specStarInterpValue,dtype=float)
    OC=specStarInterpValue-specSynthInterpValue
    chi=np.sum((OC)**2/specSynthInterpValue)
    return chi



if __name__=='__main__':
    stepAdd=2
    script, starname, spectrafile, spectratype = sys.argv
    waveSynth, specSynth = readfile('SXFor2','Synth',2)
    #waveStar, specStar = readfile('/Users/christinagilligan/RRLyraeProject/RRLyraeSpectra/curve.txt','SALT')
    #waveStar, specStar = readfile('/Users/christinagilligan/RRLyraeProject/RRLyraeSpectra/SXFor_H_uwm_shifted_cont_boxcar.txt_new','SALT')
    #waveStar, specStar = readfile('/Users/christinagilligan/RRLyraeProject/pymoogi-master/example/363sxfortot.txt','Carnegie')
    waveStar, specStar = readfile('/Users/christinagilligan/RRLyraeProject/pymoogi-master/example/'+str(spectrafile),str(spectratype),0)
    SynthInterp=interp1d(waveSynth,specSynth)
    StarInterp=interp1d(waveStar,specStar)


    velocities=np.linspace(-500,500,1001)
    OCvalues=[]
    
    if spectratype == 'SALT':
        wavelengths=np.linspace(4400,4700,200)
    else:
        wavelengths=np.linspace(4200,4300,200)
    fig = plt.figure()
    ims=[]
    for v in velocities:
        shiftedWave=velocityShiftStar(waveSynth, specSynth, v)
        ShiftedSynthInterp=interp1d(shiftedWave,specSynth)
        OCvalues.append(computeOC(wavelengths, ShiftedSynthInterp, StarInterp))
        #if v  == 0:
            #plt.figure()
            #plt.plot(shiftedWave,specSynth)
            #plt.plot(waveStar, specStar)
            #plt.xlim(4320,4370)
            #plt.ylim(0,1.2)
            #plt.show()
            #plt.title('Velocity: '+str(round(v,2)))
            #plt.savefig(str(v)+'.jpg')

    vel=velocities[np.flatnonzero(OCvalues==np.min(OCvalues))]
    shiftedWave=velocityShiftStar(waveSynth, specSynth, vel)
    ShiftedSynthInterp=interp1d(shiftedWave,specSynth)


    #now that we have the velocity shift figured out, it's time to find the right smoothing


    fwhms=np.linspace(0.05,1,20)
    OCvalues=[]
    for fwhm in fwhms:
        normalz=smoothModel(specSynth,fwhm)
        smoothed=interp1d(shiftedWave, normalz)
        OCvalues.append(computeOC(wavelengths, smoothed, StarInterp))
#        plt.figure()
#        plt.plot(shiftedWave,smoothed(shiftedWave))
#        plt.plot(waveStar, specStar)
#        plt.xlim(4320,4370)
#        plt.ylim(0,1.2)
        #plt.show()
#        plt.title('FWHM: '+str(round(fwhm,2)))
#        plt.savefig(str(fwhm)+'.jpg')

#    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
#    plt.show()

    fwhm=fwhms[np.flatnonzero(OCvalues==np.min(OCvalues))]
    normalz=smoothModel(specSynth,fwhm)
    smoothed=interp1d(shiftedWave, normalz)

    plt.figure()
    plt.plot(shiftedWave,smoothed(shiftedWave))
    plt.plot(waveStar, specStar)
    plt.xlim(4320,4370)
    plt.ylim(0,1.2)
#    plt.savefig('final'+'.jpg')
    plt.show()

    #write the output to my file
    f = open('velandfwhm.txt','a')
    f.write('\n'+str(starname)+','+str(vel[0])+','+str(fwhm[0]))
    f.close()

