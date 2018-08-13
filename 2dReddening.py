# get current dir
import os
#here = os.path.dirname(__file__)
# use matplotlib without X server
import matplotlib
matplotlib.use('Agg')

import numpy as np
import scipy.interpolate as spi
import h5py
from astropy import units as u
from astropy.coordinates import SkyCoord
from dustmaps.planck import PlanckQuery
from dustmaps.sfd import SFDQuery
import astropy.coordinates as coord
import matplotlib.pyplot as plt
from matplotlib import *
import matplotlib
import time

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx



gaiaid=[]
parallax=[]
parallaxError=[]
period=[]
w1mean=[]
w1error=[]
RA=[]
Dec=[]
G=[]
B=[]
R=[]

gaiaid=[]
parallax=[]
parallaxError=[]
period=[]
w1mean=[]
w1error=[]
RA=[]
Dec=[]
G=[]
B=[]
R=[]

file = open('WISE_v4_w1A.tbl')
for line in file:
    if '#' not in line:
        line=line.strip()
        line=line.split()
        line=np.asarray(line)
        gaiaid.append(line[0])
        parallax.append(line[1])
        parallaxError.append(line[2])
        period.append(line[3])
        w1mean.append(line[4])
        w1error.append(line[5])
        RA.append(line[7])
        Dec.append(line[8])
        G.append(line[9])
        B.append(line[10])
        R.append(line[11])

gaiaid=np.array(gaiaid,dtype=str)
parallax=np.array(parallax,dtype=float)
parallaxError=np.array(parallaxError,dtype=float)
period=np.array(period,dtype=float)
w1mean=np.array(w1mean,dtype=float)
w1error=np.array(w1error,dtype=float)
RA=np.array(RA,dtype=float)
Dec=np.array(Dec,dtype=float)
G=np.array(G,dtype=str)
B=np.array(B,dtype=str)
R=np.array(R,dtype=str)

planck=PlanckQuery()
sfd=SFDQuery()
file = open('2DReddening.dat','w')
file.write('#gaiaid parallax sig_parallax w1mean sig_w1 RA Dec G G_BP G_RP Planckred SFDred newDistance catalog flag'+'\n')
file.close()

for index in range(0, len(gaiaid)-1):
    try:
        coords=SkyCoord(RA[index]*u.deg,Dec[index]*u.deg,frame='icrs')
        reddeningPlanck=planck(coords)
        
        
        file = open('2DReddening.dat','a')
        file.write(str(gaiaid[index])+' '+str(parallax[index])+' '+str(parallaxError[index])+' '+str(period[index])+' '+str(w1mean[index])+' '+str(w1error[index])+' '+str(RA[index])+' '+str(Dec[index])+' '+str(G[index])+' '+str(B[index])+' '+str(R[index])+' '+str(reddeningPlanck)+' '+str(reddeningSFD)+'\n')
        file.close()
    except:
        file = open('2DReddening.dat','a')
        file.write(str(gaiaid[index])+' ERROR \n')
        file.close()

