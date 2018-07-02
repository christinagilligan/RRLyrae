#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 22:28:16 2017

@author: christinagilligan
"""

#creates a obj.list file that is needed for QUBA standards
#
#DOES NOT CONVERT FROM 2MASS MAGS TO IRSF MAGS
#
#usage: catalogCreation.py CS_Eri
#
#works in python 3 and 2.7

import csv
import numpy as np
from sys import argv

#give the name of the star as a an argument for the script
script,objectName = argv


#format
#_RAJ2000	_DEJ2000	RAJ2000	DEJ2000	2MASS	Jmag	e_Jmag	Hmag	e_Hmag	Kmag	e_Kmag	Qflg	Rflg	Bflg	Cflg	Xflg	Aflg


#want AAA for Qflg (12th column)
#want 222 for Aflg (13th column)
#want 111 for Bflg (14th column)

f = open('/Users/christinagilligan/Massimo/pipeline/coordinate_std/infrared/ds9.tsv')
csv_f = csv.reader(f,delimiter='\t')

next(csv_f)
catalog = np.array(list(csv_f))

goodindexQflg=[]
for x in range(0,len(catalog[:,11])):
    if catalog[x,11]=='AAA':
        goodindexQflg.append(x) 

goodindexAflg=[]
for x in range(0,len(catalog[:,12])):
    if catalog[x,12]=='222':
        goodindexAflg.append(x) 
        
goodindexBflg=[]
for x in range(0,len(catalog[:,13])):
    if catalog[x,13]=='111':
        goodindexBflg.append(x) 

a=np.intersect1d(goodindexQflg, goodindexAflg)
goodIndex=np.intersect1d(a, goodindexBflg)


#add column JH col6-col8
#HK col8-col10

JH=catalog[goodIndex,5].astype(float)-catalog[goodIndex,7].astype(float)
HK=catalog[goodIndex,7].astype(float)-catalog[goodIndex,9].astype(float)

#save columns 1,2,6,JH,JK

with open('/Users/christinagilligan/Massimo/pipeline/coordinate_std/infrared/'+str(objectName)+'.list', 'w') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    for x in range(0,len(goodIndex)):
        writer.writerow([x+1,catalog[goodIndex[x],0], catalog[goodIndex[x],1],catalog[goodIndex[x],5],JH[x],HK[x]])
