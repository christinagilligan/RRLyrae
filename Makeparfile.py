import sys
import numpy as np

def writepar(modelList):
    f = open('gridpython.par','w')
    f.write('gridsyn                                           \n')
    f.write('terminal       "x11"                              \n')
    f.write('standard_out   "out1"                             \n')
    f.write('smoothed_out   "out3"                             \n')
    f.write('atmosphere     1                                  \n')
    f.write('molecules      2                                  \n')
    f.write('lines          1                                  \n')
    f.write('strong         1                                  \n')
    f.write('flux/int       0                                  \n')
    f.write('damping        1                                  \n')
    f.write('plot           0                                  \n')
    i=0
    for model in modelList:
        i=i+1
        f.write('RUN             '+str(i)+'                                  \n')
        f.write('model_in       "testModels/'+str(model)+'"                            \n')
        f.write('summary_out    "synspec/'+str(model[:-4])+'.rawblul"                        \n')
        f.write('lines_in       "lin4537good"                      \n')
        f.write('stronglines_in "strong4537"                       \n')
        f.write('synlimits                                         \n')
        f.write('  4399.0  4676.0    0.01    0.60                  \n')
        i=i+1
        f.write('RUN            '+str(i)+'                                  \n')
        f.write('model_in       "testModels/'+str(model)+'"                            \n')
        f.write('summary_out    "synspec/'+str(model[:-4])+'.rawyel"                        \n')
        f.write('lines_in       "lin4537good"                      \n')
        f.write('stronglines_in "strong4537"                       \n')
        f.write('synlimits                                         \n')
        f.write('  5149.0  5451.0    0.01    0.60                  \n')
#    f.write('                                                  \n')

    f.close()

modelList=np.loadtxt('modellist', dtype=str)

#smallMod=[]
#smallMod.append(modelList[0])
#smallMod.append(modelList[1])
writepar(modelList)

