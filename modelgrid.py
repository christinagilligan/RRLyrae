import numpy as np
import subprocess
import shutil


#make the for loop for temperature, logg, feh, vmicro

begtemp, endtemp, steptemp = 5000, 8000, 100
beglogg, endlogg, steplogg = 1, 2, 0.25
begfeh, endfeh, stepfeh = -3, 0, 0.25
begvmicro, endvmicro, stepvmicro = 2, 4, 0.25

#real steps - ONLY RUN ON DISCOVERY OR POLARIS BECAUSE IT'LL CRASH YOUR MACHINE
#begtemp, endtemp, steptemp = 5900, 7800, 10
#beglogg, endlogg, steplogg = 1.50, 3.0, 0.05
#begfeh, endfeh, stepfeh = -2.50, 0, 0.05
#begvmicro, endvmicro, stepvmicro = 2.5, 4.7, 0.05

temp=np.arange(begtemp, endtemp, step=steptemp)
logg=np.arange(beglogg, endlogg, step=steplogg)
feh=np.arange(begfeh, endfeh, step=stepfeh)
vmicro=np.arange(begvmicro, endvmicro, step=stepvmicro)

#everything=zip(temp,logg,feh,vmicro)

for Temp in temp:
    for Logg in logg:
        for Feh in feh:
                for Vmicro in vmicro:
                    modelName='t'+str(Temp)+'g'+'{:.2f}'.format(Logg)+'mm'+'{:.2f}'.format(Feh)+'v'+str(Vmicro)+'alp'+str('000')+'.mod'
                    query1=str(Temp)+' '+str(Logg)+' '+str(Feh)+' '+str(Vmicro)
                    query2='AODFNEW'
                    concat_query="{}\n{}".format(query1, query2)
                    p=subprocess.Popen(['../modprog/makekurandy/makekurucz3.e'],stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                    p.communicate(input=concat_query.encode('utf-8'))[0]
                    shutil.copy("FINALMODEL", "./testModels/"+str(modelName))
#example query
#query1='7000 2 -1 3'
#query2='AODFNEW'
#concat_query="{}\n{}".format(query1, query2)
#p=subprocess.Popen(['./makekurucz3.e'],stdin=subprocess.PIPE, stdout=subprocess.PIPE)
#p.communicate(input=concat_query.encode('utf-8'))[0]


#then need to move this model to the normal naming convention that Chris uses
#shutil.move("path/to/current/file.foo", "path/to/new/destination/for/file.foo")
