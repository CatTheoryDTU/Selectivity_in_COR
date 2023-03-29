import os,sys
import numpy as np
from ase.io import read,write
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from general_tools import quad_fun

home=os.getcwd()

dat=[]
for files in os.listdir():
    if files[-5:] != '.traj': continue

    atoms=read(files)
    pot,ne,E = atoms.calc.results['electrode_potential'],atoms.calc.results['ne'],atoms.calc.results['energy']
    dat.append([pot,ne,E])


    plt.plot(ne,E,'o')

dat=np.array(dat)
c,d=curve_fit(quad_fun,dat[:,1],dat[:,2])
xvals=np.linspace(-1.1,0.5,30)
print(c)


plt.plot(xvals,xvals[:]**2*c[0]+xvals[:]*c[1]+c[2])
plt.xlim([-0.5,1])
plt.show()


