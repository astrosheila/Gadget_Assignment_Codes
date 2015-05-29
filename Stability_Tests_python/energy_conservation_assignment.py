# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 13:06:46 2014

@author: Erik
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pylab
pylab.ion()

# Adding the lowres filename and another color to overplot with in these lists will automatically read-in and overplot the new data
filenames= ['energy_mediumres.txt'] 
colorlist=['b'] 

plt.figure(1)
plt.clf()

for i in range(len(filenames)):    
    data = np.loadtxt(filenames[i])
    time = data[:,0]
    th_energy = data[:,1]
    pot_energy = data[:,2]
    kin_energy = data[:,3]
    
    tot_energy = kin_energy + pot_energy
    tot_energy_norm = tot_energy/tot_energy[100]

    plt.plot(time,tot_energy_norm,color=colorlist[i])
    plt.xlabel('Time (Gyr)')
    plt.ylabel('Normalize Total Energy')

plt.show()
