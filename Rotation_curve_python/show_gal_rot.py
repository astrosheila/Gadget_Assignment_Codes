# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 02:04:36 2014

@author: Erik
"""

import numpy as np
import matplotlib.pyplot as plt
from read_snapshot_galaxy import read_snapshot
from galaxy_split import gal_split

gal_1, gal_2 = gal_split()
number_snap = 7

plt.figure(1)
plt.ion()
plt.show()


for i in range(number_snap + 1):
    
    num = i
    fname = "snapshot_%03d" % num #frun + "/snapshot_"+exts
    print 'Reading',fname,'...'
    
    pos,vel,partid,time_snap,npart = read_snapshot(fname)
    print 'Time of snapshot ->',np.round(time_snap, decimals=1)
    
    NGas   = npart[0]
    NHalo  = npart[1]
    NDisk  = npart[2]
    NBulge = npart[3]
    NStars = npart[4]
    
    sorted_inds = np.argsort(partid)    

    pos1_x = pos[sorted_inds,0]*gal_1
    pos1_y = pos[sorted_inds,1]*gal_1
    pos1_z = pos[sorted_inds,2]*gal_1

    pos2_x = pos[sorted_inds,0]*gal_2
    pos2_y = pos[sorted_inds,1]*gal_2
    pos2_z = pos[sorted_inds,2]*gal_2

    if NDisk > 0:
        xdisk1=pos1_x[np.nonzero(pos1_x)]
        xdisk2=pos2_x[np.nonzero(pos2_x)]

        ydisk1=pos1_y[np.nonzero(pos1_y)]
        ydisk2=pos2_y[np.nonzero(pos2_y)]
        
        zdisk1=pos1_z[np.nonzero(pos1_z)]
        zdisk2=pos2_z[np.nonzero(pos2_z)]        
    
    plt.clf()
    plt.plot(xdisk1,ydisk1,'b.',ms=1.0, alpha=1.0)
    plt.plot(xdisk2,ydisk2,'r.',ms=1.0, alpha=1.0)
    plt.xlim(-200,200)
    plt.ylim(-200,200)
    plt.title("Time: %s Gyr" %np.round(time_snap,decimals=1))
    plt.pause(0.5)    
