# -*- coding: utf-8 -*-
"""
Created on Fri Mar 07 13:43:18 2014

@author: Erik
"""

import numpy as np
import matplotlib.pyplot as plt
from read_snapshot_galaxy import read_snapshot
#from galaxy_split_disk import gal_split
from galaxy_split_disk import gal_split
import pylab
pylab.ion()
pi=np.pi

num_snap = 11
base_name_mer = 'merger/merger'
base_name_sin = 'single/single'

gal_1_mer, gal_2_mer = gal_split('merger/merger_000',20000,20000)
gal_1_sin, gal_2_sin = gal_split('single/single_000',20000,0)

# The function gal_split takes 3 arguments, the name of the first
# file in the directory, the number of particles in the primary galaxy disk
# and the number of particles in the satellite galaxy disk

z_scale_arr_mer = np.zeros(num_snap+1)
z_scale_arr_sin = np.zeros(num_snap+1)
time_arr = np.zeros(num_snap+1)

for num in range(num_snap):
    #Read in all relevant parameters for the merger

    fname = base_name_mer + '_%03d' % num
    print 'Reading',fname,'...'
        
    pos,vel,partid,time_snap,npart = read_snapshot(fname)
    time_snap = np.round(time_snap,2)
    
    time_arr[num] = time_snap
    
    sorted_inds = np.argsort(partid)
           
    xdisk = pos[sorted_inds, 0]
    xdisk = xdisk[np.where(gal_1_mer == 'True')]
    ydisk = pos[sorted_inds, 1]
    ydisk = ydisk[np.where(gal_1_mer == 'True')]
    zdisk = pos[sorted_inds, 2]
    zdisk = zdisk[np.where(gal_1_mer == 'True')]        
    
    xdisk_2 = pos[sorted_inds, 0]
    xdisk_2 = xdisk_2[np.where(gal_2_mer == 'True')]
    ydisk_2 = pos[sorted_inds, 1]
    ydisk_2 = ydisk_2[np.where(gal_2_mer == 'True')]
    zdisk_2 = pos[sorted_inds, 2]
    zdisk_2 = zdisk_2[np.where(gal_2_mer == 'True')]

    #Here we shift our points such that the primary disk is always at the center

    x_center = np.median(xdisk)
    y_center = np.median(ydisk)
    z_center = np.median(zdisk)
    
    xdisk = xdisk - x_center
    xdisk_2 = xdisk_2 - x_center
    
    ydisk = ydisk- y_center
    ydisk_2 = ydisk_2 - y_center

    zdisk = zdisk - z_center
    zdisk_2 = zdisk_2 - z_center
    
    plt.figure(1)
    plt.clf()
    plt.plot(xdisk,ydisk, 'b.', markersize = 0.5, alpha = 0.5)
    plt.plot(xdisk_2,ydisk_2,'r.',markersize = 0.5, alpha = 0.5)
    plt.xlim(-100,100)
    plt.ylim(-100,100)
    plt.title('Time = %s' %time_snap)
    plt.xlabel('X Proj (kpc)')
    plt.ylabel('Y Proj (kpc)')
    
    # Measure disk thickening in terms of the scale height of the disk
    
    z_sorted = zdisk[np.argsort(zdisk)]
    z_scale1 = np.abs(z_sorted[np.int(0.84*len(z_sorted))])
    z_scale2 = np.abs(z_sorted[np.int(0.16*len(z_sorted))])

    zrms = 0.5*(z_scale1+z_scale2)
    z_scale_arr_mer[num] = zrms #z_scale

    #Read in all relevant parameters for the single galaxy

    fname = base_name_sin + '_%03d' % num
    print 'Reading',fname,'...'

    pos,vel,partid,time_snap,npart = read_snapshot(fname)
    time_snap = np.round(time_snap,2)
    print 'Time of snapshot ->',time_snap
    sorted_inds = np.argsort(partid)
            
    xdisk = pos[sorted_inds, 0]
    xdisk = xdisk[np.where(gal_1_sin == 'True')]
    ydisk = pos[sorted_inds, 1]
    ydisk = ydisk[np.where(gal_1_sin == 'True')]
    zdisk = pos[sorted_inds, 2]
    zdisk = zdisk[np.where(gal_1_sin == 'True')]        
    
    # Here we shift our points such that the primary disk is always at the center

    x_center = np.median(xdisk)
    y_center = np.median(ydisk)
    z_center = np.median(zdisk)
    
    xdisk = xdisk - x_center
    
    ydisk = ydisk- y_center

    zdisk = zdisk - z_center
        
    plt.figure(2)
    plt.clf()
    plt.plot(xdisk,ydisk, 'b.', markersize = 0.5, alpha = 0.5)
    plt.xlim(-100,100)
    plt.ylim(-100,100)
    plt.title('Time = %s' %time_snap)
    plt.xlabel('X Proj (kpc)')
    plt.ylabel('Y Proj (kpc)')
    
    # Measure disk thickening in terms of the scale height of the disk
    
    z_sorted = zdisk[np.argsort(zdisk)]
    z_scale1 = np.abs(z_sorted[np.int(0.84*len(z_sorted))])
    z_scale2 = np.abs(z_sorted[np.int(0.16*len(z_sorted))])

    zrms = 0.5*(z_scale1+z_scale2)
    z_scale_arr_sin[num] = zrms 

    domore = 1
    if domore:
        plt.figure(3)
        plt.clf()
        plt.plot(time_arr,z_scale_arr_mer,'r.')
        plt.plot(time_arr,z_scale_arr_sin,'b.')
        plt.xlim(0.0,4.1) 
        plt.xlabel('Time (Gyr)')
        plt.ylabel('Disk Scale Height (kpc)')
    
    raw_input("Press Enter to move to next time step...")

