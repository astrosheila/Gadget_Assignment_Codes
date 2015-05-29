# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 02:43:46 2014

@author: Erik
"""


import numpy as np
import matplotlib.pyplot as plt
from read_snapshot_galaxy import read_snapshot
from scipy.special import erf
from scipy.optimize import brentq
from galaxy_split_rot import gal_split
import pylab
pylab.ion()

num_snap = 5
base_name = 'flyby'
asym_arr = np.zeros(num_snap)
time_arr = np.zeros(num_snap)

gal_1, gal_2 = gal_split('flyby_000',20000,4000)

for num in range(num_snap + 1):
    fname = base_name + '_%03d' % num
    print 'Reading',fname,'...'
    
    #Read in all relevant parameters    
    pos,vel,partid,time_snap,npart = read_snapshot(fname)
    time_snap = np.round(time_snap,1)
    time_arr[num-1] = time_snap
    print 'Time of snapshot ->',time_snap        
    
    sorted_inds = np.argsort(partid)
            
    xdisk = pos[sorted_inds, 0]*gal_1
    xdisk = xdisk[np.nonzero(xdisk)]
    ydisk = pos[sorted_inds, 1]*gal_1
    ydisk = ydisk[np.nonzero(ydisk)]
    zdisk = pos[sorted_inds, 2]*gal_1
    zdisk = zdisk[np.nonzero(zdisk)]        
    
    xdisk_2 = pos[sorted_inds, 0]*gal_2
    xdisk_2 = xdisk_2[np.nonzero(xdisk_2)]
    ydisk_2 = pos[sorted_inds, 1]*gal_2
    ydisk_2 = ydisk_2[np.nonzero(ydisk_2)]
    zdisk_2 = pos[sorted_inds, 2]*gal_2
    zdisk_2 = zdisk_2[np.nonzero(zdisk_2)]
    
    part_a = 1
    if part_a:
        # Here we use previously defined statistics to bound our selection for particles within 
        # the disk to simulate our "slit" placement along the major axis of the galaxy.
        # Typically, this parameter would be chosen to match your observational data. Here I have 
        # chosen it (somewhat arbitrarily) to be the width of all the stars within one scale height of 
        # the disk at the beginning of the simulation. Changing the slitwidth parameter will change the
        # number of stars you sample for your rotation curve. 
            
        slitwidth = 2.5
        z_center = np.median(zdisk)
        sel = np.argwhere(np.abs(zdisk - z_center) < slitwidth/2.0)
        xdisk_sel = xdisk[sel]
        zdisk_sel = zdisk[sel]
    
        # Here we plot the positions of the galaxies  
        plt.figure(1)    
        plt.clf()
        plt.plot(xdisk,ydisk, 'b.', markersize = 0.5, alpha=0.5)
        plt.plot(xdisk_2,ydisk_2, 'c.', markersize = 0.5, alpha=0.5)
        plt.title('Time = %s' % time_snap)
        plt.xlim(-150,150)
        plt.ylim(-150,150)
        plt.xlabel('X Projection (kpc)')
        plt.ylabel('Y Projection (kpc)')

        # In a separate figure, plot the x-z projection of only the primary galaxy and 
        # overplot the position of your slit.
        plt.figure(2)
        plt.clf()
        plt.plot(xdisk,zdisk, 'b.', markersize = 0.5, alpha=0.5)
        plt.plot(xdisk_sel,zdisk_sel, 'r.', markersize = 0.5)
        plt.xlim(-70,70)
        plt.ylim(-70,70)
        plt.xlabel('X Projection (kpc)')
        plt.ylabel('Z Projection (kpc)')
    
    part_2 = 1
    if part_2:

        #Import velocity data and split according to galaxy.
        x_vel_disk = vel[sorted_inds, 0]*gal_1
        x_vel_disk = x_vel_disk[np.nonzero(x_vel_disk)]
        y_vel_disk = vel[sorted_inds, 1]*gal_1
        y_vel_disk = y_vel_disk[np.nonzero(y_vel_disk)]
        z_vel_disk = vel[sorted_inds, 2]*gal_1
        z_vel_disk = z_vel_disk[np.nonzero(z_vel_disk)]   

        #Plot raw data rotation curve
        x_center = np.median(xdisk)
        v_rot = y_vel_disk[sel] 
        rads = xdisk[sel] - x_center
        zeroline = np.zeros(2)
        rads_zero = [min(rads),max(rads)]
    
        plt.figure(3)
        plt.clf()
        plt.plot(rads,v_rot,'r.', alpha=0.4)
        plt.plot(rads_zero,zeroline,'k--')
        plt.ylim(-600,600)
        plt.xlim(-40,40)
        plt.xlabel('Radial Position (kpc)')
        plt.ylabel('Doppler Velocity (km/s)')

        # Our next step after plotting the raw rotation data is to "extract" a rotation
        # curve that you may be more familiar with. The challenge is this: in order to match
        # observations, we must first bin our velocity radially (i.e. split it into small radial
        # bins to match our spatial resolution). Next, we must decide which velocity in each bin
        # is representative of the true radial velocity of the disk in that bin (this can be tough).       # Below we've outlined a popular method to help do this, however items such as binsize, vel_errs,
        # and most probable velocity have been left for you to determine. These parameters would generally
        # get some constraints from the data set you are trying to fit, but will in every case be 
        # left open to interpretation.

    
        #Bin raw rotation curve data
        binsize = 0.5 # in kpc, completely adjustable feel free to change it around and see how it affects the results. 
        rads_binned = np.arange(-20, 20 + binsize,binsize)
        vel_binned=np.zeros(len(rads_binned))
    
        for i in range(len(rads_binned)):
            sel = np.where((rads > (rads_binned[i]-binsize/2.)) & (rads < (rads_binned[i]+binsize/2.0)))
            v_rot_binned = v_rot[sel]
            const_err = 5.0 # Constant value for velocity error to be assigned by user
            v_rot_errs = np.array([const_err]*len(v_rot_binned))
    
        # Here we implement the min/max algorithm outlined by Raychaudhury 1997
        # The value of probability is meant to signify the fraction of velocities in a given bin expected
        # to lie above the adopted velocity's absolute magnitude. 
            probability = 0.1           
            if rads_binned[i] > 0:
                vmin = np.min(v_rot_binned)
                def probmin(vs):
                    pmin = np.sum(np.log((0.5+0.5*erf((v_rot_binned - vs)/(np.sqrt(2.)*v_rot_errs))))- np.log(probability))
                    return pmin
                vel_binned[i] = brentq(probmin, 1.2*vmin, 0.4*vmin)
            if rads_binned[i] < 0:
                vmax = np.max(v_rot_binned)
                def probmax(vs):
                    pmax = np.sum(np.log((0.5+0.5*erf((vs - v_rot_binned)/(np.sqrt(2.)*v_rot_errs))))- np.log(probability))
                    return pmax
                vel_binned[i] = brentq(probmax, 0.4*vmax, 1.2*vmax)
    
        #Overplot exracted rotation curve on raw data        
    
        plt.plot(rads_binned,vel_binned,'b--',linewidth=2.0)
    
    
        #Measure Rotation Curve Asymmetry

        #Select near and far side
        sel_far = np.where(rads_binned > 0.0)
        sel_near  = np.where(rads_binned < 0.0)
    
        rads_far = rads_binned[sel_far]
        vel_rot_far = vel_binned[sel_far]
    
        rads_near  = rads_binned[sel_near]
        vel_rot_near = vel_binned[sel_near]
    
        # Invert the negative side of the rotation curve to overplot on top of the positive 
        rads_near = np.abs(rads_near)
        vel_rot_far = np.abs(vel_rot_far)
    
        plt.figure(4)
        plt.clf()
        plt.plot(rads_far,vel_rot_far,'b.')
        plt.plot(rads_near,vel_rot_near,'r.')
        plt.ylim(0,600)
        plt.xlim(0,20.5)
        plt.xlabel('Radial Position (kpc)')
        plt.ylabel('Magnitude of Dopller Velocity (km/s)')
    
        # Measure asymmetry based on method of Daly et. al 2001
        vel_rot_near_rev = vel_rot_near[::-1]
        vel_diff = np.abs((vel_rot_far - vel_rot_near_rev))
        vel_avg = (1./2.)*(vel_rot_far + vel_rot_near_rev)
        
        asym = np.sum(vel_diff)/np.sum(vel_avg)
        asym_arr[num-1] = asym

        plt.figure(5)
        plt.clf()
        plt.plot(time_arr,asym_arr,'bo')
        plt.xlim(-0.2, 1.2)
        plt.ylim(0.0, 0.6)
        plt.xlabel('Time (Gyr)')
        plt.ylabel('Fractional Area Between Curves')
    

    raw_input("Press Enter to move to next time-step...")



