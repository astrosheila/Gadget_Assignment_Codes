# -*- coding: utf-8 -*-
"""
Created on Wed Feb 05 12:16:36 2014

@author: Erik
"""

import numpy as np

def read_snapshot(filename):
    infile = open(filename,"rb")
#    name = str(filename)

    #### Begin header fortran read statement
    dummy = np.fromfile(infile,dtype=np.int32,count=1)

    #The following part should add up to 256 bytes

    npart      = np.fromfile(infile,dtype=np.int32,count=6)
    massarr    = np.fromfile(infile,dtype=np.float64,count=6)
    time_snap   = np.fromfile(infile,dtype=np.float64,count=1) #given as a in read_snapshot.py
    redshift   = np.fromfile(infile,dtype=np.float64,count=1)
    flagsfr    = np.fromfile(infile,dtype=np.int32,count=1)
    flagfeedb  = np.fromfile(infile,dtype=np.int32,count=1)
    npartTotal = np.fromfile(infile,dtype=np.int32,count=6)
    la         = np.fromfile(infile,dtype=np.int32,count=34)

    #### End header fortran read statement
    dummy = np.fromfile(infile,dtype=np.int32,count=1)
    
    #### Begin positions fortran read statement
    dummy = np.fromfile(infile,dtype=np.int32,count=1)

    ### We have particles of both type 1,2 (dm and stars)
    pos = np.fromfile(infile,dtype=np.float32,count=(npart[1]+npart[2])*3).reshape([(npart[1]+npart[2]),3])

    #### End positions fortran read statement
    dummy = np.fromfile(infile,dtype=np.int32,count=1)

    #### Begin velocities fortran read statement
    dummy = np.fromfile(infile,dtype=np.int32,count=1)

    ### Same as positions statement
    vel = np.fromfile(infile,dtype=np.float32,count=(npart[1]+npart[2])*3).reshape([(npart[1]+npart[2]),3])

    #### end velocities fortran read statement
    dummy = np.fromfile(infile,dtype=np.int32,count=1)

    #### Begin id fortran read statement
    dummy = np.fromfile(infile,dtype=np.int32,count=1)

    ####BE CAREFUL HERE!!! If dtype is not specified as np.int32, will get weird numbers...
    part_id = np.fromfile(infile,dtype=np.int32,count=(npart[1]+npart[2]))

    #### end id fortran read statement
    dummy = np.fromfile(infile,dtype=np.int32,count=1)
    
    return pos, vel, part_id, time_snap[0], npart




