# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 18:04:28 2014

@author: Erik
"""

import numpy as np
from read_snapshot_galaxy import read_snapshot

def gal_split(initfile,ND1,ND2):
    
        
    pos,vel,partid,time_snap,npart = read_snapshot(initfile)
    
    NGas   = npart[0]
    NHalo  = npart[1]
    NDisk  = npart[2]
    NBulge = npart[3]
    NStars = npart[4]
    
    empty_array = np.zeros(np.sum(npart)) 
    
    gal_1 = np.zeros(np.sum(npart)) 
    gal_1[np.int(NHalo + NGas):np.int(NHalo + NGas + ND1)] = 1
    
    gal_2 = np.zeros(np.sum(npart)) 
    gal_2[np.int(NHalo + NGas + ND1):np.int(NHalo + NGas + ND1 + ND2)] = 1
    
    sorted_ids = np.argsort(partid)
    
    gal_1 = gal_1[sorted_ids]
    gal_2 = gal_2[sorted_ids]
    
    return gal_1, gal_2