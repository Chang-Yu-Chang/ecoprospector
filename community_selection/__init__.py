#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import copy
import community_selection
from multiprocessing import Pool
from functools import partial

from community_simulator import Community

class Metacommunity(Community):
    """
    Inherited object from community-simulator package. 
    
    Changes:
    
    - Passage are Possion distributed
    
    """
    def Passage(self,f,scale=None,refresh_resource=True):
        """
        Transfer cells to a fresh plate.
        
        f = matrix specifying fraction of each old well (column) to transfer 
            to each new well (row)
            
        scale = option for using a different scale factor from the one defined 
            for the plate on initialization.
            
        refresh_resource says whether the new plate comes supplied with fresh 
            media. The resource concentrations in the media are assumed to be
            the same as the initial resource concentrations from the first plate.
            The "Reset" method can be used to adjust these concentrations.
        """
        #HOUSEKEEPING
        if scale == None:
            scale = self.scale #Use scale from initialization by default
        f = np.asarray(f) #Allow for f to be a dataframe
        self.N[self.N<0] = 0 #Remove any negative values that may have crept in
        self.R[self.R<0] = 0
        
        #DEFINE NEW VARIABLES
        N_tot = np.sum(self.N)
        R_tot = np.sum(self.R)
        N = np.zeros(np.shape(self.N))
        
        #MULTINOMIAL SAMPLING
        #(simulate transfering a finite fraction of a discrete collection of cells)
        for k in range(self.n_wells):
            for j in range(self.n_wells):
                if f[k,j] > 0 and N_tot[j] > 0:
                    N[:,k] += np.random.multinomial(np.random.poisson(scale*N_tot[j]*f[k,j]),(self.N/N_tot).values[:,j])*1./scale  
        self.N = pd.DataFrame(N, index = self.N.index, columns = self.N.keys())
        
        #In batch culture, there is no need to do multinomial sampling on the resources,
        #since they are externally replenished before they cause numerical problems
        if refresh_resource:
            self.R = pd.DataFrame(np.dot(self.R,f.T), index = self.R.index, columns = self.R.keys())
            self.R = self.R+self.R0

        #In continuous culture, it is useful to eliminate the resources that are
        #going extinct, to avoid numerical instability
        else:
            R_tot = np.sum(self.R)
            R = np.zeros(np.shape(self.R))
            for k in range(self.n_wells):
                for j in range(self.n_wells):
                    if f[k,j] > 0 and R_tot[j] > 0:
                        R[:,k] += np.random.multinomial(int(scale*R_tot[j]*f[k,j]),(self.R/R_tot).values[:,j])*1./scale
            self.R = pd.DataFrame(R, index = self.R.index, columns = self.R.keys())
