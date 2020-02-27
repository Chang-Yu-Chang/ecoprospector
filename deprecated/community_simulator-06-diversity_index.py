#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

def Shannon(N):
    p = N/np.sum(N)
    p = p[p>0]
    return np.exp(-np.sum(p*np.log(p)))

def Richness(N,thresh=0):
    return np.sum(N>thresh)
