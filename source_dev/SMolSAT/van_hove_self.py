#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:58:21 2018

@author: bruce
"""

from ._SMolSAT import *

import numpy as np
import matplotlib.pyplot as plt
import mpltex


class vhs:
    def __init__(self,system=None,trajs=None,listname=None,out=None,nbins=None,r_cut=None):
        self.system=system
        self.trajs=trajs
        self.listname=listname
        self.out=out
        self.nbins=nbins
        self.r_cut=r_cut

        self.analysis=Van_Hove_Self(system,nbins,r_cut)
        self.analysis.run(trajs,listname)
        self.analysis.write(out)
        self.read()

    def get(self):
        return self.data["ngp"]
    
    def get_t(self):
        return self.data[:,0]

    def read(self):
        with open(self.out) as f:
            line_header = (line.strip() for line in f if len(line.strip().split())==self.nbins)
            self.r_bins=np.loadtxt(line_header, delimiter='\t')
            
        with open(self.out) as f:
            lines = (line.strip() for line in f if len(line.strip().split())==self.nbins+1)
            self.data = np.loadtxt(lines, delimiter='\t')

    def plot(self,file=None):
        
        linestyle=mpltex.linestyles(hollow_styles=[True],lines=["-"],markers=[])
        #create the plot object
        fig, ax = plt.subplots(nrows=1)

        ax.plot(self.data["t"],self.data["ngp"],**next(linestyle),label=None,zorder=3)
        
        #set the plot parameters
        ax.legend(loc="best",ncol=1)
        #ax.set_yscale('log')
        ax.set_xscale('log')
        #ax.set_xlim(90,1100)
        #ax.set_ylim(1e-7,1e-3)
        ax.set_xlabel(r"time")
        ax.set_ylabel(r"ngp")
        fig.tight_layout(pad=0.1) 
        if file!=None:
            fig.savefig(file)
