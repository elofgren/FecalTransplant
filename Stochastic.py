#!/usr/bin/env python

###################################################
# Dynamic Transmission Model of C. difficile      #
# Fecal Microbiota Transplant Based Interventions #
# Gillespie Direct Method Implementation          #
# Author: Eric Lofgren (Eric.Lofgren@gmail.com)   #
###################################################

# Module Imports
import os
import fileinput # may not need this
import stochpy
import pylab as pl
import matplotlib as mpl # may not need this
import numpy as np
import requests

# Import most recent PML file from Github
pml = requests.get('https://raw.github.com/elofgren/PML/master/cdiff_FT.psc?login=elofgren&token=bf6f8b071be7c69c9b6f7afd8df216a9'
,verify=False)
PMLout = os.path.join(os.getcwd(),'cdiff_FT.psc')
f = open(PMLout,'w')
f.write(pml.content)
f.close()

CDI = stochpy.SSA()
CDI.Model(File='cdiff_FT.psc', dir=os.getcwd())

# Set up simulation parameters
start_time = 0.0
end_time = 8760
n_runs = 50

#######################################
# Baseline Scenario - No Intervention #
#######################################

