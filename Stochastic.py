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
pml = requests.get(
'https://raw.github.com/elofgren/PML/master/cdiff_FT.psc?login=elofgren&token=bf6f8b071be7c69c9b6f7afd8df216a9'
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
n_runs = 1000

# Headers for output file
header = "Treated, Incident, Recur, Level"

#######################################
# Baseline Scenario - No Intervention #
#######################################

baselineresults = np.zeros([n_runs, 4])

def BaselineRun(model,iteration):
	model.DoStochSim()
	data = model.data_stochsim.getSpecies()
	Treated = data[-1,8]
	Incident = data[-1,10]
	Recur = data[-1,14]
	baselineresults[iteration,0] = Treated
	baselineresults[iteration,1] = Incident
	baselineresults[iteration,2] = Recur
	baselineresults[iteration,3] = 0

for k in range(n_runs):
	print "Baseline Iteration %i of %i" % (k+1,n_runs)
	BaselineRun(CDI,k)
	
np.savetxt('baseline.csv',baselineresults,delimiter=',',header=header,comments='')
del baselineresults

########################################
# FMT Intervention - CDI Patients Only #
########################################
cdiresults = np.zeros([n_runs, 4])

def CDIRun(model,percent,iteration):
	model.DoStochSim()
	model.ChangeParameter('chi_postCDI',percent)	
	data = model.data_stochsim.getSpecies()
	Treated = data[-1,8]
	Incident = data[-1,10]
	Recur = data[-1,14]
	cdiresults[iteration,0] = Treated
	cdiresults[iteration,1] = Incident
	cdiresults[iteration,2] = Recur
	cdiresults[iteration,3] = percent
	
for j in range(n_runs):
	print "CDI Only (0.20) Iteration %i of %i" % (j+1,n_runs)
	CDIRun(CDI,0.20,j)

np.savetxt('cdionlyresults20.csv',cdiresults,delimiter=',',header=header,comments='')

for j in range(n_runs):
	print "CDI Only (0.40) Iteration %i of %i" % (j+1,n_runs)
	CDIRun(CDI,0.40,j)

np.savetxt('cdionlyresults40.csv',cdiresults,delimiter=',',header=header,comments='')

for j in range(n_runs):
	print "CDI Only (0.60) Iteration %i of %i" % (j+1,n_runs)
	CDIRun(CDI,0.60,j)

np.savetxt('cdionlyresults60.csv',cdiresults,delimiter=',',header=header,comments='')

for j in range(n_runs):
	print "CDI Only (0.80) Iteration %i of %i" % (j+1,n_runs)
	CDIRun(CDI,0.80,j)

np.savetxt('cdionlyresults80.csv',cdiresults,delimiter=',',header=header,comments='')

for j in range(n_runs):
	print "CDI Only (1.00) Iteration %i of %i" % (j+1,n_runs)
	CDIRun(CDI,1.00,j)

np.savetxt('cdionlyresults100.csv',cdiresults,delimiter=',',header=header,comments='')
del cdiresults

########################################
# FMT Intervention - ABX Patients Only #
########################################
abxresults = np.zeros([n_runs, 4])

def ABXRun(model,percent,iteration):
	model.DoStochSim()
	model.ChangeParameter('chi',percent)
	model.ChangeParameter('phi_c',0.001102321)	
	data = model.data_stochsim.getSpecies()
	Treated = data[-1,8]
	Incident = data[-1,10]
	Recur = data[-1,14]
	abxresults[iteration,0] = Treated
	abxresults[iteration,1] = Incident
	abxresults[iteration,2] = Recur
	abxresults[iteration,3] = percent
	
for l in range(n_runs):
	print "ABX Only (0.20) Iteration %i of %i" % (l+1,n_runs)
	ABXRun(CDI,0.20,j)

np.savetxt('abxonlyresults20.csv',abxresults,delimiter=',',header=header,comments='')

for l in range(n_runs):
	print "ABX Only (0.40) Iteration %i of %i" % (l+1,n_runs)
	ABXRun(CDI,0.40,j)

np.savetxt('abxonlyresults40.csv',abxresults,delimiter=',',header=header,comments='')

for l in range(n_runs):
	print "ABX Only (0.60) Iteration %i of %i" % (l+1,n_runs)
	ABXRun(CDI,0.60,j)

np.savetxt('abxonlyresults60.csv',abxresults,delimiter=',',header=header,comments='')

for l in range(n_runs):
	print "ABX Only (0.80) Iteration %i of %i" % (l+1,n_runs)
	ABXRun(CDI,0.80,j)

np.savetxt('abxonlyresults80.csv',abxresults,delimiter=',',header=header,comments='')

for l in range(n_runs):
	print "ABX Only (1.00) Iteration %i of %i" % (l+1,n_runs)
	ABXRun(CDI,1.00,j)

np.savetxt('abxonlyresults100.csv',abxresults,delimiter=',',header=header,comments='')
del abxresults

#########################################
# FMT Intervention - ABX & PPI Patients #
#########################################
abxppiresults = np.zeros([n_runs, 4])

def ABXPPIRun(model,percent,iteration):
	model.DoStochSim()
	model.ChangeParameter('chi',percent)
	model.ChangeParameter('phi_c',0.001689646)	
	data = model.data_stochsim.getSpecies()
	Treated = data[-1,8]
	Incident = data[-1,10]
	Recur = data[-1,14]
	abxppiresults[iteration,0] = Treated
	abxppiresults[iteration,1] = Incident
	abxppiresults[iteration,2] = Recur
	abxppiresults[iteration,3] = percent
	
for m in range(n_runs):
	print "ABX Only (0.20) Iteration %i of %i" % (l+1,n_runs)
	ABXPPIRun(CDI,0.20,j)

np.savetxt('abxppiresults20.csv',abxppiresults,delimiter=',',header=header,comments='')

for m in range(n_runs):
	print "ABX Only (0.40) Iteration %i of %i" % (l+1,n_runs)
	ABXPPIRun(CDI,0.40,j)

np.savetxt('abxppiresults40.csv',abxppiresults,delimiter=',',header=header,comments='')

for m in range(n_runs):
	print "ABX Only (0.60) Iteration %i of %i" % (l+1,n_runs)
	ABXPPIRun(CDI,0.60,j)

np.savetxt('abxppiresults60.csv',abxppiresults,delimiter=',',header=header,comments='')

for m in range(n_runs):
	print "ABX Only (0.80) Iteration %i of %i" % (l+1,n_runs)
	ABXPPIRun(CDI,0.80,j)

np.savetxt('abxppiresults80.csv',abxppiresults,delimiter=',',header=header,comments='')

for m in range(n_runs):
	print "ABX Only (1.00) Iteration %i of %i" % (l+1,n_runs)
	ABXPPIRun(CDI,1.00,j)

np.savetxt('abxppiresults100.csv',abxppiresults,delimiter=',',header=header,comments='')
del abxppiresults

###############################
# FMT Intervention - Combined #
###############################
combinedresults = np.zeros([n_runs, 4])

def CombinedRun(model,percent,iteration):
	model.DoStochSim()
	model.ChangeParameter('chi',percent)
	model.ChangeParameter('chi_postCDI',percent)
	model.ChangeParameter('phi_c',0.001689646)	
	data = model.data_stochsim.getSpecies()
	Treated = data[-1,8]
	Incident = data[-1,10]
	Recur = data[-1,14]
	combinedresults[iteration,0] = Treated
	combinedresults[iteration,1] = Incident
	combinedresults[iteration,2] = Recur
	combinedresults[iteration,3] = percent
	
for n in range(n_runs):
	print "Combined (0.20) Iteration %i of %i" % (l+1,n_runs)
	CombinedRun(CDI,0.20,j)

np.savetxt('combinedresults20.csv',combinedresults,delimiter=',',header=header,comments='')

for n in range(n_runs):
	print "Combined (0.40) Iteration %i of %i" % (l+1,n_runs)
	CombinedRun(CDI,0.40,j)

np.savetxt('combinedresults40.csv',combinedresults,delimiter=',',header=header,comments='')

for n in range(n_runs):
	print "Combined (0.60) Iteration %i of %i" % (l+1,n_runs)
	CombinedRun(CDI,0.60,j)

np.savetxt('combinedresults60.csv',combinedresults,delimiter=',',header=header,comments='')

for n in range(n_runs):
	print "Combined (0.80) Iteration %i of %i" % (l+1,n_runs)
	CombinedRun(CDI,0.80,j)

np.savetxt('combinedresults80.csv',combinedresults,delimiter=',',header=header,comments='')

for n in range(n_runs):
	print "Combined (1.00) Iteration %i of %i" % (l+1,n_runs)
	CombinedRun(CDI,1.00,j)

np.savetxt('combinedresults100.csv',combinedresults,delimiter=',',header=header,comments='')




