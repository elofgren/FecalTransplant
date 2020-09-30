#!/usr/bin/env python

###################################################
# Dynamic Transmission Model of C. difficile      #
# Fecal Microbiota Transplant Based Interventions #
# Gillespie Direct Method Implementation          #
# Author: Eric Lofgren (Eric.Lofgren@gmail.com)   #
###################################################

# Module Imports
import os
import stochpy
import pylab as pl
import numpy as np

# Set up simulation parameters
start_time = 0.0
end_time = 8760
n_runs = 1000

# Headers for output file
header = "Treated, Incident, Recur, Level"

#######################################
# Baseline Scenario - No Intervention #
#######################################
CDI = stochpy.SSA()
CDI.Model(File='cdiff_FT.psc', dir=os.getcwd())

baselineresults = np.zeros([n_runs, 4])

def BaselineRun(model,iteration):
	model.Endtime(end_time)
	model.DoStochSim()
	data = model.data_stochsim.getSpecies()
	Treated = data[-1,8]
	Incident = data[-1,10]
	Recur = data[-1,14]
	baselineresults[iteration,0] = Treated
	baselineresults[iteration,1] = Incident
	baselineresults[iteration,2] = Recur
	baselineresults[iteration,3] = 0

for i in range(n_runs):
	print "Baseline Iteration %i of %i" % (i+1,n_runs)
	BaselineRun(CDI,i)
i = 0	
np.savetxt('baseline.csv',baselineresults,delimiter=',',header=header,comments='')
del baselineresults

########################################
# FMT Intervention - CDI Patients Only #
########################################
CDI = stochpy.SSA()
CDI.Model(File='cdiff_FT.psc', dir=os.getcwd())

cdiresults = np.zeros([n_runs, 4])

def CDIRun(model,percent,iteration):
	model.Endtime(end_time)
	model.ChangeParameter('chi_postCDI',percent)	
	model.DoStochSim()
	data = model.data_stochsim.getSpecies()
	Treated = data[-1,8]
	Incident = data[-1,10]
	Recur = data[-1,14]
	cdiresults[iteration,0] = Treated
	cdiresults[iteration,1] = Incident
	cdiresults[iteration,2] = Recur
	cdiresults[iteration,3] = percent*100
	
for i in range(n_runs):
	print "CDI Only (0.20) Iteration %i of %i" % (i+1,n_runs)
	CDIRun(CDI,0.20,i)
i = 0
np.savetxt('cdionlyresults20.csv',cdiresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "CDI Only (0.40) Iteration %i of %i" % (i+1,n_runs)
	CDIRun(CDI,0.40,i)
i = 0
np.savetxt('cdionlyresults40.csv',cdiresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "CDI Only (0.60) Iteration %i of %i" % (i+1,n_runs)
	CDIRun(CDI,0.60,i)
i = 0
np.savetxt('cdionlyresults60.csv',cdiresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "CDI Only (0.80) Iteration %i of %i" % (i+1,n_runs)
	CDIRun(CDI,0.80,i)
i = 0
np.savetxt('cdionlyresults80.csv',cdiresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "CDI Only (1.00) Iteration %i of %i" % (i+1,n_runs)
	CDIRun(CDI,1.00,i)
i = 0
np.savetxt('cdionlyresults100.csv',cdiresults,delimiter=',',header=header,comments='')
del cdiresults

########################################
# FMT Intervention - ABX Patients Only #
########################################
CDI = stochpy.SSA()
CDI.Model(File='cdiff_FT.psc', dir=os.getcwd())

abxresults = np.zeros([n_runs, 4])

def ABXRun(model,percent,iteration):
	model.Endtime(end_time)
	model.ChangeParameter('chi',percent)
	model.ChangeParameter('phi_c',0.001102321)
	model.DoStochSim()	
	data = model.data_stochsim.getSpecies()
	Treated = data[-1,8]
	Incident = data[-1,10]
	Recur = data[-1,14]
	abxresults[iteration,0] = Treated
	abxresults[iteration,1] = Incident
	abxresults[iteration,2] = Recur
	abxresults[iteration,3] = percent*100
	
for i in range(n_runs):
	print "ABX Only (0.20) Iteration %i of %i" % (i+1,n_runs)
	ABXRun(CDI,0.20,i)
i = 0
np.savetxt('abxonlyresults20.csv',abxresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "ABX Only (0.40) Iteration %i of %i" % (i+1,n_runs)
	ABXRun(CDI,0.40,i)
i = 0
np.savetxt('abxonlyresults40.csv',abxresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "ABX Only (0.60) Iteration %i of %i" % (i+1,n_runs)
	ABXRun(CDI,0.60,i)
i = 0
np.savetxt('abxonlyresults60.csv',abxresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "ABX Only (0.80) Iteration %i of %i" % (i+1,n_runs)
	ABXRun(CDI,0.80,i)
i = 0
np.savetxt('abxonlyresults80.csv',abxresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "ABX Only (1.00) Iteration %i of %i" % (i+1,n_runs)
	ABXRun(CDI,1.00,i)
i = 0
np.savetxt('abxonlyresults100.csv',abxresults,delimiter=',',header=header,comments='')
del abxresults

#########################################
# FMT Intervention - ABX & PPI Patients #
#########################################
CDI = stochpy.SSA()
CDI.Model(File='cdiff_FT.psc', dir=os.getcwd())

abxppiresults = np.zeros([n_runs, 4])

def ABXPPIRun(model,percent,iteration):
	model.Endtime(end_time)
	model.ChangeParameter('chi',percent)
	model.ChangeParameter('phi_c',0.001689646)
	model.DoStochSim()
	data = model.data_stochsim.getSpecies()
	Treated = data[-1,8]
	Incident = data[-1,10]
	Recur = data[-1,14]
	abxppiresults[iteration,0] = Treated
	abxppiresults[iteration,1] = Incident
	abxppiresults[iteration,2] = Recur
	abxppiresults[iteration,3] = percent*100
	
for i in range(n_runs):
	print "ABX Only (0.20) Iteration %i of %i" % (i+1,n_runs)
	ABXPPIRun(CDI,0.20,i)
i = 0
np.savetxt('abxppiresults20.csv',abxppiresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "ABX Only (0.40) Iteration %i of %i" % (i+1,n_runs)
	ABXPPIRun(CDI,0.40,i)
i = 0
np.savetxt('abxppiresults40.csv',abxppiresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "ABX Only (0.60) Iteration %i of %i" % (i+1,n_runs)
	ABXPPIRun(CDI,0.60,i)
i = 0
np.savetxt('abxppiresults60.csv',abxppiresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "ABX Only (0.80) Iteration %i of %i" % (i+1,n_runs)
	ABXPPIRun(CDI,0.80,i)
i = 0
np.savetxt('abxppiresults80.csv',abxppiresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "ABX Only (1.00) Iteration %i of %i" % (i+1,n_runs)
	ABXPPIRun(CDI,1.00,i)
i = 0
np.savetxt('abxppiresults100.csv',abxppiresults,delimiter=',',header=header,comments='')
del abxppiresults

###############################
# FMT Intervention - Combined #
###############################
CDI = stochpy.SSA()
CDI.Model(File='cdiff_FT.psc', dir=os.getcwd())

combinedresults = np.zeros([n_runs, 4])

def CombinedRun(model,percent,iteration):
	model.Endtime(end_time)
	model.ChangeParameter('chi',percent)
	model.ChangeParameter('chi_postCDI',percent)
	model.ChangeParameter('phi_c',0.001689646)	
	model.DoStochSim()
	data = model.data_stochsim.getSpecies()
	Treated = data[-1,8]
	Incident = data[-1,10]
	Recur = data[-1,14]
	combinedresults[iteration,0] = Treated
	combinedresults[iteration,1] = Incident
	combinedresults[iteration,2] = Recur
	combinedresults[iteration,3] = percent*100
	
for i in range(n_runs):
	print "Combined (0.20) Iteration %i of %i" % (i+1,n_runs)
	CombinedRun(CDI,0.20,i)
i = 0
np.savetxt('combinedresults20.csv',combinedresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "Combined (0.40) Iteration %i of %i" % (i+1,n_runs)
	CombinedRun(CDI,0.40,i)
i = 0
np.savetxt('combinedresults40.csv',combinedresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "Combined (0.60) Iteration %i of %i" % (i+1,n_runs)
	CombinedRun(CDI,0.60,i)
i = 0
np.savetxt('combinedresults60.csv',combinedresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "Combined (0.80) Iteration %i of %i" % (i+1,n_runs)
	CombinedRun(CDI,0.80,i)
i = 0
np.savetxt('combinedresults80.csv',combinedresults,delimiter=',',header=header,comments='')

for i in range(n_runs):
	print "Combined (1.00) Iteration %i of %i" % (i+1,n_runs)
	CombinedRun(CDI,1.00,i)
i = 0
np.savetxt('combinedresults100.csv',combinedresults,delimiter=',',header=header,comments='')




