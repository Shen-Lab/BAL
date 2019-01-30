
import os
import numpy as np

class Configurations(object):
	"""
	Parameters of the TNMA program are set and stored in the Configurations class.
	
	This is the class to set specific program parameters.
	"""

	def __init__(self):
		"""
		Constructor of the Configurations class. 
		"""
		### Paths ###
		
		# Path to the output

		outputPath = os.path.dirname(__file__)
		self.outputPath = outputPath[:-5]
		# self.outputPath = "~/Github/Github_shen/cNMA/Manual/Example/Example1/"

				
		# os.environ['paramsPath'] = self.outputPath+'params.txt'
		# index = os.popen("if [ -e $paramsPath ]; then d=true; else d=false; fi; echo -n $d").read()
		# if index == 'true':
		# 	params = np.loadtxt(self.outputPath+"params.txt", dtype = str)
		# 	customHRdistance = float(params[0])
		# 	customForceConstant = float(params[1])
		# 	maxModesToCalculate = int(params[2])
		# elif index == 'false':

		# Key parameters
		customHRdistance = 12.0
		customForceConstant = 0.25
		maxModesToCalculate = 400


		# Experiment name prefix to be used to create the results output folder
		self.experimentNamePrefix = "Output"
		
		# NMAUnified investigationsOn on "Individual" or "Complex"
		self.investigationsOn = "Complex"
		
		# measures on "whole" if true, else on "interface"
		self.measuresOnWhole = True
		
		# calculate zero eigenvalue modes
		self.calculateZeroEigvalModes = True
		
		# NMAUnified align as "alpha" or "beta"
		self.align = "beta"
		
		self.complexRMSDreduction = "2k" # for complex: 1k1k or 2k or 1k1k6
		self.whichCustomHC = "HC_U1" # HC_0 or HC_U1
		self.enforceAllModesAfterProjection = True 
		
		### PDB ###
		# When parsing PDB files, filter by ProDys "protein" selection? Without this filtering (set as "None"), mismatches have occurred.
		self.filterPDB = "protein"
		
		# What atoms are subject to the matching of chains (calpha, bb or all)
		self.whatAtomsToMatch = "bb"
		
		# custom delta HR, if HC_U1, set self.customH to True, deprecated HR_A, HR_B: version A has bound structures in the second partial derivative terms, B only a penalty
		self.customH = True
		self.customHR_A = False
		self.customHR_B = False

		# self.customHRdistance = customHRdistance
		# self.customForceConstant = customForceConstant
		# self.whichCustomHIndividual = "HC_subvector"

		# Cut-off distance D for intermolecular springs
		self.customHRdistance = customHRdistance

		# Force constant gamma for intermolecular springs
		self.customForceConstant = customForceConstant

		self.whichCustomHIndividual = "HC_subvector"
		
		# Is a projection technique on the hessian (projection matrix treadment 8.27 NMa book) to be used
		self.projectHessian = True
		self.projectionStyle = "full" # "full" or "intra"
		
		# Deprecated options, Modify the HR/HL prior to projecting it
		self.modifyHDelta = False
		self.deltamultiplicatorForH = 0.5
		
		# Deprecated options, Take 1k from HR and selected 1k from HR tilde
		self.HR1kHRtilde1k = False
		self.selectKmodes = 20
		
		# Pertaining 2b, should the eigenvalues be rescaled by taking the eigenvectors and eigenvalues from the complex
		self.rescaleEigenvalues = True
		
		### RMSD Reduction ###
		# Small value to consider if the initial RMSD is > 0
		self.floatingPointThreshold = 0.000000000001
		
		# Should the RMSD reduction based on the Swarmdock betas approach be stopped upon reaching a certain number of modes
		# To not have any limit/stop, set this to a high value
		self.stopRMSDReductionAt = 400
		
		# Upper limit for mode calculation, set to very high number (1000000) to calculate 3n-6 modes
		# self.maxModesToCalculate = maxModesToCalculate
		self.maxModesToCalculate = maxModesToCalculate
		
		# Precision for RMSD beta fitting
		self.precisionBetaFitting = 1e-6
		
		# If self.guard breaks, how many extra modes to try if the RMSD reduction beta goes overdetermined AND the determinant is 0
		self.goOverdetermined = 50
		
		# How many iterations for the betas fitter
		self.maxIterBetas = 60000
		
		# guard after which beta fitting is conditioned, to disable guarded fitting, set it to the same value as self.stopRMSDReductionAt
		self.guard = 400
		
		# RMSD timeout, after how many seconds forcefully stop the beta fitting
		# It has been observed that large proteins might practically deadlock the optimizer
		self.RMSDtimeout = 120
		

