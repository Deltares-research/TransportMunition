# Input file for transport model

# Packages required
import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
import warnings
from sklearn.decomposition import PCA
from scipy.optimize import curve_fit, fsolve
from scipy.signal import find_peaks
from scipy.interpolate import interp1d, NearestNDInterpolator
from scipy.spatial import cKDTree

# Import functions
sys.path.append('/DataBodem')
sys.path.append('/Dwarsdoorsnede')
sys.path.append('/Excel')
sys.path.append('/Functies')
from Bathy import bathy
from Rivier import rivier
from Transport import transport
warnings.filterwarnings("ignore")


# Input bathy
pp = 0   # Yes or no, generation new roughness values/submergence
filesBT = {
   	 "OOlijst": "OOlijst.xlsx", 	     		# Lijst OO bommen
     "testcases": "Testcases.xlsx",        	# Database fit
     "output_file": "Ruwheid.xlsx",				# Ruwheid per bom
}
if pp == 1:
    bathy(files=filesBT)

# Input parameter transport
parametersTM = {
    "rho_w": 1000,  # Water density (kg/m^3)
    "g": 9.81,      # Gravitational acceleration (m/s^2)
    "mu": 0.001     # Dynamic viscosity (Pa·s or kg/(m·s))
}

filesTM = {
    "database": "Database.xlsx",        	# Database fit
    "OOlijst": "OOlijst.xlsx", 	     		# Lijst OO bommen
    "ruwheid": "Ruwheid.xlsx",				# Ruwheid per bom
    "testcases": "Testcases.xlsx",        	# Database fit
    "output_file": "minimumSnelheid.xlsx",  # Output per OO minimum snelheid
}
transport(params=parametersTM, files=filesTM)

# Input parameter rivier
# filesRV = {
#     "testcases": "Testcases.xlsx",        	# Database fit
#     "OOlijst": "OOlijst.xlsx", 	     		# Lijst OO bommen
#     "snelheid": "minimumSnelheid.xlsx",     # Lijst minimumsnelheden
#     "ruwheid": "Ruwheid.xlsx",				# Ruwheid per bom
#     "output_file": "beweging.txt",         	# Snelheid water, beweging object
# }
# rivier(files=filesRV)
