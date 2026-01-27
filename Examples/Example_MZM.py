import numpy as np 
import matplotlib.pyplot as plt
import pandas as pandas
import sys
import os
# Add the directory containing the project to sys.path
current_path = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_path)
from CST_Constructor import CST_Commands, Curves


# =============================================================================
# Open new MWS CST 
# =============================================================================


obj = CST_Commands()
obj.New_Project("MWS")


# =============================================================================
# Create the MZM
# =============================================================================


# Set Optical or Electrical Temptlates for the project 
# Units Properties
Parameters = {}
Parameters['Dimensions'] = "um"
Parameters['Frequency']  = "THz"
# Parameters['Frequency']  = "GHz"
Parameters['Time'] = "ns"
Parameters['Temperature'] = "degC"

# Set FreqeueWavelength Range 
Parameters["Min Wavelength"] = 1.5
Parameters["Max Wavelength"] = 1.6
# Parameters["Min Frequency"] = 1
# Parameters["Max Frequency"] = 150

# Set Background
Parameters["Type Background"] = "Normal"
Parameters["Xmin Background"] = 0
Parameters["Xmax Background"] = 0
Parameters["Ymin Background"] = 0
Parameters["Ymax Background"] = 0
Parameters["Zmin Background"] = 0
Parameters["Zmax Background"] = 0


# Set Boundary
Parameters["Xmin Boundary"] = 'open'
Parameters["Xmax Boundary"] = 'open'
Parameters["Ymin Boundary"] = 'open'
Parameters["Ymax Boundary"] = 'open'
Parameters["Zmin Boundary"] = 'open'
Parameters["Zmax Boundary"] = 'open'
Parameters["Xsymmetry Boundary"] = 'none'
Parameters["Ysymmetry Boundary"] = 'none'
Parameters["Zsymmetry Boundary"] = 'none'

# Mesh Settings
Parameters["Mesh Type"] = "PBA"
Parameters["Mesh Cells Near Object"] = 30
Parameters["Mesh Cells far Object"] = 30

obj.setOpticalSimulationProperties(Parameters)



# MZM Parameters
Parameters = {}

Parameters['Substrate Height'] = 1e-6
Parameters["Optical"] = ["Au (Gold) - CRC", "SiO2 (Glass) - Palik"]
Parameters["Electrical"] = ["Air", "Au (Gold) - CRC", "LiNbO3 semiconductor - X/Y cut (Lithium Niobate)", "SiO2 (Glass) - Sze"]
Parameters['angle'] = 30
Parameters['Slab Height'] = 0.3e-6
Parameters['WG Height'] = 0.3e-6
Parameters['WG Width'] = 1e-6
Parameters['WG Length'] = 10e-6
Parameters["GND Electrodes Width"] = 5e-6
Parameters["Signal Electrodes Width"] = 3e-6
Parameters["Electrodes Height"] = 0.8e-6
Parameters["Gap"] = 1.4e-6
Parameters["Wavelength"] = 1.55e-6


# Create MZM
obj.MZM(Parameters) 
