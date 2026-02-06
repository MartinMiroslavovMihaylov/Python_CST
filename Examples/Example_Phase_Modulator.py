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
# Create the Phase modulator
# =============================================================================


# Set Optical or Electrical Temptlates for the project 
# Units Properties
Parameters = {}
Parameters['Length'] = "um"
Parameters['Frequency']  = "THz"
Parameters['Time'] = "ns"
Parameters['Temperature'] = "degC"

# Set FreqeueWavelength Range 
Parameters["Min Wavelength"] = 1.5
Parameters["Max Wavelength"] = 1.6


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



# Create Phase Modulator
Parameters = {}
Parameters["Electrodes Lenght"] = 50
Parameters["Width GND"] = 40    
Parameters["Width Signal"] = 10 
Parameters["Width WG"] = 0.8
Parameters["Gap"] = 1.565
Parameters["Angle"] = 35
Parameters["Height Electrodes"] = 0.8
Parameters["Height WG"] = 0.4
Parameters["Height Slab"] = 0.2 
Parameters["Height Substrate"] = 2
Parameters["Angle X"] = 0
Parameters["Angle Y"] = 90
Parameters["Angle Z"] = 0 

obj.PhaseModulator(Parameters)