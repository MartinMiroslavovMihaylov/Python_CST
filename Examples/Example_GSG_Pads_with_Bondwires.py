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
# Create the GSG PADS with Bondwires
# =============================================================================

# Set Units
Parameters = {}
Parameters["Unit Length"] = "um"
Parameters["Unit Frequency"] = "GHz"
Parameters["Unit Time"] = "ns"
Parameters["Unit Temperature"] = "K"
obj.set_Units(Parameters)


# Set Boundary and background
Params = {}
Params["Type Background"] = "Normal"
Params["Xmin Background"] = 0
Params["Xmax Background"] = 0
Params["Ymin Background"] = 0
Params["Ymax Background"] = 0
Params["Zmin Background"] = 0
Params["Zmax Background"] = 0

obj.setBackground(Params)


BoundaryParams = {}
BoundaryParams["Xmin Boundary"] = 'open'
BoundaryParams["Xmax Boundary"] = 'open'
BoundaryParams["Ymin Boundary"] = 'open'
BoundaryParams["Ymax Boundary"] = 'open'
BoundaryParams["Zmin Boundary"] = 'open'
BoundaryParams["Zmax Boundary"] = 'open'
BoundaryParams["Xsymmetry Boundary"] = 'none'
BoundaryParams["Ysymmetry Boundary"] = 'none'
BoundaryParams["Zsymmetry Boundary"] = 'none'

obj.setBoundary(BoundaryParams)


# Set Frequency Range for the simulation
Parameters = {}
Parameters["Min Frequency"] = 0     
Parameters["Max Frequency"] = 200
obj.setSimFreqeuncy(Parameters)


# Set Time Solver with TLM Mesh 
Parameters= {}
Parameters["Accuracy"] = 40
Parameters["Caclculate Modes Only"] = False
Parameters["Auto Impedance"] = True
Parameters["Impedance"] = 50
Parameters["Source Port"]  = 1
Parameters["Solver Mesh Type"] = "TLM"
obj.setTimeSolver(Parameters)

# Set parameters 
Parameters = {}
Parameters["PAD Width GND"] = 100
Parameters["PAD Width Signal"] = 40
Parameters["PAD Length"] = 100
Parameters["PAD Thickness"] = 2.8


Parameters["PADs Distance"] = 219
Parameters["Bonwire Height"] = 20
Parameters["Bonwire Radius"] = 25/2
Parameters["Glue Thickness"] = 50
Parameters["Floating Shield Thickness"]= 35

Parameters["Port Y Span"] = 70
Parameters["Port Z Span"] = 200
Parameters["Accuracy"] = 40



# Parameters["Probes"] is set to False. If True two Probes without Ports will be set!
Parameters["Probes"] = False
obj.GSG_Bondwire_ChipToChip_connection(Parameters)

        