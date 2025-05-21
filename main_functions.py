import numpy as np 
import matplotlib.pyplot as plt
import pandas as pandas
import sys
import os
# Add the directory containing the project to sys.path
current_path = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_path)


# Import Curves Function 
from Curves_Functions import Curves
# Import VBA.py File that containes all the VBA CST codes
import Functions as VBA
# Import Components.py file needed to create the MZMs and some other structures
import Components as Components




# ##################################################################################
# # Call  CST 
# ##################################################################################

# # Call CST as Object and Open existing Project 

# # Local path to CST project file 
# cst_path = r'C:/Users/..../CST_Project/' # path


# # Open an existing CST Tamplate Project
# cst_project2 = 'MZM_Template' # CST project
# cst_project_path = cst_path + cst_project2 + '.cst'


# # Call CST from CST Lib. You can find the CST Lib in ".../CST Studio Suite 2025/AMD64/python_cst_libraries"
# sys.path.append(r"C:/Program Files (x86)/CST Studio Suite 2025/AMD64/python_cst_libraries")
# import cst
# print(cst.__file__)
# # Open CST as software
# from cst.interface import DesignEnvironment



# mycst = DesignEnvironment.new()
# mycst.quiet_mode_enabled()
# mycst.in_quiet_mode()


# # Create new Empty MWS CST Project
# proj = mycst.new_mws()


# # Open Existing Project. 
# proj = mycst.open_project(cst_project_path)


# # save project into Patbh
# Name = "Save_CST"
# proj.save(os.path.join(cst_path, Name), allow_overwrite=  True)


# #Close CST Project
# proj.close()


# # Close CST
# mycst.close()



# ##################################################################################
# # Set Units in CST
# ################################################################################## 


# # Set Units
# UnitParams = {}
# UnitParams['Dimensions'] = "um"
# UnitParams['Frequency'] = "THz"
# UnitParams['Time'] = "ns"
# UnitParams['Temperature'] = "degC"
# # UnitParams['Voltage'] = "V"
# # UnitParams['Current'] = "A"
# # UnitParams['Resistance'] = "Ohm"
# # UnitParams['Conductance'] = "Siemens"
# # UnitParams['Capacitance'] = "PikoF"
# # UnitParams['Inductance'] = "NanoH"

# Units = VBA.SetUnits(UnitParams)
# proj.schematic.execute_vba_code(Units, timeout=None)


# ##################################################################################
# # Set Frequency or Wavelength in CST
# ################################################################################## 

# # Set Wavelength of operation
# Parameters = {}
# Parameters["Min Wavelength"] = 1.5
# Parameters["Max Wavelength"] = 1.7


# data = VBA.SetSimWavelength(Parameters)
# proj.schematic.execute_vba_code(data, timeout=None)


# # Set Freqeuency of operation
# Parameters = {}
# Parameters["Min Frequency"] = 1
# Parameters["Max Frequency"] = 150

# data = VBA.SetSimFreqeuncy(Parameters)
# proj.schematic.execute_vba_code(data, timeout=None)




# ##################################################################################
# # Add and Remove Global Variables to CST Simulation
# ################################################################################## 



# # Add Global Parameters
# GlobalParam = {}
# GlobalParam['Length'] = 80
# GlobalParam['Width_GND'] = 10
# GlobalParam['Width_Signal'] = 5
# GlobalParam['Width_Gap'] = 5
# GlobalParam['Hight'] = 2
# TotalWidth = 2*GlobalParam['Width_GND'] + GlobalParam['Width_Signal'] + 2*GlobalParam['Width_Gap']

# Globals  = VBA.AddGlobalParameter(GlobalParam)
# proj.schematic.execute_vba_code(Globals, timeout=None)


# # Delete Global Parameters
# GlobalParamDel = {}
# GlobalParamDel['Length'] = 10

# delete = VBA.DeleteGlobalParameter(GlobalParamDel)
# proj.schematic.execute_vba_code(delete, timeout=None)



# ##################################################################################
# # Define Background 
# ################################################################################## 


# # Define Background
# Params = {}
# Params["Type Background"] = "Normal"
# Params["Xmin Background"] = 5
# Params["Xmax Background"] = 5
# Params["Ymin Background"] = 5
# Params["Ymax Background"] = 5
# Params["Zmin Background"] = 5
# Params["Zmax Background"] = 5

# BackObj = VBA.BackgroundSet(Params)
# proj.schematic.execute_vba_code(BackObj, timeout=None)


# ##################################################################################
# # Define Boundary 
# ################################################################################## 



# # Set Boundary
# BoundaryParams = {}
# BoundaryParams["Xmin Boundary"] = 'open'
# BoundaryParams["Xmax Boundary"] = 'open'
# BoundaryParams["Ymin Boundary"] = 'open'
# BoundaryParams["Ymax Boundary"] = 'open'
# BoundaryParams["Zmin Boundary"] = 'open'
# BoundaryParams["Zmax Boundary"] = 'open'
# BoundaryParams["Xsymmetry Boundary"] = 'none'
# BoundaryParams["Ysymmetry Boundary"] = 'none'
# BoundaryParams["Zsymmetry Boundary"] = 'none'


# BoundarySet = VBA.SetBoundary(BoundaryParams)
# proj.schematic.execute_vba_code(BoundarySet, timeout=None)



# ##################################################################################
# # Pick face/ Clear Pick  
# ################################################################################## 


# # Pick Face or Centerpoint
# PicParams = {}
# PicParams["Option"] = "Face"
# # PicParams["Option"] = "Centerpoint"
# PicParams["Object"] = "Electrode"
# PicParams["Face Number"] = 1
# PickFace = VBA.Pick(PicParams)
# proj.schematic.execute_vba_code(PickFace, timeout=None)


# # Clear Picks
# ClearPick = VBA.ClearAllPicks()
# proj.schematic.execute_vba_code(ClearPick, timeout=None)





# ##################################################################################
# # Create and Plot Bezier, Cosinus and Euler Curves
# ################################################################################## 

# # Define Curves Parameters and Data
# Lenght = 100
# Offset = 40
# points = 100



# # Generate the Bezier and Cos points
# ObjCurves = Curves(Lenght, Offset, points)
# BezierCuve = ObjCurves.Bezier_Curve()
# CosinusCurve = ObjCurves.Cosinus_Curve()
# EulerCurve = ObjCurves.Euler_Curve()


# # Plot the Curves
# plt.figure()
# plt.plot(BezierCuve[:,0], BezierCuve[:,1], color = "red", label = "Bezier Curve")
# plt.plot(CosinusCurve[:,0], CosinusCurve[:,1], color = "blue", label = "Cosinus Curve")
# plt.plot(EulerCurve[:,0], EulerCurve[:,1], color = "green", label = "Euler Curve")
# plt.xlabel("Length S-Bend/$\mu m$")
# plt.ylabel("Offset S-Bends/ $\mu m$")
# plt.legend(loc = "best")
# plt.grid()
# plt.show()



# ##################################################################################
# # Create Bondwire and Curves with given trajectories 
# ################################################################################## 


# Parameters = {}
# Parameters['X1'] = 0
# Parameters['Y1'] = 0
# Parameters['Z1'] = 0
# Parameters['X2'] = 5
# Parameters['Y2'] = 5
# Parameters['Z2'] = 0

# Points = {}
# x = []
# y = []
# for i in range(0, 100):
#     x.append(i)
#     y.append(i*4)
# x = np.array(x)
# y = np.array(y)
# Points['X'] = x
# Points['Y'] = y

# Points['X'] = CosinusCurve[:,0]
# Points['Y'] = CosinusCurve[:,1]




# BondWire = Components.BondWire(NameWire = "TestWire" ,Coordinates = Parameters, Height = 1, Radius = 0.5 , BondwireType = "Spline", Material = "Copper (annealed)",  NameFolder = "BondWire")
# proj.schematic.execute_vba_code(BondWire, timeout=None)
# ToSolid1 = VBA.ToSolid(SolidName = "TestSolid", CurveName = "TestWire", NameFolder = "BondWire", Material = "Copper (annealed)")
# proj.schematic.execute_vba_code(ToSolid1, timeout=None)


# Curve = VBA.Curve(CurveName = "TestCurve", Points = Points)
# proj.schematic.execute_vba_code(Curve, timeout=None)
# DataCurveWire = Components.CurveWire(NameWire = "CurveWire2", Radius = 0.5 , Points = Points, Material = "PEC", CurveFolderName = "TestCurve", CurveName = "TestCurve")
# proj.schematic.execute_vba_code(DataCurveWire, timeout=None)
# ToSolid2 = VBA.ToSolid(SolidName = "TestSolidCurve2", CurveName = "CurveWire2", NameFolder = "CurveWire2", Material = "Copper (annealed)")
# proj.schematic.execute_vba_code(ToSolid2, timeout=None)



# ##################################################################################
# # Create an optical waveguide
# ################################################################################## 



# #Create Waveguide with Ports 
# Parameters = {}
# Parameters["Lenght WG"] = 5
# Parameters["Hight WG"] = 0.3
# Parameters["Width WG"] = 1
# Parameters["Substrate Height"] = 1
# Parameters["Slab Heigh"] = 0

# WG = Components.Squere_Waveguide(Parameters, proj)




# ##################################################################################
# # Create an LNOI MZM
# ################################################################################## 



# # MZM Parameters
# Parameters = {}
# Parameters["Electrodes Lenght"] = 50
# Parameters["Width GND"] = 40    
# Parameters["Width Signal"] = 10 
# Parameters["Width WG"] = 0.8
# Parameters["Gap"] = 1.565
# Parameters["angle"] = 35
# Parameters["High Electrodes"] = 0.8
# Parameters["High WG"] = 0.4
# Parameters["High Slab"] = 0.2 
# Parameters["High Substrate"] = 2
# Parameters["Angle X"] = 0
# Parameters["Angle Y"] = 90
# Parameters["Angle Z"] = 0 


# # Ganze Mach-Zehnder Modulator
# Components.MZM(Parameters, proj)





# # Ports on Electrodes for MZM Modulator
# Parameters["Electrodes_Names"] = ["Electrode_Left:Electrode_Left1", "Electrode_Right:Electrode_Right1" , "Signal:Signal1"]
# Parameters["Orientation"] = ["Positive", "Positive"]
# Parameters["Coordinates"] = "Picks"
# Parameters["Span"] = [[[3,3],[3,3]], [[3,3],[3,3]]]
# Parameters["Potential"] = [1,2]
# Parameters["Port Number"] = [1,2]
# Parameters["Polarity"] = ["Positive", "Negative"]
# Parameters["Face ID"] = [4,6]

# VBA.WaveguidePorts_on_Electrodes_MZM(Parameters, proj)




# # Optical Waveguides Port for MZM Modulator
# Parameters["Electrodes_Names"] = ["Waveguide_Left:Waveguide_Left", "Waveguide_Right:Waveguide_Right"]
# Parameters["Orientation"] = ["Positive", "Positive"]
# Parameters["Coordinates"] = "Picks"
# Parameters["Span"] = [[[3,3],[3,3]], [[3,3],[3,3]]]
# Parameters["Potential"] = [1,2]
# Parameters["Port Number"] = [1,2]
# Parameters["Polarity"] = ["Positive", "Negative"]
# Parameters["Face ID"] = [5,6]

# VBA.Optical_WaveguidePorts_MZM(Parameters, proj)


# ##################################################################################
# # Create an LNOI Phase Modulator
# ################################################################################## 


# # Phase modulator parameters
# Parameters = {}
# Parameters["Electrodes Lenght"] = 50
# Parameters["Width GND"] = 40    
# Parameters["Width Signal"] = 10 #50
# Parameters["Width WG"] = 0.8
# Parameters["Gap"] = 1.565
# Parameters["angle"] = 35
# Parameters["High Electrodes"] = 0.8
# Parameters["High WG"] = 0.4
# Parameters["High Slab"] = 0.2 
# Parameters["High Substrate"] = 2
# Parameters["Angle X"] = 0
# Parameters["Angle Y"] = 90
# Parameters["Angle Z"] = 0 



# # Simple EO Phase Modulator
# Components.PhaseModulator(Parameters, proj)






# # Ports on Optical Waveguides for Phase Modulator
# Parameters["Electrodes_Names"] = ["Waveguide_Left:Waveguide_Left"]
# Parameters["Orientation"] = ["Positive", "Positive"]
# Parameters["Coordinates"] = "Picks"
# Parameters["Span"] = [[[3,3],[3,3]], [[3,3],[3,3]]]
# Parameters["Potential"] = [1,2]
# Parameters["Port Number"] = [1,2]
# Parameters["Polarity"] = ["Positive", "Negative"]
# Parameters["Face ID"] = [5,6]


# VBA.Optical_WaveguidePorts_MZM(Parameters, proj)



# # Ports on Electrodes for Phase Modulator
# Parameters["Electrodes_Names"] = ["Electrode:Electrode1", "Signal:Signal1"]
# Parameters["Orientation"] = ["Positive", "Positive"]
# Parameters["Coordinates"] = "Picks"
# Parameters["Span"] = [[[3,3],[3,3]], [[3,3],[3,3]]]
# Parameters["Potential"] = [1,2]
# Parameters["Port Number"] = [1,2]
# Parameters["Polarity"] = ["Positive", "Negative"]
# Parameters["Face ID"] = [4,6]

# VBA.Optical_WaveguidePorts_MZM(Parameters, proj)




# # Set Discrete Ports Phase Modulator
# Parameters["Electrodes_Names"] = ["Signal:Signal1", "Electrode:Electrode1"]
# Parameters["Face ID"] = [4,6]
# Parameters["Discrete Port Type"] = "Voltage"
# Parameters["Port Impedance"] = 50
# Parameters["Port Voltage"] = 2
# Parameters["Port Current"] = 1
# Parameters["Port Radius"] = 0
# Parameters["Port Number"] = [3,4]

# VBA.Discrete_Port(Parameters, proj)





# ##################################################################################
# # Set Time Solver in CST
# ################################################################################## 


# # Set Time Solver
# Parameters= {}
# Parameters["Accuracy"] = 30
# Parameters["Caclculate Modes Only"] = False
# Parameters["Auto Impedance"] = True
# Parameters["Impedance"] = 50
# Parameters["Source Port"]  = 1
# Parameters["Solver Mesh Type"] = "TLM"


# Solver = VBA.SetTimeSolver(Parameters)
# proj.schematic.execute_vba_code(Solver, timeout=None)


# ##################################################################################
# # Set Frequency Solver in CST
# ################################################################################## 


# # Set Freq Solver
# Parameters= {}
# Parameters["Accuracy"] = 30
# Parameters["Caclculate Modes Only"] = False
# Parameters["Auto Impedance"] = True
# Parameters["Impedance"] = 50
# Parameters["Source Port"]  = 1

# Solver = VBA.SetFreqSolver(Parameters)
# proj.schematic.execute_vba_code(Solver, timeout=None)



# ##################################################################################
# # Set Monitors for different wavelengths
# ################################################################################## 

# # Set Monitors
# Parameters = {}
# Parameters["Wavelength"] = 1.55e-6
# Parameters["Monitor Type"] = "Efield"
# Monitor = VBA.CreateEfieldMonitor(Parameters)
# proj.schematic.execute_vba_code(Monitor, timeout=None)


# Parameters = {}
# Parameters["Wavelength"] = 1.55e-6
# Parameters["Monitor Type"] = "Powerflow"
# Monitor = VBA.CreateEfieldMonitor(Parameters)
# proj.schematic.execute_vba_code(Monitor, timeout=None)


# ##################################################################################
# # Set Mesh
# ################################################################################## 

# # Set Mesh 
# Parameters = {}
# Parameters["Mesh Type"] = "PBA"
# Parameters["Mesh Cells Near Object"] = 30
# Parameters["Mesh Cells far Object"] = 30


# Mesh = VBA.SetMesh(Parameters)
# proj.schematic.execute_vba_code(Mesh, timeout=None)

    
# ##################################################################################
# # Change Solver Type
# ################################################################################## 

# # Change Solver type, for now only Time solver can be changed 
# Solver = VBA.ChangeSolverType("Time")
# proj.schematic.execute_vba_code(Solver, timeout=None)




# ##################################################################################
# # Run Time Solver
# ################################################################################## 

# # Start Time Solver
# Start = VBA.StartTimeSolver()
# proj.schematic.execute_vba_code(Start, timeout=None)


# #################################################################################
# # Extract Results after Simulation to CSV Files
# ################################################################################## 



# # Extract Results 
# # ExportResults(Path to save file into, Name of extracted file, Exact name of the CST Result path)
# Data = VBA.ExportResults("C:\\Users\\Martin\\Desktop\\CST_Project", "S21", "1D Results\\S-Parameters\\S2,1")
# proj.schematic.execute_vba_code(Data, timeout=None)



