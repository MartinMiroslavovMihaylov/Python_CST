import numpy as np 
import matplotlib.pyplot as plt
import pandas as pandas
import sys
import os
# Add the directory containing the project to sys.path
current_path = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_path)
from Curves_Functions import Curves
import Functions as VBA
import Components as Components



# Define Curves Parameters and Data
Lenght = 100
Offset = 40
points = 100



# Generate the Bezier and Cos points
ObjCurves = Curves(Lenght, Offset, points)
BezierCuve = ObjCurves.Bezier_Curve()
CosinusCurve = ObjCurves.Cosinus_Curve()
EulerCurve = ObjCurves.Euler_Curve()


# plt.figure()
# plt.plot(BezierCuve[:,0], BezierCuve[:,1], color = "red", label = "Bezier Curve")
# plt.plot(CosinusCurve[:,0], CosinusCurve[:,1], color = "blue", label = "Cosinus Curve")
# plt.plot(EulerCurve[:,0], EulerCurve[:,1], color = "green", label = "Euler Curve")
# plt.xlabel("Length S-Bend/$\mu m$")
# plt.ylabel("Offset S-Bends/ $\mu m$")
# plt.legend(loc = "best")
# plt.grid()
# plt.show()



# Call CST as Object and Open existing Project 

# Local path to CST project file --> Please adapt
#C:/Users/Martin/Desktop/CST_Project/
#C:/Users/marti/Desktop/UPB Kursen/CST with Python/
#/homes/lift/martinmi/Desktop/CST_Python/
cst_path = r'C:/Users/Martin/Desktop/CST_Project/' # path
cst_project2 = 'MZM_TestSim' # CST project
# cst_project2 = 'Optical_Sim' # CST project
cst_project_path = cst_path + cst_project2 + '.cst'



import sys
# # Server
# sys.path.append(r"/opt/sct/eda/cst/2025/LinuxAMD64/python_cst_libraries")
# Windows
sys.path.append(r"C:/Program Files (x86)/CST Studio Suite 2025/AMD64/python_cst_libraries")
import cst
print(cst.__file__)
# Open CST as software
from cst.interface import DesignEnvironment



mycst = DesignEnvironment.new()
mycst.quiet_mode_enabled()
mycst.in_quiet_mode()
# mycst =cst.interface.DesignEnvironment()



# Open Existing Project
proj = mycst.open_project(cst_project_path)



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




# # Set Freqeuency of operation
# Parameters = {}
# Parameters["Min Wavelength"] = 1.5
# Parameters["Max Wavelength"] = 1.7

# data = VBA.SetSimWavelength(Parameters)
# proj.schematic.execute_vba_code(data, timeout=None)



# # Add Global Parameters
# GlobalParam = {}
# GlobalParam['Length'] = 80
# GlobalParam['Width_GND'] = 10
# GlobalParam['Width_Signal'] = 5
# GlobalParam['Width_Gap'] = 5
# GlobalParam['Hight'] = 2
# TotalWidth = 2*GlobalParam['Width_GND'] + GlobalParam['Width_Signal'] + 2*GlobalParam['Width_Gap']

# Globals  = AddGlobalParameter(GlobalParam)
# proj.schematic.execute_vba_code(Globals, timeout=None)


# # Delete Global Parameters
# GlobalParamDel = {}
# GlobalParamDel['Length'] = 10
# delete = DeleteGlobalParameter(GlobalParamDel)
# proj.schematic.execute_vba_code(delete, timeout=None)


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





# Create MZM Object
Parameters = {}
Parameters["Electrodes Lenght"] = 15
Parameters["Width GND"] = 5
Parameters["Width Signal"] = 2
Parameters["Width WG"] = 1
Parameters["Gap"] = 0.5
Parameters["angle"] = 30
Parameters["High Electrodes"] = 1
Parameters["High WG"] = 0.4
Parameters["High Slab"] = 0.2 
Parameters["High Substrate"] = 2
Parameters["Angle X"] = 0
Parameters["Angle Y"] = 90
Parameters["Angle Z"] = 0 


# Ganze Mach-Zehnder Modulator
Components.MZM(Parameters, proj)


# # Simple EO Phase Modulator
# Components.PhaseModulator(Parameters, proj)




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





# Ports on Electrodes for MZM Modulator
Parameters["Electrodes_Names"] = ["Electrode_Left:Electrode_Left1", "Electrode_Right:Electrode_Right1" , "Signal:Signal1"]
Parameters["Orientation"] = ["Positive", "Positive"]
Parameters["Coordinates"] = "Picks"
Parameters["Span"] = [[[3,3],[3,3]], [[3,3],[3,3]]]
Parameters["Potential"] = [1,2]
Parameters["Port Number"] = [1,2]
Parameters["Polarity"] = ["Positive", "Negative"]
Parameters["Face ID"] = [4,6]

VBA.WaveguidePorts_on_Electrodes_MZM(Parameters, proj)








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






# #Create Waveguide with Ports 
# Parameters = {}
# Parameters["Lenght WG"] = 5
# Parameters["Hight WG"] = 0.3
# Parameters["Width WG"] = 1
# Parameters["Substrate Height"] = 1
# Parameters["Slab Heigh"] = 0

# WG = Components.Squere_Waveguide(Parameters, proj)




# # Ports on Electrodes for Phase Modulator
# Parameters["Electrodes_Names"] = ["WG:WG1"]
# Parameters["Orientation"] = ["Positive", "Positive"]
# Parameters["Coordinates"] = "Picks"
# Parameters["Span"] = [[[3,3],[3,3]], [[3,3],[3,3]]]
# Parameters["Potential"] = [1,2]
# Parameters["Port Number"] = [1,2]
# Parameters["Polarity"] = ["Positive", "Positive"]
# Parameters["Face ID"] = [4,6]

# VBA.Optical_WaveguidePorts_MZM(Parameters, proj)








# # Clear Picks
# ClearPick = VBA.ClearAllPicks()
# proj.schematic.execute_vba_code(ClearPick, timeout=None)




# # Set Time Solver
# Parameters= {}
# Parameters["Accuracy"] = 35
# Parameters["Caclculate Modes Only"] = False
# Parameters["Auto Impedance"] = False
# Parameters["Impedance"] = 50
# Parameters["Source Port"] = 1


# Solver = VBA.SetTimeSolver(Parameters)
# proj.schematic.execute_vba_code(Solver, timeout=None)






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






# Set Mesh 
Parameters = {}
Parameters["Mesh Type"] = "PBA"
Parameters["Mesh Cells Near Object"] = 2
Parameters["Mesh Cells far Object"] = 2


Mesh = VBA.SetMesh(Parameters)
proj.schematic.execute_vba_code(Mesh, timeout=None)

    



# Set optical solver enviroment
# Units Properties
Parameters['Dimensions'] = "um"
# Parameters['Frequency']  = "THz"
Parameters['Frequency']  = "GHz"
Parameters['Time'] = "ns"
Parameters['Temperature'] = "degC"

# Set FreqeueWavelength Range 
# Parameters["Min Wavelength"] = 1.5
# Parameters["Max Wavelength"] = 1.6
Parameters["Min Frequency"] = 1
Parameters["Max Frequency"] = 150

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
Parameters["Mesh Cells Near Object"] = 4
Parameters["Mesh Cells far Object"] = 4


# Env = VBA.SetOpticalSimulationProperties(Parameters)
Env = VBA.SetElectricalSimulationProperties(Parameters)
proj.schematic.execute_vba_code(Env, timeout=None)

# Solver = VBA.ChangeSolverType("Time")
# proj.schematic.execute_vba_code(Solver, timeout=None)



# Set Time Solver
Parameters= {}
Parameters["Accuracy"] = 20
Parameters["Caclculate Modes Only"] = False
Parameters["Auto Impedance"] = True
Parameters["Impedance"] = 50
Parameters["Source Port"]  = 1
Parameters["Solver Mesh Type"] = "TLM"


Solver = VBA.SetTimeSolver(Parameters)
proj.schematic.execute_vba_code(Solver, timeout=None)


# # Set Freq Solver
# Parameters= {}
# Parameters["Accuracy"] = 30
# Parameters["Caclculate Modes Only"] = False
# Parameters["Auto Impedance"] = True
# Parameters["Impedance"] = 50
# Parameters["Source Port"]  = 1

# Solver = VBA.SetFreqSolver(Parameters)
# proj.schematic.execute_vba_code(Solver, timeout=None)


# # Set Monitors Optical Domain
# Parameters = {}
# Parameters["Wavelength"] = 1.55
# Parameters["Monitor Type"] = "Efield"
# Monitor = VBA.CreateEfieldMonitor(Parameters)
# proj.schematic.execute_vba_code(Monitor, timeout=None)


# # Set Monitors Electrical Domain
# Parameters = {}
# Parameters["Wavelength"] = 1.55
# Parameters["Monitor Type"] = "Efield"
# Monitor = VBA.CreateEfieldMonitor(Parameters)
# proj.schematic.execute_vba_code(Monitor, timeout=None)


# Parameters = {}
# Parameters["Wavelength"] = 1.55
# Parameters["Monitor Type"] = "Powerflow"
# Monitor = VBA.CreateEfieldMonitor(Parameters)
# proj.schematic.execute_vba_code(Monitor, timeout=None)



# Start Time Solver
Start = VBA.StartTimeSolver()
proj.schematic.execute_vba_code(Start, timeout=None)



# # Extract Results 
# Data = VBA.ExportResults("C:\\Users\\Martin\\Desktop\\CST_Project", "S21", "1D Results\\S-Parameters\\S1(1),2(1)")
# proj.schematic.execute_vba_code(Data, timeout=None)



# # Dim V As Variant
# # V = Monitor.GetNumberOfMonitors
# # Debug.Print V

# # Dim G As Variant
# # G = Monitor.GetSubVolumeSampling ("Powerflow (wl=1.55)")
# # Debug.Print G(0)
# # Debug.Print G(1)
# # Debug.Print G(2)



# # Sub Main()
# # Dim G As Variant
# # G = Port.GetModeType ("1", "1")
# # Debug.Print G
# # Debug.Print G
# # Debug.Print G
# # End Sub



# # def Port_Data(Parameters):

# #     PortNumber = Parameters["Port Number"]
# #     ModeNumber = Parameters["Mode Number"]

# #     Port = 'Sub Main () ' \
# #         '\n.GetFcutoff ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
# #         '\n.GetFrequency ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
# #         '\n.GetModeType ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
# #         '\n.GetBeta ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
# #         '\n.GetAlpha ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
# #         '\n.GetWaveImpedance ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
# #         '\n.GetLineImpedance ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
# #         '\n.GetNumberOfModes ' + '"' + str(PortNumber) + \
# #     '\nEnd Sub'
# #     data = ''.join(Port)
# #     return Port
    


# # Parameters = {}
# # Parameters["Port Number"] = 1
# # Parameters["Mode Number"] = 1

# # Data = Port_Data(Parameters)
# # proj.schematic.execute_vba_code(Data, timeout=None)





# # # save project into Patbh
# # proj.save(cst_project_path,allow_overwrite=  True)

# # #Close CST Project
# # proj.close()
# # # Close CST
# # mycst.close()



# # Sub Main
# # 	Pick.PickFaceFromId "Waveguide_Left:Waveguide_Left" "4"
# # End Sub'




# import numpy as np
# from scipy.special import fresnel
# import matplotlib.pyplot as plt


# def proper_euler_s_bend(Span_X, Span_Y, num_pts=1000):
#     """
#     Returns a smooth, continuous-curvature Euler S-bend from (0, 0) to (Span_X, Span_Y)
#     using Fresnel integrals and linear curvature.
#     """
#     # Arc length parameter
#     s = np.arange(-1, 1, 2/num_pts)
    
#     # Fresnel integrals (Euler spiral)
#     S, C = fresnel(s)

#     # Normalize to range [-1, 1]
#     C = C - C[0]
#     S = S - S[0]
#     C = C / (C[-1] - C[0])
#     S = S / (S[-1] - S[0])

#     # Euler spiral from -90° to +90° → makes an S-bend
#     x = C * Span_X
#     y = S * Span_Y

#     curve = np.vstack((y, x)).T

#     return curve



# # Example usage
# curve = proper_euler_s_bend(Span_X=100, Span_Y=20)

# plt.figure(figsize=(6, 3))
# plt.plot(curve[:, 0], curve[:, 1], label="Euler S-Bend")
# plt.axis("equal")
# plt.grid(True)
# plt.title("Proper Smooth Euler S-Bend")
# plt.xlabel("X (µm)")
# plt.ylabel("Y (µm)")
# plt.legend()
# plt.tight_layout()
# plt.show()