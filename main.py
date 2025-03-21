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
points = 1000



# Generate the Bezier and Cos points
ObjCurves = Curves(Lenght, Offset, points)
BezierCuve = ObjCurves.Bezier_Curve()
CosinusCurve = ObjCurves.Cosinus_Curve()



# plt.figure()
# plt.plot(BezierCuve[:,0], BezierCuve[:,1], color = "red", label = "Bezier Curve")
# plt.plot(CosinusCurve[:,0], CosinusCurve[:,1], color = "blue", label = "Cosinus Curve")
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
# cst_project2 = 'MZM_TestSim' # CST project
cst_project2 = 'Optical_Sim' # CST project
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

#
# Parameters = {}
# Parameters['X1'] = 0
# Parameters['Y1'] = 0
# Parameters['Z1'] = 0
# Parameters['X2'] = 5
# Parameters['Y2'] = 5
# Parameters['Z2'] = 0
#
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
# Parameters["Max Wavelength"] = 1.6

# data = VBA.SetSimWavelength(Parameters)
# proj.schematic.execute_vba_code(data, timeout=None)


# BondWire = BondWire(NameWire = "TestWire" ,Coordinates = Parameters, Height = 1, Radius = 0.5 , BondwireType = "Spline", Material = "Copper (annealed)",  NameFolder = "BondWire")
# proj.schematic.execute_vba_code(BondWire, timeout=None)
# ToSolid1 = ToSolid(SolidName = "TestSolid", CurveName = "TestWire", NameFolder = "BondWire", Material = "Copper (annealed)")
# proj.schematic.execute_vba_code(ToSolid1, timeout=None)


# Curve = Curve(CurveName = "TestCurve", Points = Points)
# proj.schematic.execute_vba_code(Curve, timeout=None)
# DataCurveWire = CurveWire(NameWire = "CurveWire", Radius = 0.5 , Points = Points, Material = "Copper (annealed)", CurveFolderName = "TestCurve", CurveName = "TestCurve")
# proj.schematic.execute_vba_code(DataCurveWire, timeout=None)
# ToSolid2 = ToSolid(SolidName = "TestSolidCurve", CurveName = "CurveWire", NameFolder = "CurveWire", Material = "Copper (annealed)")
# proj.schematic.execute_vba_code(ToSolid2, timeout=None)


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




# Create MZM Object
Parameters = {}
Parameters["Electrodes Lenght"] = 40
Parameters["Width GND"] = 5
Parameters["Width Signal"] = 2
Parameters["Width WG"] = 1
Parameters["Gap"] = 0.5
Parameters["angle"] = 30
Parameters["High Electrodes"] = 1
Parameters["High WG"] = 0.4
Parameters["High Slab"] = 0.2 
Parameters["High Substrate"] = 2


# Components.MZM(Parameters, proj)
Components.PhaseModulator(Parameters, proj)





# # Loop for go thrue electrodes and set ports
# Names = ["Electrode_Left:Electrode_Left1", "Electrode_Right:Electrode_Right1" , "Signal:Signal1"]
# # Ports = [[1,2], [3,4], [5,6]]
# PortNum = [1,2]
# # Pick Face
# for i in range(len(PortNum)):
#     Parameters["Orientation"] = ["Positive", "Positive"]
#     Parameters["Coordinates"] = "Picks"
#     Parameters["Span"] = [[[3,3],[3,3]], [[3,3],[3,3]]]
#     Parameters["Potential"] = [1,2]
#     Parameters["Port Number"] = PortNum
#     Parameters["Polarity"] = ["Positive", "Negative"]
#     Parameters["Face ID"] = [4,6]
    
#     for j in range(len(Names)):
#         PicParams = {}
#         PicParams["Option"] = "Face"
#         PicParams["Object"] = Names[j]
#         PicParams["Face Number"] = Parameters["Face ID"][i]
#         PickFace = VBA.Pick(PicParams)
#         Parameters["Solid Name"] = Names[j]
#         proj.schematic.execute_vba_code(PickFace, timeout=None)

#     Port = VBA.WaveguidePort(Parameters)
#     proj.schematic.execute_vba_code(Port[str(i+1)], timeout=None)


# for i in range(len(PortNum)):
#     Parameters["Distance"] = Parameters["Electrodes Lenght"] 
#     WidthObject = (2*Parameters["Width WG"] + Parameters["Width Signal"] + 2*Parameters["Width GND"] + 4*Parameters["Gap"])/2
#     Hight = (Parameters["High Slab"] + Parameters["High Slab"] + Parameters["High Electrodes"]) / 2
#     Parameters["Span"] = [[[-WidthObject - 3, -WidthObject - 3],[WidthObject + 3, WidthObject + 3]], [[-1.2, -1.2],[4.2,4.2]]]
#     Parameters["Port Number"] = PortNum

#     Port = VBA.MoveWaveguidePorts(Parameters)
#     proj.schematic.execute_vba_code(Port[str(i+1)], timeout=None)




# Loop for go thrue electrodes and set ports
Names = ["Electrode:Electrode1", "Signal:Signal1"]
# Ports = [[1,2], [3,4], [5,6]]
PortNum = [1,2]
# Pick Face
for i in range(len(PortNum)):
    Parameters["Orientation"] = ["Positive", "Positive"]
    Parameters["Coordinates"] = "Picks"
    Parameters["Span"] = [[[3,3],[3,3]], [[3,3],[3,3]]]
    Parameters["Potential"] = [1,2]
    Parameters["Port Number"] = PortNum
    Parameters["Polarity"] = ["Positive", "Negative"]
    Parameters["Face ID"] = [4,6]
    
    for j in range(len(Names)):
        PicParams = {}
        PicParams["Option"] = "Face"
        PicParams["Object"] = Names[j]
        PicParams["Face Number"] = Parameters["Face ID"][i]
        PickFace = VBA.Pick(PicParams)
        Parameters["Solid Name"] = Names[j]
        proj.schematic.execute_vba_code(PickFace, timeout=None)

    Port = VBA.WaveguidePort(Parameters)
    proj.schematic.execute_vba_code(Port[str(i+1)], timeout=None)


for i in range(len(PortNum)):
    Parameters["Distance"] = Parameters["Electrodes Lenght"] 
    WidthObject = (2*Parameters["Width WG"] + Parameters["Width Signal"] + 2*Parameters["Width GND"] + 4*Parameters["Gap"])/2
    Hight = (Parameters["High Slab"] + Parameters["High Slab"] + Parameters["High Electrodes"]) / 2
    Parameters["Span"] = [[[-WidthObject - 3, -WidthObject - 3],[WidthObject + 3, WidthObject + 3]], [[-1.2, -1.2],[4.2,4.2]]]
    Parameters["Port Number"] = PortNum

    Port = VBA.MoveWaveguidePorts(Parameters)
    proj.schematic.execute_vba_code(Port[str(i+1)], timeout=None)



# # # Clear Picks
# # ClearPick = VBA.ClearAllPicks()
# # proj.schematic.execute_vba_code(ClearPick, timeout=None)



# # Set discrete Ports 
# Parameters["Discrete Port Number"] = 3
# PortNum = [3,4,5,6]
# Parameters["Discrete Port Type"] = "Voltage"
# Parameters["Port Impedance"] = 50
# Parameters["Port Voltage"] = 2
# Parameters["Port Current"] = 1
# Parameters["Port Radius"] = 0
# Coordinates = {}
# # Coordinates["X1"] = Parameters["Electrodes Lenght"] / 2
# # Coordinates["X2"] = Parameters["Electrodes Lenght"] / 2
# # Coordinates["Y1"] = 0
# # Coordinates["Y2"] = Parameters["Width Signal"]/2 + 2*Parameters["Gap"] + Parameters["Width WG"] + Parameters["Width GND"]/2
# # Coordinates["Z1"] = Parameters["High Substrate"]/2 + Parameters["High Slab"] + (Parameters["High Electrodes"] / 2)
# # Coordinates["Z2"] = Parameters["High Substrate"]/2 + Parameters["High Slab"] + Parameters["High Electrodes"] / 2
# # Parameters["Discrete Port Coordinates"] = Coordinates


# Cor = {}
# Cor["X1"] = [Parameters["Electrodes Lenght"] / 2, Parameters["Electrodes Lenght"] / 2, -Parameters["Electrodes Lenght"] / 2, -Parameters["Electrodes Lenght"] / 2]
# Cor["X2"] = [Parameters["Electrodes Lenght"] / 2, Parameters["Electrodes Lenght"] / 2, -Parameters["Electrodes Lenght"] / 2, -Parameters["Electrodes Lenght"] / 2]
# Cor["Y1"] = [0, 0, 0, 0]
# Cor["Y2"] = [(Parameters["Width Signal"]/2 + 2*Parameters["Gap"] + Parameters["Width WG"] + Parameters["Width GND"]/2), -(Parameters["Width Signal"]/2 + 2*Parameters["Gap"] + Parameters["Width WG"] + Parameters["Width GND"]/2), (Parameters["Width Signal"]/2 + 2*Parameters["Gap"] + Parameters["Width WG"] + Parameters["Width GND"]/2), -(Parameters["Width Signal"]/2 + 2*Parameters["Gap"] + Parameters["Width WG"] + Parameters["Width GND"]/2)]
# Cor["Z1"] = [(Parameters["High Substrate"]/2 + Parameters["High Slab"] + (Parameters["High Electrodes"] / 2)), (Parameters["High Substrate"]/2 + Parameters["High Slab"] + (Parameters["High Electrodes"] / 2)), (Parameters["High Substrate"]/2 + Parameters["High Slab"] + (Parameters["High Electrodes"] / 2)), (Parameters["High Substrate"]/2 + Parameters["High Slab"] + (Parameters["High Electrodes"] / 2))]
# Cor["Z2"] = [(Parameters["High Substrate"]/2 + Parameters["High Slab"] + (Parameters["High Electrodes"] / 2)), (Parameters["High Substrate"]/2 + Parameters["High Slab"] + (Parameters["High Electrodes"] / 2)), (Parameters["High Substrate"]/2 + Parameters["High Slab"] + (Parameters["High Electrodes"] / 2)), (Parameters["High Substrate"]/2 + Parameters["High Slab"] + (Parameters["High Electrodes"] / 2))]



# for i in range(len(PortNum)):
#     Parameters["Discrete Port Number"] = PortNum[i]
#     Coordinates["X1"] = Cor["X1"][i]
#     Coordinates["X2"] = Cor["X2"][i]
#     Coordinates["Y1"] = Cor["Y1"][i]
#     Coordinates["Y2"] = Cor["Y2"][i]
#     Coordinates["Z1"] = Cor["Z1"][i]
#     Coordinates["Z2"] = Cor["Z2"][i]
#     Parameters["Discrete Port Coordinates"] = Coordinates

#     DiscretePort = VBA.SetDiscretePort(Parameters)
#     proj.schematic.execute_vba_code(DiscretePort, timeout=None)





# #Create Waveguide with Ports 
# Parameters = {}
# Parameters["Lenght WG"] = 5
# Parameters["Hight WG"] = 0.3
# Parameters["Width WG"] = 1
# Parameters["Substrate Height"] = 1
# Parameters["Slab Heigh"] = 0

# WG = Components.Squere_Waveguide(Parameters, proj)


# FaceID = [4,6]
# for index, i in enumerate(FaceID):
#     PicParams = {}
#     PicParams["Option"] = "Face"
#     PicParams["Object"] = "WG:WG1"
#     PicParams["Face Number"] = i
#     PickFace = VBA.Pick(PicParams)
#     proj.schematic.execute_vba_code(PickFace, timeout=None)




#     Parameters = {}
#     Parameters["Orientation"] = ["Positive", "Positive"]
#     Parameters["Coordinates"] = "Picks"
#     Parameters["Span"] = [[[3, 3],[3, 3]], [[3, 3],[3, 3]]]
#     Parameters["Potential"] = [1,2]
#     Parameters["Port Number"] = [1,2]
#     Parameters["Polarity"] = ["Positive", "Positive"]
#     Parameters["Solid Name"] = "WG:WG1"
#     Parameters["Face ID"] = [2,4]

#     Port = VBA.WaveguidePort(Parameters)
#     proj.schematic.execute_vba_code(Port[str(index + 1)], timeout=None)






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






# # Set Mesh 
# Parameters = {}
# Parameters["Mesh Type"] = "PBA"
# Parameters["Mesh Cells Near Object"] = 8
# Parameters["Mesh Cells far Object"] = 2


# Mesh = VBA.SetMesh(Parameters)
# proj.schematic.execute_vba_code(Mesh, timeout=None)


# '





# Set optical solver enviroment
# Units Properties
Parameters['Dimensions'] = "um"
Parameters['Frequency']  = "THz"
Parameters['Time'] = "ns"
Parameters['Temperature'] = "degC"

# Set FreqeueWavelength Range 
Parameters["Min Wavelength"] = 1.5
Parameters["Max Wavelength"] = 1.6
# Parameters["Min Frequency"] = 1
# Parameters["Max Frequency"] = 100

# Set Background
Parameters["Type Background"] = "Normal"
Parameters["Xmin Background"] = 5
Parameters["Xmax Background"] = 5
Parameters["Ymin Background"] = 5
Parameters["Ymax Background"] = 5
Parameters["Zmin Background"] = 5
Parameters["Zmax Background"] = 5


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
Parameters["Mesh Cells far Object"] = 1


Env = VBA.SetOpticalSimulationProperties(Parameters)
proj.schematic.execute_vba_code(Env, timeout=None)






# Set Time Solver
Parameters= {}
Parameters["Accuracy"] = 40
Parameters["Caclculate Modes Only"] = False
Parameters["Auto Impedance"] = True
Parameters["Impedance"] = 34
Parameters["Source Port"]  = 1


Solver = VBA.SetTimeSolver(Parameters)
proj.schematic.execute_vba_code(Solver, timeout=None)



# Set Monitors
Parameters = {}
Parameters["Wavelength"] = 1.55
Parameters["Monitor Type"] = "Efield"
Monitor = VBA.CreateEfieldMonitor(Parameters)
proj.schematic.execute_vba_code(Monitor, timeout=None)


Parameters = {}
Parameters["Wavelength"] = 1.55
Parameters["Monitor Type"] = "Powerflow"
Monitor = VBA.CreateEfieldMonitor(Parameters)
proj.schematic.execute_vba_code(Monitor, timeout=None)



# # # Start Time Solver
# # Start = VBA.StartTimeSolver()
# # proj.schematic.execute_vba_code(Start, timeout=None)


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