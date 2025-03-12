import numpy as np 
import matplotlib.pyplot as plt
import pandas as pandas
import sys
import os
# Add the directory containing the project to sys.path
current_path = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_path)
from Curves_Functions import Curves
from VBA_Test_Syntax import *
from Components import *



# Define Curves Parameters and Data
Lenght = 100
Offset = 40
points = 1000



# Generate the Bezier and Cos points
ObjCurves = Curves(Lenght, Offset, points)
BezierCuve = ObjCurves.Bezier_Curve()
CosinusCurve = ObjCurves.Cosinus_Curve()



plt.figure()
plt.plot(BezierCuve[:,0], BezierCuve[:,1], color = "red", label = "Bezier Curve")
plt.plot(CosinusCurve[:,0], CosinusCurve[:,1], color = "blue", label = "Cosinus Curve")

plt.xlabel("Length S-Bend/$\mu m$")
plt.ylabel("Offset S-Bends/ $\mu m$")
plt.legend(loc = "best")
plt.grid()
plt.show()



# Call CST as Object and Open existing Project 

# Local path to CST project file --> Please adapt
cst_path = r'C:/Users/Martin/Desktop/CST_Project/' # path
cst_project2 = 'CST_ProgrammTest' # CST project
cst_project_path = cst_path + cst_project2 + '.cst'



import sys
sys.path.append(r"C:/Program Files (x86)/CST Studio Suite 2025/AMD64/python_cst_libraries")
import cst
print(cst.__file__)
# Open CST as software
from cst.interface import DesignEnvironment

mycst = DesignEnvironment.new()
mycst.quiet_mode_enabled()
# mycst =cst.interface.DesignEnvironment()

# Create Project in Project directory and name it
proj = mycst.open_project(cst_project_path)


Parameters = {}
Parameters['X1'] = 0
Parameters['Y1'] = 0
Parameters['Z1'] = 0
Parameters['X2'] = 5
Parameters['Y2'] = 5
Parameters['Z2'] = 0

Points = {}
x = []
y = []
for i in range(0, 100):
    x.append(i)
    y.append(i*4)
x = np.array(x)
y = np.array(y)
Points['X'] = x
Points['Y'] = y


# Set Units
UnitParams = {}
UnitParams['Dimensions'] = "um"
UnitParams['Frequency'] = "GHz"
UnitParams['Time'] = "ns"
UnitParams['Temperature'] = "degC"
# UnitParams['Voltage'] = "V"
# UnitParams['Current'] = "A"
# UnitParams['Resistance'] = "Ohm"
# UnitParams['Conductance'] = "Siemens"
# UnitParams['Capacitance'] = "PikoF"
# UnitParams['Inductance'] = "NanoH"
#


# def SetUnits():
#     Units = 'Sub Main () ' \
#                     '\n With Units' + \
#                         '\n.Frequency "GHz"' + \
#                         '\n.Time "ms"' + \
#                     '\n End With' + \
#                 '\n End Sub'
#     return Units

Units = SetUnits(UnitParams)
proj.schematic.execute_vba_code(Units, timeout=None)




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


# #
# Parameters = {}
# Parameters['X1'] = -5
# Parameters['X2'] = 5
# Parameters['Y1'] = -5
# Parameters['Y2'] = 5
# Parameters['Z1'] = 0
# Parameters['Z2'] = 2
# TestBrick = Brick('TestBrick', Parameters)
# proj.schematic.execute_vba_code(TestBrick, timeout=None)




# Add Global Parameters
GlobalParam = {}
GlobalParam['Length'] = 80
GlobalParam['Width_GND'] = 10
GlobalParam['Width_Signal'] = 5
GlobalParam['Width_Gap'] = 5
GlobalParam['Hight'] = 2
TotalWidth = 2*GlobalParam['Width_GND'] + GlobalParam['Width_Signal'] + 2*GlobalParam['Width_Gap']

Globals  = AddGlobalParameter(GlobalParam)
proj.schematic.execute_vba_code(Globals, timeout=None)


# # Delete Global Parameters
# GlobalParamDel = {}
# GlobalParamDel['Length'] = 10
# delete = DeleteGlobalParameter(GlobalParamDel)
# proj.schematic.execute_vba_code(delete, timeout=None)





# Define Background
Params = {}
Params["Xmin"] = 1
Params["Xmax"] = 1
Params["Ymin"] = 1
Params["Ymax"] = 1
Params["Zmin"] = 1
Params["Zmax"] = 1

BackObj = BackgroundSet(Params)
proj.schematic.execute_vba_code(BackObj, timeout=None)





# Port = WaveguidePort()
# proj.schematic.execute_vba_code(Port, timeout=None)





# Parameters = {}
# Parameters["Lenght_Electrodes"] = 40
# Parameters["Lenght_Electrodes"] = 40
# Parameters["Width_GND"] = 5
# Parameters["Width_Signal"] = 2
# Parameters["WG_Width"] = 1
# Parameters["Gap"] = 0.5
# Parameters["angle"] = 30
# Parameters["High_Electrodes"] = 1
# Parameters["High_WG"] = 0.4
# Parameters["High_Slab"] = 0.2 
# Parameters["High_Substrate"] = 2



# MZM(Parameters, proj)





# # Set Boundary
# BoundaryParams = {}
# BoundaryParams["Xmin"] = 'open'
# BoundaryParams["Xmax"] = 'open'
# BoundaryParams["Ymin"] = 'open'
# BoundaryParams["Ymax"] = 'open'
# BoundaryParams["Zmin"] = 'open'
# BoundaryParams["Zmax"] = 'open' 
# BoundaryParams["Xsymmetry"] = 'none'
# BoundaryParams["Ysymmetry"] = 'none'
# BoundaryParams["Zsymmetry"] = 'none'



# BoundarySet = SetBoundary(BoundaryParams)
# proj.schematic.execute_vba_code(BoundarySet, timeout=None)






# # Loop for go thrue electrodes and set ports
# Names = ["Electrode_Left:Electrode_Left1", "Electrode_Right:Electrode_Right1" , "Signal:Signal1"]
# Ports = [[1,2], [3,4], [5,6]]
# for index, j in enumerate(Names, start=1):
#     # print(Names[index-1])
#     Parameters = {}
#     Parameters["Orientation"] = ["Positive", "Negative"]
#     Parameters["Coordinates"] = "Picks"
#     Parameters["Span"] = [[[1,1],[0.5,0.5]], [[1,1],[1,1]]]
#     Parameters["Potential"] = [1,2]
#     Parameters["Port Number"] = Ports[index-1]
#     Parameters["Polarity"] = ["Positive", "Negative"]
#     Parameters["Solid Name"] = Names[index-1]
#     Parameters["Face ID"] = [4,6]
#     PortNum = [1,2]

#     PicParams = {}
#     PicParams["Option"] = "Face"
#     PicParams["Object"] = Names[index-1]
#     Port = WaveguidePort(Parameters)
    
#     for i in range(len(Parameters["Face ID"])):
#         PicParams["Face Number"] = Parameters["Face ID"][i]
#         PickFace = Pick(PicParams)
#         proj.schematic.execute_vba_code(PickFace, timeout=None)
#         proj.schematic.execute_vba_code(Port[str(i+1)], timeout=None)




# # Loop for go thrue electrodes and set ports
# Names = ["Waveguide_Left:Waveguide_Left", "Waveguide_Right:Waveguide_Right"]
# Ports = [[7, 8], [9, 10]]
# for index, j in enumerate(Names, start=1):
#     # print(Names[index-1])
#     Parameters = {}
#     Parameters["Orientation"] = ["Positive", "Negative"]
#     Parameters["Coordinates"] = "Picks"
#     Parameters["Span"] = [[[0.3, 0.3],[0.3, 0.3]], [[0.3, 0.3],[0.3, 0.3]]]
#     Parameters["Potential"] = [1,2]
#     Parameters["Port Number"] = Ports[index-1]
#     Parameters["Polarity"] = ["Positive", "Negative"]
#     Parameters["Solid Name"] = Names[index-1]
#     Parameters["Face ID"] = [2,4]
#     PortNum = [1,2]

#     PicParams = {}
#     PicParams["Option"] = "Face"
#     PicParams["Object"] = Names[index-1]
#     Port = WaveguidePort(Parameters)
    
#     for i in range(len(Parameters["Face ID"])):
#         PicParams["Face Number"] = Parameters["Face ID"][i]
#         PickFace = Pick(PicParams)
#         print(PortNum[i])
#         proj.schematic.execute_vba_code(PickFace, timeout=None)
#         proj.schematic.execute_vba_code(Port[str(i+1)], timeout=None)




Parameters = {}
Parameters["Lenght WG"] = 5
Parameters["Hight WG"] = 0.3
Parameters["Width WG"] = 1
Parameters["Substrate Height"] = 1
Parameters["Slab Heigh"] = 0

WG = Squere_Waveguide(Parameters, proj)


PicParams = {}
PicParams["Option"] = "Face"
PicParams["Object"] = "WG:WG1"
PicParams["Face Number"] = "6"
PickFace = Pick(PicParams)
proj.schematic.execute_vba_code(PickFace, timeout=None)




Parameters = {}
Parameters["Orientation"] = ["Positive", "Negative"]
Parameters["Coordinates"] = "Picks"
Parameters["Span"] = [[[0.5, 0.5],[0.5, 0.5]], [[0.5, 0.5],[0.5, 0.5]]]
Parameters["Potential"] = [1,2]
Parameters["Port Number"] = [1,2]
Parameters["Polarity"] = ["Positive", "Positive"]
Parameters["Solid Name"] = "WG:WG1"
Parameters["Face ID"] = [2,4]

Port = OpticalWaveguidePort(Parameters)
proj.schematic.execute_vba_code(Port['2'], timeout=None)







def Port_Data(Parameters):

    PortNumber = Parameters["Port Number"]
    ModeNumber = Parameters["Mode Number"]

    Port = 'Sub Main () ' \
        '\n.GetFcutoff ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetFrequency ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetModeType ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetBeta ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetAlpha ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetWaveImpedance ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetLineImpedance ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetNumberOfModes ' + '"' + str(PortNumber) + \
    '\nEnd Sub'
    data = ''.join(Port)
    return Port
    


Parameters = {}
Parameters["Port Number"] = 1
Parameters["Mode Number"] = 1

Data = Port_Data(Parameters)
proj.schematic.execute_vba_code(Data, timeout=None)





# # save project into Patbh
# proj.save(cst_project_path,allow_overwrite=  True)

# #Close CST Project
# proj.close()
# # Close CST
# mycst.close()



# Sub Main
# 	Pick.PickFaceFromId "Waveguide_Left:Waveguide_Left" "4"
# End Sub