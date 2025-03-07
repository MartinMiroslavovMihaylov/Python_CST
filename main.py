import numpy as np 
import matplotlib.pyplot as plt
import pandas as pandas
import sys
import os
# Add the directory containing 'Curves_Functions' to sys.path
module_dir = os.path.abspath("C:/Users/Martin/Desktop/CST_Project")
sys.path.append(module_dir)
from Curves_Functions import Curves
from VBA_Test_Syntax import *




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





# # Define Background
# Params = {}
# Params["Xmin"] = 15
# Params["Xmax"] = 25
# Params["Ymin"] = 35
# Params["Ymax"] = 45
# Params["Zmin"] = 55
# Params["Zmax"] = 65

# BackObj = BackgroundSet(Params)
# proj.schematic.execute_vba_code(BackObj, timeout=None)



# # Set Boundary
# BoundarySet = SetBoundary()
# proj.schematic.execute_vba_code(BoundarySet, timeout=None)


# Port = WaveguidePort()
# proj.schematic.execute_vba_code(Port, timeout=None)








# MZM Design

Parameters = {}
Parameters['X1'] = -GlobalParam['Length']/2
Parameters['X2'] = GlobalParam['Length']/2
Parameters['Y1'] = -TotalWidth/2
Parameters['Y2'] = -TotalWidth/2 + GlobalParam['Width_GND']
Parameters['Z1'] = 0
Parameters['Z2'] = 2
TestBrick = Brick('GND_Left_Electrode', Parameters)
proj.schematic.execute_vba_code(TestBrick, timeout=None)



Parameters = {}
Parameters['X1'] = -GlobalParam['Length']/2
Parameters['X2'] = GlobalParam['Length']/2
Parameters['Y1'] = -TotalWidth/2 + GlobalParam['Width_GND'] + GlobalParam['Width_Gap']
Parameters['Y2'] = -TotalWidth/2 + GlobalParam['Width_GND'] + GlobalParam['Width_Gap'] +  GlobalParam['Width_Signal']
Parameters['Z1'] = 0
Parameters['Z2'] = 2
TestBrick = Brick('Signal_Electrode', Parameters)
proj.schematic.execute_vba_code(TestBrick, timeout=None)



Parameters = {}
Parameters['X1'] = -GlobalParam['Length']/2
Parameters['X2'] = GlobalParam['Length']/2
Parameters['Y1'] = -TotalWidth/2 + GlobalParam['Width_GND'] + 2*GlobalParam['Width_Gap'] + GlobalParam['Width_Signal']
Parameters['Y2'] = -TotalWidth/2 + 2*GlobalParam['Width_GND'] + 2*GlobalParam['Width_Gap'] + GlobalParam['Width_Signal']
Parameters['Z1'] = 0
Parameters['Z2'] = 2
TestBrick = Brick('GND_right_Electrode', Parameters)
proj.schematic.execute_vba_code(TestBrick, timeout=None)


Points = {}
Points["X"] = [GlobalParam['Length']/2, -GlobalParam['Length']/2]
Points["Y"] = [-TotalWidth/2 + GlobalParam['Width_GND'], -TotalWidth/2 + GlobalParam['Width_GND'] + 1]
Points["Z"] = [0,2]



def RibWaveguide(WGName, Points):

    WGName  = WGName
    # Extract the first points as starting Points
    Start_PointX = Points['X'][0]
    Start_PointY = Points['Y'][0]
    Start_PointZ = Points['Z'][0]
    tmp = []
    
    tmp.append('Sub Main () ' + '\nWith Polygon3D' + '\n.Reset' + '\n.Name ' + '"' + WGName  + '"' + '\n.Curve ' + '"' + WGName  + '"' + '\n.Point ' + '"' + str(Start_PointX) + '"' + ',' + '"' + str(Start_PointY) + '"'',' + '"' + str(Start_PointZ) + '"')
    for i in range(1, len(Points['X'])):
        b = '\n.LineTo ' + '"' + str(Points['X'][i]) + '"' + ',' + '"' + str(Points['Y'][i]) + '"'',' + '"' + str(Points['Z'][i]) + '"'
        tmp.append(b)
    tmp.append('\n.Create' + '\nEnd With' + '\nEnd Sub')
    CurveData = ''.join(tmp)

    return CurveData



Sub Main
With Polygon3D
.Reset
.Name "WG"
.Curve "WG"
.Point "40.0","-7.5","0"
.Point "40.0","-7.4","2"
.Point "40.0","-6.5","2"
.Point "40.0","-6.4","0"
.Point "-40.0","-6.4","0"
.Point "-40.0","-6.5","2"
.Point "-40.0","-7.4","2"
.Point "-40.0","-7.5","0"
.Create
End With
End Sub


RibWG = RibWaveguide("WG",Points)
proj.schematic.execute_vba_code(RibWG, timeout=None)

ToSolid2 = ToSolid(SolidName = "TestSolidWG", CurveName = "WG", NameFolder = "WG", Material = "Copper (annealed)")

# # save project into Patbh
# proj.save(cst_project_path,allow_overwrite=  True)

# #Close CST Project
# proj.close()
# # Close CST
# mycst.close()