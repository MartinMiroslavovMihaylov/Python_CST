import numpy as np 
import matplotlib.pyplot as plt
import pandas as pandas
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
cst_path = r'/homes/lift/martinmi/Desktop/CST Python/Test_Object/' # path
cst_project2 = 'BondWire2' # CST project
cst_project_path = cst_path + cst_project2 + '.cst'


# Open CST as software
mycst =cst.interface.DesignEnvironment()


# Create Project in Project directory and name it
obj =cst.interface.DesignEnvironment.open_project(mycst,cst_project_path)

# cst.interface.Project.project_type(obj)




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
UnitParams['Dimensions'] = "mm"
UnitParams['Frequency'] = "ghz"
UnitParams['Time'] = "ms"
UnitParams['TemperatureUnit'] = "kelvin"
# UnitParams['Voltage'] = "V"
# UnitParams['Current'] = "A"
# UnitParams['Resistance'] = "Ohm"
# UnitParams['Conductance'] = "Siemens"
# UnitParams['Capacitance'] = "PikoF"
# UnitParams['Inductance'] = "NanoH"
#

Units = SetUnits(UnitParams)
obj.schematic.execute_vba_code(Units, timeout=None)





BondWire = BondWire(NameWire = "TestWire" ,Coordinates = Parameters, Height = 1, Radius = 0.5 , BondwireType = "Spline", Material = "Copper (annealed)",  NameFolder = "BondWire")
obj.schematic.execute_vba_code(BondWire, timeout=None)
ToSolid1 = ToSolid(SolidName = "TestSolid", CurveName = "TestWire", NameFolder = "BondWire", Material = "Copper (annealed)")
obj.schematic.execute_vba_code(ToSolid1, timeout=None)

# Curve = Curve(CurveName = "TestCurve", Points = Points)
# obj.schematic.execute_vba_code(Curve, timeout=None)
# DataCurveWire = CurveWire(NameWire = "CurveWire", Radius = 0.5 , Points = Points, Material = "Copper (annealed)", CurveFolderName = "TestCurve", CurveName = "TestCurve")
# obj.schematic.execute_vba_code(DataCurveWire, timeout=None)
# ToSolid2 = ToSolid(SolidName = "TestSolidCurve", CurveName = "CurveWire", NameFolder = "CurveWire", Material = "Copper (annealed)")
# obj.schematic.execute_vba_code(ToSolid2, timeout=None)

# #
# Parameters = {}
# Parameters['X1'] = -5
# Parameters['X2'] = 5
# Parameters['Y1'] = -5
# Parameters['Y2'] = 5
# Parameters['Z1'] = 0
# Parameters['Z2'] = 2
# TestBrick = Brick('TestBrick', Parameters)
# obj.schematic.execute_vba_code(TestBrick, timeout=None)
#
# # Add Global Parameters
# GlobalParam = {}
# GlobalParam['Length'] = 10
# GlobalParam['Width'] = 5
# GlobalParam['Hight'] = 2
#
# Globals  = AddGlobalParameter(GlobalParam)
# obj.schematic.execute_vba_code(Globals, timeout=None)
#
# # Delete Global Parameters
# GlobalParamDel = {}
# GlobalParamDel['Length'] = 10
# delete = DeleteGlobalParameter(GlobalParamDel)
# obj.schematic.execute_vba_code(delete, timeout=None)
#
#

# Define Background
Params = {}
Params["Xmin"] = 15
Params["Xmax"] = 15
Params["Ymin"] = 15
Params["Ymax"] = 15
Params["Zmin"] = 15
Params["Zmax"] = 15

BackObj = BackgroundSet(Params)
obj.schematic.execute_vba_code(BackObj, timeout=None)



# Set Boundary
BoundarySet = SetBoundary()
obj.schematic.execute_vba_code(BoundarySet, timeout=None)

#
# Port = WaveguidePort()
# obj.schematic.execute_vba_code(Port, timeout=None)



# save project into Patbh
cst.interface.Project.save(obj, cst_project_path, False)



#Close CST
cst.interface.Project.close(obj)
