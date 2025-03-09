import numpy as np 
import matplotlib.pyplot as plt
import pandas as pandas
import sys
import os
# Add the directory containing 'Curves_Functions' to sys.path
module_dir = os.path.abspath("C:/Users/marti/Desktop/UPB Kursen/CST with Python")
sys.path.append(module_dir)
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
cst_path = r'C:/Users/marti/Desktop/UPB Kursen/CST with Python/' # path
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








Parameters = {}
Parameters["Lenght_Electrodes"] = 40
Parameters["Lenght_Electrodes"] = 40
Parameters["Width_GND"] = 5
Parameters["Width_Signal"] = 2
Parameters["WG_Width"] = 1
Parameters["Gap"] = 0.5
Parameters["angle"] = 30
Parameters["High_Electrodes"] = 1
Parameters["High_WG"] = 0.4
Parameters["High_Slab"] = 0.2 
Parameters["High_Substrate"] = 2



MZM(Parameters, proj)





# Set Boundary
BoundaryParams = {}
BoundaryParams["Xmin"] = 'electric'
BoundaryParams["Xmax"] = 'electric'
BoundaryParams["Ymin"] = 'electric'
BoundaryParams["Ymax"] = 'electric'
BoundaryParams["Zmin"] = 'electric'
BoundaryParams["Zmax"] = 'electric' 
BoundaryParams["Xsymmetry"] = 'none'
BoundaryParams["Ysymmetry"] = 'none'
BoundaryParams["Zsymmetry"] = 'none'



BoundarySet = SetBoundary(BoundaryParams)
proj.schematic.execute_vba_code(BoundarySet, timeout=None)



# Set Waveguide Ports

PicParams = {}
PicParams["Option"] = "Face"
PicParams["Object"] = "Electrode_Left:Electrode_Left1"
PicParams["Face Number"] = 4

PickFace = Pick(PicParams)
proj.schematic.execute_vba_code(PickFace, timeout=None)




def WaveguidePortTest(Parameters):

    # Parameters to determin port position 
    Orientation = Parameters["Orientation"]
    Coordinates = Parameters["Coordinates"]
    Xmin = Parameters["Xmin"]
    Ymin = Parameters["Ymin"]
    Zmin = Parameters["Zmin"]   



    Port = 'Sub Main () ' \
                '\nWith Port' + \
                '\n.Reset' + \
                '\n.PortNumber (1)' + \
                '\n.NumberOfModes (2)' + \
                '\n.AdjustPolarization (False)' + \
                '\n.PolarizationAngle (0.0)' + \
                '\n.ReferencePlaneDistance (0)' + \
                '\n.Coordinates ' + '"' + str(Coordinates) + '"' + \
                '\n.Orientation ' +  '"' + str(Orientation) + '"' + \
                '\n.PortOnBound (False)' + \
                '\n.ClipPickedPortToBound (False)' + \
                '\n.Xrange ' + '"' + str(Xmin) + '"' + \
                '\n.Yrange ' + '"' + str(Ymin) + '"' + \
                '\n.Zrange ' + '"' + str(Zmin) + '"' + \
                '\n.Create' + \
            '\nEnd With'  + \
            '\nEnd Sub'

    return Port



Parameters = {}
Parameters["Orientation"] = "Positive"
Parameters["Coordinates"] = "Use Picks"
Parameters["Xmin"] = 5
Parameters["Ymin"] = 5
Parameters["Zmin"] = 5



Port = WaveguidePortTest(Parameters)
proj.schematic.execute_vba_code(Port, timeout=None)


# .AdjustPolarization (False)
# .PolarizationAngle (0.0)
# .ReferencePlaneDistance (0)
# .Orientation "Positive"
# .PortOnBound (False)
# .ClipPickedPortToBound (False)
# .Xrange "0", "5"
# .Yrange "0", "5"
# .Zrange "0", "5"



    




# # Material Definition
# Data = {}
# Data["X"] = 27
# Data["Y"] = 43
# Data["Z"] = 43

# # LiNbO3
# MaterialDAta = Material("LiNbO3", Data)
# proj.schematic.execute_vba_code(MaterialDAta, timeout=None)

# # Silicon
# dat = Material_Silicon("Silicon (lossy)")
# proj.schematic.execute_vba_code(dat, timeout=None)

# Parameters = {}
# Parameters['X1'] = -GlobalParam['Length']/2
# Parameters['X2'] = GlobalParam['Length']/2
# Parameters['Y1'] = TotalWidth/2
# Parameters['Y2'] = -TotalWidth/2 
# Parameters['Z1'] = 0
# Parameters['Z2'] = 1
# TestBrick = Brick('Substrate', Parameters, Material = "Silicon (lossy)")
# proj.schematic.execute_vba_code(TestBrick, timeout=None)



# Parameters = {}
# Parameters['X1'] = -GlobalParam['Length']/2
# Parameters['X2'] = GlobalParam['Length']/2
# Parameters['Y1'] = -TotalWidth/2
# Parameters['Y2'] = -TotalWidth/2 + GlobalParam['Width_GND']
# Parameters['Z1'] = 0
# Parameters['Z2'] = 2
# TestBrick = Brick('GND_Left_Electrode', Parameters)
# proj.schematic.execute_vba_code(TestBrick, timeout=None)



# Parameters = {}
# Parameters['X1'] = -GlobalParam['Length']/2
# Parameters['X2'] = GlobalParam['Length']/2
# Parameters['Y1'] = -TotalWidth/2 + GlobalParam['Width_GND'] + GlobalParam['Width_Gap']
# Parameters['Y2'] = -TotalWidth/2 + GlobalParam['Width_GND'] + GlobalParam['Width_Gap'] +  GlobalParam['Width_Signal']
# Parameters['Z1'] = 0
# Parameters['Z2'] = 2
# TestBrick = Brick('Signal_Electrode', Parameters)
# proj.schematic.execute_vba_code(TestBrick, timeout=None)



# Parameters = {}
# Parameters['X1'] = -GlobalParam['Length']/2
# Parameters['X2'] = GlobalParam['Length']/2
# Parameters['Y1'] = -TotalWidth/2 + GlobalParam['Width_GND'] + 2*GlobalParam['Width_Gap'] + GlobalParam['Width_Signal']
# Parameters['Y2'] = -TotalWidth/2 + 2*GlobalParam['Width_GND'] + 2*GlobalParam['Width_Gap'] + GlobalParam['Width_Signal']
# Parameters['Z1'] = 0
# Parameters['Z2'] = 2
# TestBrick = Brick('GND_right_Electrode', Parameters)
# proj.schematic.execute_vba_code(TestBrick, timeout=None)


# Points = {}
# Points["X"] = [GlobalParam['Length']/2, -GlobalParam['Length']/2]
# Points["Y"] = [-TotalWidth/2 + GlobalParam['Width_GND'], -TotalWidth/2 + GlobalParam['Width_GND'] + 1]
# Points["Z"] = [0,2]


# WGLeft_Edge = -TotalWidth/2 + GlobalParam['Width_GND']
# WGTopWidth = 4

# Points = {}
# Points["X"] = [-GlobalParam['Length']/2, -GlobalParam['Length']/2, GlobalParam['Length']/2, GlobalParam['Length']/2, -GlobalParam['Length']/2]
# # Points['X'] = np.array([-10, -10, 10, 10, -10])
# Points["Y"] = [WGLeft_Edge + WGTopWidth/2 , WGLeft_Edge + WGTopWidth, WGLeft_Edge + WGTopWidth, WGLeft_Edge + WGTopWidth/2, WGLeft_Edge + WGTopWidth/2]
# # Points['Y'] = np.array([-5, 5, 5, -5, -5])



# # Waveguide and Waveguide to solid
# WG = RibWG(WGName = "Waveguide_Left", Points = Points)
# proj.schematic.execute_vba_code(WG, timeout=None)
# RibWG_Test = RibWaveguide_ToSolid("Waveguide_Left", WG_Hight = 2, Angle = 10, WGFolderName = "Waveguide_Left", WGName = "Waveguide_Left", Material="LiNbO3")
# proj.schematic.execute_vba_code(RibWG_Test, timeout=None)



# WGLeft_Edge = TotalWidth/2 - GlobalParam['Width_GND']
# WGTopWidth = 4



# Points = {}
# Points["X"] = [-GlobalParam['Length']/2, -GlobalParam['Length']/2, GlobalParam['Length']/2, GlobalParam['Length']/2, -GlobalParam['Length']/2]
# # Points['X'] = np.array([-10, -10, 10, 10, -10])
# Points["Y"] = [-WGLeft_Edge - WGTopWidth/2 , -WGLeft_Edge - WGTopWidth, -WGLeft_Edge - WGTopWidth, -WGLeft_Edge - WGTopWidth/2, -WGLeft_Edge - WGTopWidth/2]
# # Points['Y'] = np.array([-5, 5, 5, -5, -5])




# WG = RibWG(WGName = "Waveguide_Right", Points = Points)
# proj.schematic.execute_vba_code(WG, timeout=None)
# RibWG_Test = RibWaveguide_ToSolid("Waveguide_Right", WG_Hight = -2, Angle = 10, WGFolderName = "Waveguide_Right", WGName = "Waveguide_Right", Material="LiNbO3")
# proj.schematic.execute_vba_code(RibWG_Test, timeout=None)












# # save project into Patbh
# proj.save(cst_project_path,allow_overwrite=  True)

# #Close CST Project
# proj.close()
# # Close CST
# mycst.close()