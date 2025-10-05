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
import numpy as np



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






from CST_Python_Wrapper import CST_Commands

obj = CST_Commands()
obj.New_Project("MWS")
# obj.Save_Project("C:/Users/marti/Desktop/UPB Kursen/CST with Python/CST_Bonddrahtmodell", "Wrapper_Save", False)
# obj.Open_Project("C:/Users/marti/Desktop/UPB Kursen/CST with Python/CST_Bonddrahtmodell/Wrapper_Save.cst")



# Parameters = {}
# Parameters["Brick Lenght Max"] = 10/2
# Parameters["Brick Lenght Min"] = -10/2
# Parameters["Brick Width Max"] = 4/2
# Parameters["Brick Width Min"] = -4/2
# Parameters["Brick Hight Max"] = 2/2
# Parameters["Brick Hight Min"] = -2/2
# Parameters["Brick Name"] = "Brick Test"
# Parameters["Component Name"] = "Brick Component"
# Parameters["Material"] = "Vacuum"
# obj.Brick(Parameters)


Parameters = {}
Parameters["Unit Lenght"] = "um"
Parameters["Unit Frequency"] = "GHz"
Parameters["Unit Time"] = "ns"
Parameters["Unit Temperature"] = "K"
obj.set_Units(Parameters)



# obj.add_Si()
# obj.add_SiO2()
# obj.add_Au()
obj.add_material("Si")




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



# obj.Curve("Curve_Python", Points)
# obj.ToSolid("Curve_python_Solid", CurveName = "Curve_Python", NameFolder = "Curve_Python", Material = "Si" )
# obj.Poligon_2D()
# obj.RibWaveguide_ToSolid()


# GlobalParam = {}
# GlobalParam['Length'] = 80
# GlobalParam['Width_GND'] = 10
# GlobalParam['Width_Signal'] = 5
# GlobalParam['Width_Gap'] = 5
# GlobalParam['Hight'] = 2

# obj.AddGlobalParameter(GlobalParam)


# params = {}
# params['Width_Gap'] = 5
# params['Hight'] = 2
# obj.DeleteGlobalParameter(params)




Params = {}
Params["Type Background"] = "Normal"
Params["Xmin Background"] = 6
Params["Xmax Background"] = 5
Params["Ymin Background"] = 7
Params["Ymax Background"] = 8
Params["Zmin Background"] = 9
Params["Zmax Background"] = 11

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



Parameters = {}
Parameters["Min Frequency"] = 0     
Parameters["Max Frequency"] = 200
obj.setSimFreqeuncy(Parameters)




# # Set Time Solver
# Parameters= {}
# Parameters["Accuracy"] = 30
# Parameters["Caclculate Modes Only"] = False
# Parameters["Auto Impedance"] = True
# Parameters["Impedance"] = 50
# Parameters["Source Port"]  = 1
# Parameters["Solver Mesh Type"] = "TLM"
# obj.setTimeSolver(Parameters)



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
# Parameters["High Substrate"] = 1
# Parameters["Angle X"] = 0
# Parameters["Angle Y"] = 90
# Parameters["Angle Z"] = 0 

# obj.MZM(Parameters)




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

# obj.PhaseModulator(Parameters)


# #Create Waveguide with Ports 
# Parameters = {}
# Parameters["Lenght WG"] = 5
# Parameters["Hight WG"] = 0.3
# Parameters["Width WG"] = 1
# Parameters["Substrate Height"] = 1
# Parameters["Slab Heigh"] = 0.3
# obj.Squere_Waveguide(Parameters)






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


# obj.BondWire(NameWire = "TestWire" ,Coordinates = Parameters, Height = 1, Radius = 0.5 , BondwireType = "Spline", Material = "Copper (annealed)",  NameFolder = "BondWire")







# obj.Curve(CurveName = "TestCurve", Points = Points)
# obj.CurveWire(NameWire = "CurveWire2", Radius = 0.5 , Points = Points, Material = "PEC", CurveFolderName = "TestCurve", CurveName = "TestCurve")
# obj.ToSolid(SolidName = "TestSolidCurve2", CurveName = "CurveWire2", NameFolder = "CurveWire2", Material = "Si")




obj.setDomainSolverType("Freq")




# Parameters = {}
# Parameters["Cylinder Name"] = "Test Cylinder"
# Parameters["Component Name"] = "Cylinder"
# Parameters["Material"] = "Au"
# Parameters["Outer Radius"] = 5
# Parameters["Inner Radius"] = 2
# Parameters["Orentation Axis"] = "X"
# Parameters["X min"] = 0
# Parameters["X max"] = 10
# Parameters["Z center"] = 0
# Parameters["Y center"] = 0

# obj.add_Au()
# obj.Cylinder(Parameters)


# Parameters = {}
# Parameters["Cylinder Name"] = "Test Cylinder Y"
# Parameters["Component Name"] = "Cylinder Y"
# Parameters["Material"] = "Au"
# Parameters["Outer Radius"] = 5
# Parameters["Inner Radius"] = 2
# Parameters["Orentation Axis"] = "Y"
# Parameters["Y min"] = 0
# Parameters["Y max"] = 10
# Parameters["X center"] = 0
# Parameters["Z center"] = 0


# obj.Cylinder(Parameters)



# Parameters = {}
# Parameters["Cylinder Name"] = "Test Cylinder Z2"
# Parameters["Component Name"] = "Cylinder Z"
# Parameters["Material"] = "Au"
# Parameters["Outer Radius"] = 5
# Parameters["Inner Radius"] = 2
# Parameters["Orentation Axis"] = "Z"
# Parameters["Z min"] = 0
# Parameters["Z max"] = 10
# Parameters["Y center"] = 0
# Parameters["X center"] = 0


# obj.Cylinder(Parameters)




number_wires = 6
gap_pairs = 5
gap_wires = 1
Names = [["Sig", "Sig"], ["LeftGND", "LeftGND"], ["RightGND", "RightGND"]]
Names_Components = [ "Sig", "LeftGND","RightGND"]
Wire_Radius = 1
Lenght_Wire = 100
Dis_to_GND = 2
GND_Thickness = 1
Clad_Thickness = 2




Parameters = {}
obj.add_Au()
Parameters["Material"] = "Au"
Parameters["Outer Radius"] = Wire_Radius
Parameters["Inner Radius"] = 0
Parameters["Orentation Axis"] = "Z"
Parameters["Z min"] = 0
Parameters["Z max"] = Lenght_Wire
Parameters["Y center"] = Parameters["Outer Radius"] 


midL = -(Parameters["Outer Radius"] + gap_wires/2)
midR =  (Parameters["Outer Radius"] + gap_wires/2)

PosX = [midL, midR]
Tranlsate_pos = [0, -(2*Parameters["Outer Radius"]*2 + gap_wires + gap_pairs ) , (2*Parameters["Outer Radius"]*2 + gap_wires + gap_pairs)]

for i in range(len(Names_Components)):
    Parameters["Component Name"] = Names_Components[i]
    for j in range(2):
        Parameters["Cylinder Name"] = Names[i][j] + str(j)
        Parameters["X center"] = PosX[j]
        obj.Cylinder(Parameters)
    Trans = {}
    Trans["Name Object"] = Names_Components[i]
    Trans["Position X"] = Tranlsate_pos[i]
    Trans["Position Y"] = 0
    Trans["Position Z"] = 0
    Trans["Structure Type"] = "Shape"
    obj.Translation(Trans)


# Create metal bot 
Parameters["Brick Lenght Max"] = (max(Tranlsate_pos) * 2) + (Parameters["Outer Radius"]*2 + gap_wires)*2
Parameters["Brick Lenght Min"] = (min(Tranlsate_pos) * 2) - (Parameters["Outer Radius"]*2 + gap_wires)*2
Parameters["Brick Width Max"] =   - Dis_to_GND*2 - GND_Thickness*2
Parameters["Brick Width Min"] =  - Dis_to_GND*2 
Parameters["Brick Hight Max"] = Lenght_Wire*2 
Parameters["Brick Hight Min"] = 0
Parameters["Brick Name"] = "GND Plate"
Parameters["Component Name"] = "GND"
Parameters["Material"] = Parameters["Material"]
obj.Brick(Parameters)



# Create Cladding
Parameters["Brick Lenght Max"] = (max(Tranlsate_pos) * 2) + (Parameters["Outer Radius"]*2 + gap_wires)*2
Parameters["Brick Lenght Min"] = (min(Tranlsate_pos) * 2) - (Parameters["Outer Radius"]*2 + gap_wires)*2
Parameters["Brick Width Max"] =  Parameters["Outer Radius"]*2 +  Clad_Thickness*2
Parameters["Brick Width Min"] =  - Dis_to_GND*2 
Parameters["Brick Hight Max"] = Lenght_Wire*2 
Parameters["Brick Hight Min"] = 0
Parameters["Brick Name"] = "Cladding_1"
Parameters["Component Name"] = "Cladding"
Parameters["Material"] = "Vacuum"
obj.Brick(Parameters)




# Pick faces for solvers 

Names = [["Sig", "Sig"], ["LeftGND", "LeftGND"], ["RightGND", "RightGND"]]
Names_Components = [ "Sig", "LeftGND","RightGND"]
Parameters["Port Number"] = [1 ,2]
Parameters["Face ID"] = [1 ,3]


Parameters["Orientation"] = "Positive"
Parameters["Coordinates"] = "Picks"
Parameters["Port Number"]
Parameters["Span"] = [5, 5, 0]
Parameters["Picked Port Polarity"] = "positive", "positive", "negative", "negative", "negative", "negative"
Parameters["Picked Component Name"] = []
Parameters["Picked Port Number"]  = 1
Parameters["Number of picks"] = 6


for i in range(len(Names_Components)):
    for j in range(len(Parameters["Port Number"])):
        PicParams = {}
        PicParams["Option"] = "Face"
        PicParams["Face Number"] = Parameters["Face ID"][j]
        # PicParams["Object"] = f"""{Names_Components[i]}:{Names[i][j]+str(j)}"""
        for k in range(len(Names[0])):
            PicParams["Object"] = f"""{Names_Components[i]}:{Names[i][j]+str(k)}"""
            Parameters["Picked Component Name"].append(PicParams["Object"])
            PickFace = obj.Pick(PicParams)

obj.WaveguidePortWithPins(Parameters)
        
# lines = []
# if Parameters["Number of picks"] > 2:
#     for i in range(len(Parameters["Picked Port Polarity"])):
#         lines.append(f'.AddPotentialPicked "{1}", "{Parameters["Picked Port Polarity"][i]}", "{Parameters["Picked Component Name"][i]}", "1"')

#     # join them with newlines
#     line_block = "\n".join(lines)

#     vba_code = f"""
#                 With Port
#                 .Reset
#                 .PortNumber "{PortNumber[0]}"
#                 .NumberOfModes "5"
                
                
#                 .Orientation "{Parameters["Orientation"]}"
#                 .Coordinates "{Parameters["Coordinates"]}"
#                 .PortOnBound "False"
#                 .ClipPickedPortToBound "False"
#                 .XrangeAdd "{Parameters["Span"][0]}", "{Parameters["Span"][0]}"
#                 .YrangeAdd "{Parameters["Span"][1]}", "{Parameters["Span"][1]}"
#                 .ZrangeAdd "{Parameters["Span"][2]}", "{Parameters["Span"][2]}"
#                 .AdjustPolarization "True"
#                 .PolarizationAngle "0"
#                 .SingleEnded "False"
#                 {line_block}
#                 .Create
#                 End With
#                 """
#     obj.prj.model3d.add_to_history("create waveguide port", vba_code)
