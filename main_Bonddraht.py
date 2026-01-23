import numpy as np 
import matplotlib.pyplot as plt
import pandas as pandas
import sys
import os
# Add the directory containing the project to sys.path
current_path = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_path)
# from Curves_Functions import Curves
# import Functions as VBA
# import Components as Components
from CST_Python_Wrapper import CST_Commands, Curves
import numpy as np


# Define Curves Parameters and Data
Lenght = 100
Offset = 40
points = 100
# Generate the Bezier and Cos points
ObjCurves = Curves(Lenght, Offset, points)
CosinusCurve = ObjCurves.Cosinus_Curve()


# # Define Curves Parameters and Data
# Lenght = 100
# Offset = 40
# points = 100



# # Generate the Bezier and Cos points
# ObjCurves = Curves(Lenght, Offset, points)
# BezierCuve = ObjCurves.Bezier_Curve()
# CosinusCurve = ObjCurves.Cosinus_Curve()
# EulerCurve = ObjCurves.Euler_Curve()


# plt.figure()
# plt.plot(BezierCuve[:,0], BezierCuve[:,1], color = "red", label = "Bezier Curve")
# plt.plot(CosinusCurve[:,0], CosinusCurve[:,1], color = "blue", label = "Cosinus Curve")
# plt.plot(EulerCurve[:,0], EulerCurve[:,1], color = "green", label = "Euler Curve")
# plt.xlabel("Length S-Bend/$\mu m$")
# plt.ylabel("Offset S-Bends/ $\mu m$")
# plt.legend(loc = "best")
# plt.grid()
# plt.show()






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



Parameters = {}
Parameters["Min Frequency"] = 0     
Parameters["Max Frequency"] = 200
obj.setSimFreqeuncy(Parameters)


# Set Frequency Solver
Parameters = {}

obj.setFreqSolver(Parameters)


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




# obj.setDomainSolverType("Freq")




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


# ###################################################################################
# # Cylinder Electrodes
# ###################################################################################


# number_wires = 6
# gap_pairs = 5
# gap_wires = 1
# Names = [["Sig", "Sig"], ["LeftGND", "LeftGND"], ["RightGND", "RightGND"]]
# Names_Components = [ "Sig", "LeftGND","RightGND"]
# Wire_Radius = 1
# Lenght_Wire = 100
# Dis_to_GND = 2
# GND_Thickness = 1
# Clad_Thickness = 2




# Parameters = {}
# obj.add_Au()
# Parameters["Material"] = "Au"
# Parameters["Outer Radius"] = Wire_Radius
# Parameters["Inner Radius"] = 0
# Parameters["Orentation Axis"] = "Z"
# Parameters["Z min"] = 0
# Parameters["Z max"] = Lenght_Wire
# Parameters["Y center"] = Parameters["Outer Radius"] 


# midL = -(Parameters["Outer Radius"] + gap_wires/2)
# midR =  (Parameters["Outer Radius"] + gap_wires/2)

# PosX = [midL, midR]
# Tranlsate_pos = [0, -(2*Parameters["Outer Radius"]*2 + gap_wires + gap_pairs ) , (2*Parameters["Outer Radius"]*2 + gap_wires + gap_pairs)]

# for i in range(len(Names_Components)):
#     Parameters["Component Name"] = Names_Components[i]
#     for j in range(2):
#         Parameters["Cylinder Name"] = Names[i][j] + str(j)
#         Parameters["X center"] = PosX[j]
#         obj.Cylinder(Parameters)
#     Trans = {}
#     Trans["Name Object"] = Names_Components[i]
#     Trans["Position X"] = Tranlsate_pos[i]
#     Trans["Position Y"] = 0
#     Trans["Position Z"] = 0
#     Trans["Structure Type"] = "Shape"
#     Trans["Translate Type"] = "Translate"
#     obj.Translation(Trans)


# # Create metal bot 
# Parameters["Brick Lenght Max"] = (max(Tranlsate_pos) * 2) + (Parameters["Outer Radius"]*2 + gap_wires)*2
# Parameters["Brick Lenght Min"] = (min(Tranlsate_pos) * 2) - (Parameters["Outer Radius"]*2 + gap_wires)*2
# Parameters["Brick Width Max"] =   - Dis_to_GND*2 - GND_Thickness*2
# Parameters["Brick Width Min"] =  - Dis_to_GND*2 
# Parameters["Brick Hight Max"] = Lenght_Wire*2 
# Parameters["Brick Hight Min"] = 0
# Parameters["Brick Name"] = "GND Plate"
# Parameters["Component Name"] = "GND"
# Parameters["Material"] = Parameters["Material"]
# obj.Brick(Parameters)



# # Create Cladding
# Parameters["Brick Lenght Max"] = (max(Tranlsate_pos) * 2) + (Parameters["Outer Radius"]*2 + gap_wires)*2
# Parameters["Brick Lenght Min"] = (min(Tranlsate_pos) * 2) - (Parameters["Outer Radius"]*2 + gap_wires)*2
# Parameters["Brick Width Max"] =  Parameters["Outer Radius"]*2 +  Clad_Thickness*2
# Parameters["Brick Width Min"] =  - Dis_to_GND*2 
# Parameters["Brick Hight Max"] = Lenght_Wire*2 
# Parameters["Brick Hight Min"] = 0
# Parameters["Brick Name"] = "Cladding_1"
# Parameters["Component Name"] = "Cladding"
# Parameters["Material"] = "Vacuum"
# obj.Brick(Parameters)




# # Pick faces for solvers 

# Names = [["Sig", "Sig"], ["LeftGND", "LeftGND"], ["RightGND", "RightGND"]]
# Names_Components = [ "Sig", "LeftGND","RightGND"]
# Parameters["Port Number"] = [1 ,2]
# Parameters["Face ID"] = [1 ,3]


# Parameters["Orientation"] = "Positive"
# Parameters["Coordinates"] = "Picks"
# Parameters["Port Number"]
# Parameters["Span"] = [5, 5, 0]
# Parameters["Picked Port Polarity"] = "positive", "positive", "negative", "negative", "negative", "negative"
# Parameters["Picked Component Name"] = []
# Parameters["Picked Port Number"]  = 1
# Parameters["Number of picks"] = 6


# for i in range(len(Names_Components)):
#     for j in range(len(Parameters["Port Number"])):
#         PicParams = {}
#         PicParams["Option"] = "Face"
#         PicParams["Face Number"] = Parameters["Face ID"][j]
#         # PicParams["Object"] = f"""{Names_Components[i]}:{Names[i][j]+str(j)}"""
#         for k in range(len(Names[0])):
#             PicParams["Object"] = f"""{Names_Components[i]}:{Names[i][j]+str(k)}"""
#             Parameters["Picked Component Name"].append(PicParams["Object"])
#             PickFace = obj.Pick(PicParams)

# obj.WaveguidePortWithPins(Parameters)
      




#####################################################################################
# PAD IHP 
#####################################################################################


# Parameters = {}

# PAD_Width = 100
# SIG_PAD_Width = 80
# PAD_Length = 100
# PAD_Thickness = 5
# Name = ["GND_L", "Sig_L", "GND_Mid", "Sig_R", "GND_R"]
# Component_Name = ["_input", "_output"]
# Names_Chips = ["Chip_L", "Chip_R"]
# Dist = np.arange(-400, 500, 200)
# PAD_Dist = [-319, 319]
# Material = "Al"

# obj.add_Al()
# obj.add_Au()
# obj.add_Glue()
# obj.add_Si()
# obj.add_SiO2()

# # pad_spacing = 100
# # _dist = 2*pad_spacing + PAD_Width/2 + 2*PAD_Width
# # Dist = np.arange(-_dist, _dist + pad_spacing*2, pad_spacing*2)
        
# for j in range(len(Component_Name)):
#     for i in range(len(Name)):
#         if Name[i].split("_")[0] == "GND":
#             # Create squere Electrodes 
#             Parameters["Brick Lenght Max"] = PAD_Dist[j] + PAD_Length
#             Parameters["Brick Lenght Min"] = PAD_Dist[j] - PAD_Length
#             Parameters["Brick Width Max"] = Dist[i] + PAD_Width
#             Parameters["Brick Width Min"] = Dist[i] - PAD_Width
#             Parameters["Brick Hight Max"] = PAD_Thickness * 2
#             Parameters["Brick Hight Min"] = 0 
#             Parameters["Brick Name"] = Name[i] + Component_Name[j]
#             Parameters["Component Name"] = Name[i] 
#             Parameters["Material"] = "Al"
#             obj.Brick(Parameters)
#         else:
#             # Create squere Electrodes 
#             Parameters["Brick Lenght Max"] = PAD_Dist[j] + PAD_Length
#             Parameters["Brick Lenght Min"] = PAD_Dist[j] - PAD_Length
#             Parameters["Brick Width Max"] = Dist[i] + SIG_PAD_Width
#             Parameters["Brick Width Min"] = Dist[i] - SIG_PAD_Width
#             Parameters["Brick Hight Max"] = PAD_Thickness * 2
#             Parameters["Brick Hight Min"] = 0 
#             Parameters["Brick Name"] = Name[i] + Component_Name[j]
#             Parameters["Component Name"] = Name[i] 
#             Parameters["Material"] = "Al"
#             obj.Brick(Parameters)


# # for j in range(len(Component_Name)):
# #     for i in range(len(Name)):
# #         # Create squere Electrodes 
# #         Parameters["Brick Lenght Max"] = PAD_Dist[j] + PAD_Length
# #         Parameters["Brick Lenght Min"] = PAD_Dist[j] - PAD_Length
# #         Parameters["Brick Width Max"] = Dist[i] + PAD_Width
# #         Parameters["Brick Width Min"] = Dist[i] -PAD_Width
# #         Parameters["Brick Hight Max"] = PAD_Thickness * 2
# #         Parameters["Brick Hight Min"] = 0 
# #         Parameters["Brick Name"] = Name[i] + Component_Name[j]
# #         Parameters["Component Name"] = Name[i] 
# #         Parameters["Material"] = Material
# #         obj.Brick(Parameters)




# # Create two Chips Left and Right 
# for i in range(len(Names_Chips)):
#     # Create SiO2 Layer
#     # Parameters["Brick Lenght Max"] = PAD_Dist[1] + PAD_Length
#     # Parameters["Brick Lenght Min"] = PAD_Dist[0] - PAD_Length
#     if PAD_Dist[i] > 0:
#         Parameters["Brick Lenght Max"] = PAD_Dist[i] + PAD_Length 
#         Parameters["Brick Lenght Min"] = PAD_Dist[i] - PAD_Length - 50
#     else:
#         Parameters["Brick Lenght Max"] = PAD_Dist[i] - PAD_Length 
#         Parameters["Brick Lenght Min"] = PAD_Dist[i] + PAD_Length + 50
#     Parameters["Brick Width Max"] = max(Dist) + PAD_Width*3
#     Parameters["Brick Width Min"] = min(Dist) - PAD_Width*3
#     Parameters["Brick Hight Max"] = 0
#     Parameters["Brick Hight Min"] = -19.92
#     Parameters["Brick Name"] = "SiO2_Layer" + Component_Name[i]
#     Parameters["Component Name"] = Names_Chips[i]
#     Parameters["Material"] = "SiO2"
#     obj.Brick(Parameters)

#     # Create Substrate
#     # Parameters["Brick Lenght Max"] = PAD_Dist[1] + PAD_Length
#     # Parameters["Brick Lenght Min"] = PAD_Dist[0] - PAD_Length
#     if PAD_Dist[i] > 0:
#         Parameters["Brick Lenght Max"] = PAD_Dist[i] + PAD_Length 
#         Parameters["Brick Lenght Min"] = PAD_Dist[i] - PAD_Length - 50
#     else:
#         Parameters["Brick Lenght Max"] = PAD_Dist[i] - PAD_Length 
#         Parameters["Brick Lenght Min"] = PAD_Dist[i] + PAD_Length + 50
#     Parameters["Brick Width Max"] = max(Dist) + PAD_Width*3
#     Parameters["Brick Width Min"] = min(Dist) - PAD_Width*3
#     Parameters["Brick Hight Max"] = -19.92
#     Parameters["Brick Hight Min"] = -298.86
#     Parameters["Brick Name"] = "Substrate_Chip" + Component_Name[i]
#     Parameters["Component Name"] = Names_Chips[i]
#     Parameters["Material"] = "Si"
#     obj.Brick(Parameters)




# # Create Bomd Wires
# Names_Wires = ["Wire_L_GND_1", "Wire_L_GND_2", "Wire_Mid_GND_1", "Wire_Mid_GND_2", "Wire_R_GND_1", "Wire_R_GND_2"]
# # Dist2 = [-400, 0, 400]
# Dist2 = [-225, -175, -25, 25, 175, 225]
# Bond_Wire_Height = 20
# # Bond_Wire_Radius = 17.5/2
# Bond_Wire_Radius = 25/2
# Offset = [-25, 25]
# Offset_name = ["_1", "_2"]

# for k in range(len(Dist2)):

#     Parameters = {}
#     Parameters['X1'] = PAD_Dist[0]/2 
#     Parameters['Y1'] = Dist2[k]
#     Parameters['Z1'] = PAD_Thickness
#     Parameters['X2'] = PAD_Dist[1]/2 
#     Parameters['Y2'] = Dist2[k]
#     Parameters['Z2'] = PAD_Thickness

#     Points = {}
#     x = []
#     y = []
#     for i in range(0, 100):
#         x.append(i)
#         y.append(i*4)
#     x = np.array(x)
#     y = np.array(y)
#     Points['X'] = x
#     Points['Y'] = y


#     Points['X'] = CosinusCurve[:,0]
#     Points['Y'] = CosinusCurve[:,1]



#     obj.BondWire(NameWire = Names_Wires[k], Coordinates = Parameters, Height = Bond_Wire_Height, Radius = Bond_Wire_Radius , BondwireType = "Spline", Termination= "rounded", Material = "Al",  NameFolder = Names_Wires[k] + "_BondWire")
#     obj.ToSolid(SolidName = Names_Wires[k], CurveName = Names_Wires[k], NameFolder = Names_Wires[k] + "_BondWire", Material = "Al")





# # Create Bomd Wires
# Names_Wires = ["Wire_L_Sig", "Wire_R_Sig"]
# Dist2 = [-100, 100]
# Bond_Wire_Height = 20
# # Bond_Wire_Radius = 17.5/2
# Bond_Wire_Radius = 25/2

# for k in range(len(Dist2)):
#     Parameters = {}
#     Parameters['X1'] = PAD_Dist[0]/2 
#     Parameters['Y1'] = Dist2[k]
#     Parameters['Z1'] = PAD_Thickness
#     Parameters['X2'] = PAD_Dist[1]/2 
#     Parameters['Y2'] = Dist2[k]
#     Parameters['Z2'] = PAD_Thickness

#     Points = {}
#     x = []
#     y = []
#     for i in range(0, 100):
#         x.append(i)
#         y.append(i*4)
#     x = np.array(x)
#     y = np.array(y)
#     Points['X'] = x
#     Points['Y'] = y


#     Points['X'] = CosinusCurve[:,0]
#     Points['Y'] = CosinusCurve[:,1]



#     obj.BondWire(NameWire = Names_Wires[k] ,Coordinates = Parameters, Height = Bond_Wire_Height, Radius = Bond_Wire_Radius , BondwireType = "Spline", Termination= "rounded", Material = "Al",  NameFolder = Names_Wires[k] + "_BondWire")
#     obj.ToSolid(SolidName = Names_Wires[k], CurveName = Names_Wires[k], NameFolder = Names_Wires[k] + "_BondWire", Material = "Al")




# # Create Top Plate for Bond Wires and cladding
# Parameters = {}

# PAD_Width = 80
# PAD_Length = 80
# PAD_Thickness = 2.8
# Name_Shield = "Floating_Shield"
# Component_Name = "Floating_Shield"
# PAD_Dist_to_Wires = 50
# Material = "Au"
# Material_Clad = "DAF_Glue"


# # Create squere floating shield
# Parameters["Brick Lenght Max"] = PAD_Dist[0]/2 - PAD_Length/2
# Parameters["Brick Lenght Min"] = PAD_Dist[1]/2 + PAD_Length/2
# Parameters["Brick Width Max"] = max(Dist) + PAD_Width*2
# Parameters["Brick Width Min"] = min(Dist) - PAD_Width*2
# Parameters["Brick Hight Max"] = PAD_Thickness*2 + 100 + PAD_Dist_to_Wires*2
# Parameters["Brick Hight Min"] = PAD_Thickness*2 + 100
# # Parameters["Brick Hight Max"] = PAD_Thickness*2 + Bond_Wire_Height*2 + Bond_Wire_Radius*2 + PAD_Dist_to_Wires*2
# # Parameters["Brick Hight Min"] = PAD_Thickness*2 + Bond_Wire_Height*2 + Bond_Wire_Radius*2 + PAD_Dist_to_Wires
# Parameters["Brick Name"] = Name_Shield
# Parameters["Component Name"] = Component_Name
# Parameters["Material"] = Material
# obj.Brick(Parameters)



# # Create cladding
# Parameters["Brick Lenght Max"] = PAD_Dist[0]/2 - PAD_Length/2
# Parameters["Brick Lenght Min"] = PAD_Dist[1]/2 + PAD_Length/2
# Parameters["Brick Width Max"] = max(Dist) + PAD_Width*2
# Parameters["Brick Width Min"] = min(Dist) - PAD_Width*2
# # Parameters["Brick Hight Max"] = PAD_Thickness*2 + Bond_Wire_Height*2 + Bond_Wire_Radius*2 + PAD_Dist_to_Wires
# Parameters["Brick Hight Max"] = PAD_Thickness*2 + 100
# Parameters["Brick Hight Min"] = PAD_Thickness*2 
# Parameters["Brick Name"] = "Floating_Shield_Clad"
# Parameters["Component Name"] = Component_Name
# Parameters["Material"] = Material_Clad
# obj.Brick(Parameters)
# 
# 
# 
# 
# 
# 
# # Pick Faces for Input Waveguide Port
# Parameters = {}
# Facer_input = []
# Facer_output = []
# Parameters["Option"] = "Face"

# for i in range(len(Name)):
#     Facer_input.append(Name[i] + ":" + Name[i] + "_input")

# # pick input faces
# for i in range(len(Name)):
#     Parameters["Face Number"] = 4
#     Parameters["Object"] = Facer_input[i]
#     obj.Pick(Parameters)

# # Create Port Input
# Parameters = {}
# Parameters["Port Number"] = 1
# Parameters["Coordinates"] = "Picks"
# Parameters["Orientation"] = "Positive"
# Parameters["Span"] = [3, PAD_Width/2, PAD_Width/2]
# # Parameters["Picked Port Number"] = 1
# Parameters["Picked Port Number"] = [1,1,1,2,2,2]
# Parameters["Number of picks"] = 3
# Parameters["Picked Port Polarity"] = ["negative", "positive", "negative", "negative", "positive","negative"]
# Parameters["Picked Component Name"] = [ "GND_L:GND_L_input", "Sig_L:Sig_L_input", "GND_Mid:GND_Mid_input", "GND_Mid:GND_Mid_input", "Sig_R:Sig_R_input", "GND_R:GND_R_input"]
# Parameters["Face Number"] = [4,4,4,4,4,4]

# obj.WaveguidePortWithPins(Parameters)






# # Pick output faces
# Parameters = {}
# Facer_input = []
# Facer_output = []
# Parameters["Option"] = "Face"
# for i in range(len(Name)):
#     Facer_output.append(Name[i] + ":" + Name[i] + "_output")

# for i in range(len(Name)):
#     Parameters["Face Number"] = 6
#     Parameters["Object"] = Facer_output[i]
#     obj.Pick(Parameters)

# #Create Port Output
# Parameters = {}
# Parameters["Port Number"] = 2
# Parameters["Coordinates"] = "Picks"
# Parameters["Orientation"] = "Positive"
# Parameters["Span"] = [3, PAD_Width/2, PAD_Width/2]
# # Parameters["Picked Port Number"] = 1
# Parameters["Picked Port Number"] = [1,1,1,2,2,2]
# Parameters["Number of picks"] = 5
# Parameters["Picked Port Polarity"] = ["negative", "positive", "negative", "negative", "positive","negative"]
# Parameters["Picked Component Name"] = [ "GND_L:GND_L_output", "Sig_L:Sig_L_output", "GND_Mid:GND_Mid_output", "GND_Mid:GND_Mid_output", "Sig_R:Sig_R_output", "GND_R:GND_R_output"]
# Parameters["Face Number"] = [6,6,6,6,6,6]

# obj.WaveguidePortWithPins(Parameters)






Parameters = {}
Parameters["PAD Width GND"] = 100
Parameters["PAD Width Signal"] = 80
Parameters["PAD Length"] = 100
Parameters["PAD Thickness"] = 2.8


Parameters["PADs Distance"] = 219
Parameters["Bonwire height"] = 60
Parameters["Bonwire radius"] = 25/2
Parameters["Glue Thickness"] = 50
Parameters["Floating Shield Thickness"]= 35

Parameters["Port Y Span"] = 70
Parameters["Port Z Span"] = 200
Parameters["Accuracy"] = 40


obj.GSGSG_Bondwire_ChipToChip_connection(Parameters)
# 
# Parameters["Probes"] = True
# obj.GSG_Bondwire_ChipToChip_connection(Parameters)






# # Insert Monitors
# Field_Freq = np.arange(0,220,20)


# Parameters = {}
# Parameters["Monitor Type"] = "Efield"
# Parameters["Domain"] = "Frequency"


# for i in range(len(Field_Freq)):
#     Parameters["Monitor Frequency"] = Field_Freq[i]
#     obj.CreateEfieldMonitor(Parameters)






# # create_GGB_Probe()
# Parameters = {}
# Parameters["Component Name"] = "Probe_Left"
# Parameters["Orientation Angle"] = -30
# Parameters["Name"] = "GGB_L"
# obj.GGB_Probe(Parameters)



# # Move Probes Left
# Parameters = {}
# Parameters["Translate Type"] = "Translate"
# Parameters["Name Object"] = "Probe_Left"
# Parameters["Position X"] = -320
# Parameters["Position Y"] = 0
# Parameters["Position Z"] = 5


# obj.Translation(Parameters)




# # create_GGB_Probe()
# Parameters = {}
# Parameters["Component Name"] = "Probe_Right"
# Parameters["Orientation Angle"] = 30
# obj.GGB_Probe(Parameters)


# # Move Probes Right
# Parameters = {}
# Parameters["Translate Type"] = "Translate"
# Parameters["Name Object"] = "Probe_Right"
# Parameters["Position X"] = 320
# Parameters["Position Y"] = 0
# Parameters["Position Z"] = 5


# obj.Translation(Parameters)





# # Pick Faces for Input Waveguide Port
# Parameters = {}
# Parameters["Option"] = "Face"
# Parameters["Face Number"] = [3,9]
# Names = ["Probe_Left:Signal_Probe", "Probe_Right:Signal_Probe"]


# # Port Parmaeters
# Parameters["Coordinates"] = "Picks"
# Orientation = ["positive", "positive"]
# Parameters["Span"] = [3, PAD_Width/2, PAD_Width/2]
# Parameters["Picked Port Number"] = 1
# Parameters["Number of picks"] = 1
# Parameters["Picked Port Polarity"] = ["positive", "negative"]
# Componentname = ["Probe_Left", "Probe_Right"]
# PartName = ["Signal_Probe", "Body_Probe" ]
# Parameters["Number of picks"] = 5
# FaceNum =  [ 3, 9 ]
# Parameters["Picked Port Number"] 




# for i in range(len(Names)): 
#     Parameters["Orientation"] = Orientation[i]

#     for j in range(len(FaceNum )): 
#         Parameters["Object"] = f"{Componentname[i]}:{PartName[j]}"
       
        
#         Parameters["Face Number"] =  FaceNum[j]
        
#         obj.Pick(Parameters)
    
    
#     # Create Port Input
#     Parameters["Picked Component Name"] = Names[i]
#     Parameters["Port Number"] = i+1
#     Parameters["Picked Port Number"] = 1
#     Parameters["Picked Component Name"] = []
#     Parameters["Face Number"] = [3,9]
#     for k in range(len(FaceNum)):
#         Parameters["Picked Component Name"].append(f"{Componentname[i]}:{PartName[k]}")
#     obj.WaveguidePortWithPins(Parameters)

    


    

  

