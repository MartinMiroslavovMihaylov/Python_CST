# Components Main Builder

import numpy as np 
import sys 
import os 
# Add the directory containing the project to sys.path
current_path = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_path)
from VBA_Test_Syntax import *





# MZM Design
def MZM(Parameters, CST):
    
    Length_MZM = Parameters["Lenght_Electrodes"]
    Length_WG = Parameters["Lenght_Electrodes"] + 10
    Width_Electrodes = Parameters["Width_GND"]
    Width_Signal = Parameters["Width_Signal"]
    Width_WG = Parameters["WG_Width"]
    Gap = Parameters["Gap"]
    Angle = Parameters["angle"]
    Height_Electrodes = Parameters["High_Electrodes"]
    Height_WG = Parameters["High_WG"]
    Heigh_Slab = Parameters["High_Slab"]
    Height_Substrate = Parameters["High_Substrate"]


    # Set optical Material Properties Eps X,Y,Z
    Data = {}
    Data["X"] = 27
    Data["Y"] = 43
    Data["Z"] = 43

    MaterialDAta = Material("LiNbO3", Data)
    CST.schematic.execute_vba_code(MaterialDAta, timeout=None)

    MaterialDAta = Material_Silicon("Silicon (lossy)")
    CST.schematic.execute_vba_code(MaterialDAta, timeout=None)

    MaterialDAta = Material_Au("Au")
    CST.schematic.execute_vba_code(MaterialDAta, timeout=None)



    # Calculate Top and Bottom Widths of the Waveguide
    x = abs(Height_WG / (np.cos((Angle) * np.pi / 180)))  # in Radians
    extention = np.sqrt(x ** 2 - Height_WG ** 2)
    Width_WG_New = Width_WG + 2 * extention 

    # Global Parameters
    WidthObject = 2*Width_WG_New + Width_Signal + 2*Width_Electrodes + 4*Gap


    # Substrate Definition 
    Parameters = {}
    Parameters['X1'] = Length_WG/2
    Parameters['X2'] = -Length_WG/2
    Parameters['Y1'] = WidthObject/2
    Parameters['Y2'] = -WidthObject/2 
    Parameters['Z1'] = Height_Substrate/2
    Parameters['Z2'] = -Height_Substrate/2
    TestBrick = Brick('LNOI_Substrate', Parameters, Material = "Silicon (lossy)")
    CST.schematic.execute_vba_code(TestBrick, timeout=None)



    # Slab definition
    if Heigh_Slab == 0:
        HeightZ = Height_Substrate/2
    else:
        Parameters = {}
        Parameters['X1'] = Length_WG/2
        Parameters['X2'] = -Length_WG/2
        Parameters['Y1'] = WidthObject/2
        Parameters['Y2'] = -WidthObject/2 
        Parameters['Z1'] = Height_Substrate/2
        Parameters['Z2'] = Height_Substrate/2 + Heigh_Slab
        TestBrick = Brick('LNOI_Slab', Parameters, Material = "LiNbO3")
        CST.schematic.execute_vba_code(TestBrick, timeout=None)
        HeightZ = Height_Substrate/2 + Heigh_Slab

        
    # Calculate Top and Bottom Widths of the Waveguide
    x = abs(Height_WG / (np.cos((Angle) * np.pi / 180)))  # in Radians
    extention = np.sqrt(x ** 2 - Height_WG ** 2)
    Width_WG_New = Width_WG + 2 * extention 

    #  Electrodes 
    NamesElectrodes = "Electrode_Left", "Electrode_Right", "Signal"

    # Plase Electrodes
    Parameters = {}
    Parameters['X1'] = Length_MZM/2
    Parameters['X2'] = -Length_MZM/2
    Parameters['Y1'] = -WidthObject/2
    Parameters['Y2'] = -WidthObject/2 + Width_Electrodes
    Parameters['Z1'] = HeightZ
    Parameters['Z2'] = HeightZ + Height_Electrodes
    TestBrick = Brick(NamesElectrodes[0], Parameters, Material = "Au")
    CST.schematic.execute_vba_code(TestBrick, timeout=None)


    Parameters = {}
    Parameters['X1'] = Length_MZM/2
    Parameters['X2'] = -Length_MZM/2
    Parameters['Y1'] = -WidthObject/2 + Width_Electrodes + Width_WG_New + Gap*2
    Parameters['Y2'] = -WidthObject/2 + Width_Electrodes + Width_WG_New + Gap*2 + Width_Signal
    Parameters['Z1'] = HeightZ
    Parameters['Z2'] = HeightZ + Height_Electrodes
    TestBrick = Brick(NamesElectrodes[2], Parameters, Material = "Au")
    CST.schematic.execute_vba_code(TestBrick, timeout=None)



    Parameters = {}
    Parameters['X1'] = Length_MZM/2
    Parameters['X2'] = -Length_MZM/2
    Parameters['Y1'] = WidthObject/2 
    Parameters['Y2'] = WidthObject/2 - Width_Electrodes 
    Parameters['Z1'] = HeightZ
    Parameters['Z2'] = HeightZ+Height_Electrodes
    TestBrick = Brick(NamesElectrodes[1], Parameters, Material = "Au")
    CST.schematic.execute_vba_code(TestBrick, timeout=None)


    # Plase Waveguides
    GND_Left_Corner = (-WidthObject/2 + Width_Electrodes) + Gap + Width_WG_New/2
    GND_Right_Corner = (WidthObject/2 - Width_Electrodes) - Gap - Width_WG_New/2

    PosLeft = [round((GND_Left_Corner  + Width_WG_New/2),2), round((GND_Left_Corner  - Width_WG_New/2),2)] 
    PosRight = [round((GND_Right_Corner - Width_WG_New/2),2), round((GND_Right_Corner + Width_WG_New/2),2)] 



    PointsLeft = {}
    PointsLeft["X"] = [-Length_WG/2, -Length_WG/2, Length_WG/2, Length_WG/2, -Length_WG/2]
    PointsLeft["Y"] = [PosLeft[0], PosLeft[1], PosLeft[1], PosLeft[0], PosLeft[0]]
    PointsLeft["Z"] = [HeightZ, HeightZ, HeightZ, HeightZ, HeightZ]

    # Waveguide and Waveguide to solid
    WG = RibWG_Z(WGName = "Waveguide_Left", Points = PointsLeft)
    CST.schematic.execute_vba_code(WG, timeout=None)
    RibWG_Test = RibWaveguide_ToSolid("Waveguide_Left", WaveguideName = "Waveguide_Left", WG_Hight = Height_WG, Angle = -Angle, WGFolderName = "Waveguide_Left", WGName = "Waveguide_Left", Material="LiNbO3")
    CST.schematic.execute_vba_code(RibWG_Test, timeout=None)



    PointsRight = {}
    PointsRight["X"] = [-Length_WG/2, -Length_WG/2, Length_WG/2, Length_WG/2, -Length_WG/2]
    PointsRight["Y"] = [PosRight[0], PosRight[1], PosRight[1], PosRight[0], PosRight[0]]
    PointsRight["Z"] = [HeightZ, HeightZ, HeightZ, HeightZ, HeightZ]

    # Waveguide and Waveguide to solid
    # WG = RibWG(WGName = "Waveguide_Right", Points = PointsRight)
    WG = RibWG_Z(WGName = "Waveguide_Right", Points = PointsRight)
    CST.schematic.execute_vba_code(WG, timeout=None)
    RibWG_Test = RibWaveguide_ToSolid("Waveguide_Right", WaveguideName = "Waveguide_Right", WG_Hight = -Height_WG, Angle = -Angle, WGFolderName = "Waveguide_Right", WGName = "Waveguide_Right", Material="LiNbO3")
    CST.schematic.execute_vba_code(RibWG_Test, timeout=None)





# # Ports and Solvers and Background
# def MZM_Solver():

