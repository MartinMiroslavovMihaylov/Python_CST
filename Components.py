
import numpy as np 
import sys
import os
# Add the directory containing the project to sys.path
current_path = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_path)
import Functions as VBA





# MZM Design
def MZM(Parameters, CST):
    """Create an MZM Modulator. Materials used:
                                            Gold - For the electrodes
                                            LiNbO3 - For Optical Waveguides
                                            SiO2 - For Substrate
    Args:
        Parameters (dict): Dictionary with all the needed values
                Parameters["Lenght_Electrodes"] : Length of the Electrodes. The Waveguides will be 2 (Units) longer then the electrodes.
                Parameters["Width GND"] : Width of the GND electrodes
                Parameters["Width Signal"] : Width of the Signal Electrode
                Parameters["Width WG"] : Top wWidth of the optical waveguide. It is an Rib waveguide
                Parameters["Gap"] : Gap between Signal and optical Waveguide
                Parameters["angle"] : Angle of the side wall of the optical waveguide
                Parameters["High Electrodes"] : Hight of the Electodes
                Parameters["High WG"] : Hight of the optical Waveguide
                Parameters["High Slab"] : Hight of the Slab. When choosen "0" no Slab will be implemented. 
                Parameters["High Substrate"] : Hight of the substrate

        CST (Object): The CST Obejct that you use to load your project. 
    """

    Length_MZM = Parameters["Electrodes Lenght"]
    Length_WG = Parameters["Electrodes Lenght"] + 2
    Width_Electrodes = Parameters["Width GND"]
    Width_Signal = Parameters["Width Signal"]
    Width_WG = Parameters["Width WG"]
    Gap = Parameters["Gap"]
    Angle = Parameters["angle"]
    Height_Electrodes = Parameters["High Electrodes"]
    Height_WG = Parameters["High WG"]
    Heigh_Slab = Parameters["High Slab"]
    Height_Substrate = Parameters["High Substrate"]


    # Set optical Material Properties Eps X,Y,Z
    Data = {}
    Data["X"] = 4.906
    Data["Y"] = 4.584095
    Data["Z"] = 4.906

    MaterialDAta = VBA.Material("LiNbO3", Data)
    CST.schematic.execute_vba_code(MaterialDAta, timeout=None)

    # MaterialDAta = VBA.Material_Silicon("Silicon (lossy)")
    # CST.schematic.execute_vba_code(MaterialDAta, timeout=None)

    MaterialDAta = VBA.Material_SiO2("SiO2")
    CST.schematic.execute_vba_code(MaterialDAta, timeout=None)

    MaterialDAta = VBA.Material_Au("Au")
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
    TestBrick = Brick('LNOI_Substrate', Parameters, Material = "SiO2")
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
    WG = VBA.Poligon_3D(WGName = "Waveguide_Left", Points = PointsLeft)
    CST.schematic.execute_vba_code(WG, timeout=None)
    RibWG_Test = VBA.RibWaveguide_ToSolid("Waveguide_Left", WaveguideName = "Waveguide_Left", WG_Hight = Height_WG, Angle = -Angle, WGFolderName = "Waveguide_Left", WGName = "Waveguide_Left", Material="LiNbO3")
    CST.schematic.execute_vba_code(RibWG_Test, timeout=None)



    PointsRight = {}
    PointsRight["X"] = [-Length_WG/2, -Length_WG/2, Length_WG/2, Length_WG/2, -Length_WG/2]
    PointsRight["Y"] = [PosRight[0], PosRight[1], PosRight[1], PosRight[0], PosRight[0]]
    PointsRight["Z"] = [HeightZ, HeightZ, HeightZ, HeightZ, HeightZ]

    # Waveguide and Waveguide to solid
    WG = VBA.Poligon_3D(WGName = "Waveguide_Right", Points = PointsRight)
    CST.schematic.execute_vba_code(WG, timeout=None)
    RibWG_Test = VBA.RibWaveguide_ToSolid("Waveguide_Right", WaveguideName = "Waveguide_Right", WG_Hight = -Height_WG, Angle = -Angle, WGFolderName = "Waveguide_Right", WGName = "Waveguide_Right", Material="LiNbO3")
    CST.schematic.execute_vba_code(RibWG_Test, timeout=None)




# MZM Design
def PhaseModulator(Parameters, CST):
    """Create an Phase Modulator. Materials used:
                                            Gold - For the electrodes
                                            LiNbO3 - For Optical Waveguides
                                            SiO2 - For Substrate
    Args:
        Parameters (dict): Dictionary with all the needed values
                Parameters["Lenght_Electrodes"] : Length of the Electrodes. The Waveguides will be 2 (Units) longer then the electrodes.
                Parameters["Width GND"] : Width of the GND electrodes
                Parameters["Width Signal"] : Width of the Signal Electrode
                Parameters["Width WG"] : Top wWidth of the optical waveguide. It is an Rib waveguide
                Parameters["Gap"] : Gap between Signal and optical Waveguide
                Parameters["angle"] : Angle of the side wall of the optical waveguide
                Parameters["High Electrodes"] : Hight of the Electodes
                Parameters["High WG"] : Hight of the optical Waveguide
                Parameters["High Slab"] : Hight of the Slab. When choosen "0" no Slab will be implemented. 
                Parameters["High Substrate"] : Hight of the substrate

        CST (Object): The CST Obejct that you use to load your project. 
    """

    Length_MZM = Parameters["Electrodes Lenght"]
    Length_WG = Parameters["Electrodes Lenght"] + 2
    Width_Electrodes = Parameters["Width GND"]
    Width_Signal = Parameters["Width Signal"]
    Width_WG = Parameters["Width WG"]
    Gap = Parameters["Gap"]
    Angle = Parameters["angle"]
    Height_Electrodes = Parameters["High Electrodes"]
    Height_WG = Parameters["High WG"]
    Heigh_Slab = Parameters["High Slab"]
    Height_Substrate = Parameters["High Substrate"]


    # Set optical Material Properties Eps X,Y,Z
    Data = {}
    Data["X"] = 4.906
    Data["Y"] = 4.584095
    Data["Z"] = 4.906

    MaterialDAta = VBA.Material("LiNbO3", Data)
    CST.schematic.execute_vba_code(MaterialDAta, timeout=None)

    # MaterialDAta = VBA.Material_Silicon("Silicon (lossy)")
    # CST.schematic.execute_vba_code(MaterialDAta, timeout=None)

    MaterialDAta = VBA.Material_SiO2("SiO2")
    CST.schematic.execute_vba_code(MaterialDAta, timeout=None)

    MaterialDAta = VBA.Material_Au("Au")
    CST.schematic.execute_vba_code(MaterialDAta, timeout=None)



    # Calculate Top and Bottom Widths of the Waveguide
    x = abs(Height_WG / (np.cos((Angle) * np.pi / 180)))  # in Radians
    extention = np.sqrt(x ** 2 - Height_WG ** 2)
    Width_WG_New = Width_WG + 2 * extention 

    # Global Parameters
    WidthObject = 2*Width_WG_New + Width_Signal + Width_Electrodes + 2*Gap


    # Substrate Definition 
    Parameters = {}
    Parameters['X1'] = Length_WG/2
    Parameters['X2'] = -Length_WG/2
    Parameters['Y1'] = (WidthObject+2)/2
    Parameters['Y2'] = -(WidthObject+2)/2 
    Parameters['Z1'] = Height_Substrate/2
    Parameters['Z2'] = -Height_Substrate/2
    TestBrick = Brick('LNOI_Substrate', Parameters, Material = "SiO2")
    CST.schematic.execute_vba_code(TestBrick, timeout=None)



    # Slab definition
    if Heigh_Slab == 0:
        HeightZ = Height_Substrate/2
    else:
        Parameters = {}
        Parameters['X1'] = Length_WG/2
        Parameters['X2'] = -Length_WG/2
        Parameters['Y1'] = (WidthObject+2)/2
        Parameters['Y2'] = -(WidthObject+2)/2 
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
    NamesElectrodes = "Electrode", "Signal"

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
    TestBrick = Brick(NamesElectrodes[1], Parameters, Material = "Au")
    CST.schematic.execute_vba_code(TestBrick, timeout=None)




    # Plase Waveguides
    GND_Left_Corner = (-WidthObject/2 + Width_Electrodes) + Gap + Width_WG_New/2
    
    PosLeft = [round((GND_Left_Corner  + Width_WG_New/2),2), round((GND_Left_Corner  - Width_WG_New/2),2)] 
  

    PointsLeft = {}
    PointsLeft["X"] = [-Length_WG/2, -Length_WG/2, Length_WG/2, Length_WG/2, -Length_WG/2]
    PointsLeft["Y"] = [PosLeft[0], PosLeft[1], PosLeft[1], PosLeft[0], PosLeft[0]]
    PointsLeft["Z"] = [HeightZ, HeightZ, HeightZ, HeightZ, HeightZ]

    # Waveguide and Waveguide to solid
    WG = VBA.Poligon_3D(WGName = "Waveguide_Left", Points = PointsLeft)
    CST.schematic.execute_vba_code(WG, timeout=None)
    RibWG_Test = VBA.RibWaveguide_ToSolid("Waveguide_Left", WaveguideName = "Waveguide_Left", WG_Hight = Height_WG, Angle = -Angle, WGFolderName = "Waveguide_Left", WGName = "Waveguide_Left", Material="LiNbO3")
    CST.schematic.execute_vba_code(RibWG_Test, timeout=None)






# # Ports and Solvers and Background
def Squere_Waveguide(Parameters, CST):
    """This function generate and simple straight waveguide. Materials used:
                                            Gold - For the electrodes
                                            LiNbO3 - For Optical Waveguides
                                            SiO2 - For Substrate

    Args:
        Parameters (dict): Dictionary with all the needed values
                Parameters["Lenght WG"] : Length of the Waveguide.
                Parameters["High_WG"] : Hight of the optical Waveguide
                Parameters["Width WG"] : : Top Width of the optical waveguide. It is an Rib waveguide.
                Parameters["Substrate Height"] : Hight of the substrate.
                Parameters["Slab Heigh"] : Hight of the Slab. When choosen "0" no Slab will be implemented. 
      
    """
    Length_WG = Parameters["Lenght WG"]
    Hight_WG = Parameters["Hight WG"]
    Width_WG = Parameters["Width WG"]
    Height_Substrate = Parameters["Substrate Height"]
    Heigh_Slab = Parameters["Slab Heigh"]
    # Material = Parameters["Materials"]



    # Set optical Material Properties Eps X,Y,Z
    Data = {}
    Data["X"] = 4.906
    Data["Y"] = 4.584095
    Data["Z"] = 4.906

    MaterialDAta = VBA.Material("LiNbO3", Data)
    CST.schematic.execute_vba_code(MaterialDAta, timeout=None)

    # MaterialDAta = VBA.Material_Silicon("Silicon (lossy)")
    # CST.schematic.execute_vba_code(MaterialDAta, timeout=None)

    
    MaterialDAta = VBA.Material_SiO2("SiO2")
    CST.schematic.execute_vba_code(MaterialDAta, timeout=None)

    

    # Global Parameters for the WG
    WidthObject = 4*Width_WG 
 

    # Substrate Definition 
    Parameters = {}
    Parameters['X1'] = -Length_WG/2
    Parameters['X2'] = Length_WG/2
    Parameters['Y1'] = -WidthObject/2
    Parameters['Y2'] = WidthObject/2 
    Parameters['Z1'] = -Height_Substrate/2
    Parameters['Z2'] = Height_Substrate/2
    TestBrick = Brick('LNOI_Substrate', Parameters, Material = "SiO2")
    CST.schematic.execute_vba_code(TestBrick, timeout=None)


    # Slab definition
    if Heigh_Slab == 0:
        HeightZ = Height_Substrate/2
    else:
        Parameters = {}
        Parameters['X1'] = -Length_WG/2
        Parameters['X2'] = Length_WG/2
        Parameters['Y1'] = -WidthObject/2
        Parameters['Y2'] = WidthObject/2 
        Parameters['Z1'] = Height_Substrate/2
        Parameters['Z2'] = Height_Substrate/2 + Heigh_Slab
        TestBrick = Brick('LNOI_Slab', Parameters, Material = "LiNbO3")
        CST.schematic.execute_vba_code(TestBrick, timeout=None)
        HeightZ = Height_Substrate/2 + Heigh_Slab
    
    # Create squere Waveguide
    Parameters = {}
    Parameters['X1'] = -Length_WG/2
    Parameters['X2'] = Length_WG/2
    Parameters['Y1'] = -Width_WG/2 
    Parameters['Y2'] = Width_WG/2  
    Parameters['Z1'] = HeightZ
    Parameters['Z2'] = HeightZ + Hight_WG
    TestBrick = Brick("WG", Parameters, Material = "LiNbO3")
    CST.schematic.execute_vba_code(TestBrick, timeout=None)



    



def BondWire(NameWire, Coordinates, Height, Radius, BondwireType = "Spline", CenterPosition = 0.5, alpha = None, beta = None, Material = None, SolidWireModel = True, Termination = None, NameFolder = None):
    """Create Bond Wire

    Args:
        NameWire (str): Name of the Bondwire
        Coordinates (dict): Dictionary with Coordinates in X,Y,Z plane to create the Bondwire:  
                            For Example 
                                Parameters = {}
                                Parameters['X1'] = 0
                                Parameters['Y1'] = 0
                                Parameters['Z1'] = 0
                                Parameters['X2'] = 5
                                Parameters['Y2'] = 5
                                Parameters['Z2'] = 0
        Height (int/float): Hight of the middle point of the Bondwire
        Radius (int/float): Radius of the bond wire. Bond wire is an cylinder type of object.
        BondwireType (str, optional): The type of bond wire. Defaults to "Spline".
        CenterPosition (float, optional): The center Position of the Height. This can be moved to make 
                                           an object that dont have the top height in the middle. Defaults to 0.5.
        alpha (_type_, optional): _description_. Defaults to None.
        beta (_type_, optional): _description_. Defaults to None.
        Material (str, optional): Material for of the Bond wire. For now you need to load the material 
                                  in your simulation and then use this function. Otherwise the material 
                                  will not be found. Defaults to "PCE".
        SolidWireModel (bool, optional): _description_. Defaults to True.
        Termination (_type_, optional): _description_. Defaults to None.
        NameFolder (str, optional): The name of the folder. Defaults to name of the wire.

    Raises:
        ValueError: Error massage
        ValueError: Error massage
        ValueError: Error massage

    Returns:
        str: String with VBA Code 
    """
    # Check the BondwireType
    if BondwireType in ["Spline", "JEDEC4", "JEDEC5"]:
        BondwireType = BondwireType
        if BondwireType == "JEDEC5":
            if alpha == None or beta == None:
                raise ValueError(
                    'When using "JEDEC5" an alpha and betta parameters need to be defined. See CST documentation!! ')
            else:
                pass
    else:
        raise ValueError("BondwireType can be on of ['Spline', 'JEDEC4', 'JEDEC5']")


    # Material Def
    if Material == None:
        Material = "PCE"
    else:
        Materila = Material


    # Check Termination
    if Termination != None:
        if Termination in ["natural" , "rounded" , "extended"]:
            Termination = Termination
        else:
            raise ValueError('Termination can be only one of ["natural" , "rounded" , "extended"]')
    else:
        Termination = "natural"

    # Folder Name
    if NameFolder == None:
        NameFolder = NameWire
    else:
        NameFolder = NameFolder


    if BondwireType in ["Spline", "JEDEC4"]:
        component = 'Sub Main () ' \
                        '\nWith Wire' + \
                            '\n.Reset' + \
                            '\n.Folder ' + '"' + NameFolder + '"' + \
                            '\n.Type "Bondwire"' + \
                            '\n.Name ' + '"' + NameWire + '"' + \
                            '\n.BondWireType ' + '"' + BondwireType + '"' + \
                            '\n.Point1 ' + '"'+ str(Coordinates["X1"]) + '"' + ',' + '"' + str(Coordinates["Y1"]) + '"' + ',' + '"' + str(Coordinates["Z1"])  + '"' + ', "False"' + \
                            '\n.Point2 ' + '"'+ str(Coordinates["X2"]) + '"' + ',' + '"' + str(Coordinates["Y2"]) + '"' + ',' + '"' + str(Coordinates["Z2"])  + '"' + ', "False"' + \
                            '\n.Height ' + '"' + str(Height) + '"' + \
                            '\n.Radius ' + '"' + str(Radius) + '"' + \
                            '\n.RelativeCenterPosition ' + '"' + str(CenterPosition) + '"' + \
                            '\n.Material ' + '"' + Material + '"' + \
                            '\n.SolidWireModel ' + '"' + str(SolidWireModel) + '"' + \
                            '\n.Termination ' + '"' + Termination + '"' +\
                            '\n.Add' + \
                        '\nEnd With' + \
                    '\nEnd Sub'
        return component
    # '\n.Folder ' + '"' + NameFolder + '"' + \

    else:
        component = 'Sub Main () ' \
                        '\nWith Wire' + \
                            '\n.Reset' + \
                            '\n.Folder ' + '"' + NameFolder + '"' + \
                            '\n.Type "Bondwire"' + \
                            '\n.Name ' + '"' + NameWire + '"' + \
                            '\n.BondWireType ' + '"' + BondwireType + '"' + \
                            '\n.Point1 ' + '"' + str(Coordinates["X1"]) + '"' + ',' + '"' + str(Coordinates["Y1"]) + '"' + ',' + '"' + str(Coordinates["Z1"]) + '"' + ', "False"' + \
                            '\n.Point2 ' + '"' + str(Coordinates["X2"]) + '"' + ',' + '"' + str(Coordinates["Y2"]) + '"' + ',' + '"' + str(Coordinates["Z2"]) + '"' + ', "False"' + \
                            '\n.Height ' + '"' + str(Height) + '"' + \
                            '\n.Radius ' + '"' + str(Radius) + '"' + \
                            '\n.alpha ' + '"' + str(alpha) + '"' + \
                            '\n.beta ' + '"' + str(beta) + '"' + \
                            '\n.Material ' + '"' + Material + '"' + \
                            '\n.SolidWireModel ' + '"' + str(SolidWireModel) + '"' + \
                            '\n.Termination ' + '"' + Termination + '"' + \
                            '\n.Add' + \
                        '\nEnd With' + \
                    '\nEnd Sub'



        return component
    




def CurveWire(NameWire, Radius, Points = None, Material = None, SolidWireModel = True, Termination = None, NameFolder = None, CurveFolderName = None, CurveName = None):
    """_summary_

    Args:
        NameWire (str): Name of the wire
        Radius (int/float): Radius of the curve
        Points (dict, optional): Dictionary of X and  Y points for the curve. Defaults to None.
        Material (str, optional): Material for of the Bond wire. For now you need to load the material 
                                  in your simulation and then use this function. Otherwise the material 
                                  will not be found. Defaults to "PCE".
        SolidWireModel (bool, optional): _description_. Defaults to True.
        Termination (_type_, optional): _description_. Defaults to None.
        NameFolder (str, optional): Name of the Folder. Defaults to None.
        CurveFolderName (str, optional): Curve folder name. Defaults to None.
        CurveName (str, optional): Curve name. Defaults to None.

    Raises:
        ValueError: Error massage
        ValueError: Error massage
        ValueError: Error massage

    Returns:
        str: String with VBA Code 
    """
    # Check of Dict of Points for X and Y are given
    if Points == None:
        raise ValueError('To create an Curve you need an Dictionary with two Variable "X" and "Y". The Values for Dict["X"] can be an array of points')
    else:
        pass


    # Folder Name
    if NameFolder == None:
        NameFolder = NameWire
    else:
        NameFolder = NameFolder


    # Check Termination
    if Termination != None:
        if Termination in ["natural", "rounded", "extended"]:
            Termination = Termination
        else:
            raise ValueError('Termination can be only one of ["natural" , "rounded" , "extended"]')
    else:
        Termination = "natural"

    # Material Def
    if Material == None:
        Material = "PCE"
    else:
        Materila = Material

    # Check if Curve Name and Folder Name are given
    if CurveFolderName == None or CurveName == None:
        raise ValueError("To create an Curve Wire you need to give the Curve Folder Name (CurveFolderName) and the Curve Name (CurveName) that you already gave!")


    component = 'Sub Main () ' \
                    '\nWith Wire' + \
                        '\n.Reset' + \
                        '\n.Folder ' + '"' + NameFolder + '"' + \
                        '\n.Type "Curvewire"' + \
                        '\n.Name ' + '"' + NameWire + '"' + \
                        '\n.Curve ' + '"' + CurveFolderName + ':' + CurveName + '"' + \
                        '\n.Radius ' + '"' + str(Radius) + '"' + \
                        '\n.SolidWireModel ' + '"' + str(SolidWireModel) + '"' + \
                        '\n.Material ' + '"' + Material + '"' + \
                        '\n.Termination ' + '"' + Termination + '"' + \
                        '\n.Add' + \
                    '\nEnd With' + \
                '\nEnd Sub'
    return component




def Brick(BrickName, Coordinates, NameComponent = None, Material = None ):
    """Create an squere Brick object.

    Args:
        BrickName (str): Name of thr brick object
        Coordinates (dict): Dictionary with Min/max Values for the coordinates of the brick in the room.
                            Coordinates["X1"] - X min
                            Coordinates["X2"] - X max
                            Coordinates["Y1"] - Y min
                            Coordinates["Y2"] - Y max
                            Coordinates["Z1"] - Z min 
                            Coordinates["Z2"] - Z max
        NameComponent (str, optional): Name fo the Component. Defaults to None.
        Material (str, optional): Material for of the Bond wire. For now you need to load the material 
                                  in your simulation and then use this function. Otherwise the material 
                                  will not be found. Defaults to None.

    Returns:
        str: String with VBA Code 
    """

    # Folder Name
    if NameComponent != None:
        NameComponent = NameComponent
    else:
        NameComponent = BrickName

    # Material Def
    if Material == None:
        Material = "Copper (annealed)"
    else:
        Material = Material


    SolidWire = 'Sub Main () ' \
                    '\nWith Brick' + \
                        '\n.Reset' + \
                        '\n.Name ' + '"' + BrickName + '1"' + \
						'\n.Component ' + '"' + NameComponent + '"' + \
                        '\n.Material ' + '"' + Material + '"' + \
                        '\n.Xrange ' + '"' + str(Coordinates["X1"]) + '"' + ',' + '"' + str(Coordinates["X2"]) + '"' + \
						'\n.Yrange ' + '"' + str(Coordinates["Y1"]) + '"' + ',' + '"' + str(Coordinates["Y2"]) + '"' + \
						'\n.Zrange ' + '"' + str(Coordinates["Z1"]) + '"' + ',' + '"' + str(Coordinates["Z2"]) + '"' + \
						'\n.Create' + \
                    '\nEnd With' + \
                '\nEnd Sub'
    return SolidWire




