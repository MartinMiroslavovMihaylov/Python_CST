import numpy as np
import scipy.cluster
import scipy.constants
from Components import *



def Material(Name, Values):
    """Add Anisotropic material to Material Libs

    Args:
        Name (str): Name of the Material
        Values (Dict): Dictionary with the material Values:
                        Values["X"] : X Epsilon Value
                        Values["Y"] : Y Epsilon Value
                        Values["Z"] : Z Epsilon Value
    Returns:
        str: String with the VBA code
    """

    # Add Meterial . Values is array with X, Y and Z Values for the Anisotropic Material Permittivity
    component = 'Sub Main () ' \
                    '\nWith Material' + \
                        '\n.Reset' + \
                        '\n.Name ' + '"' + Name + '"' + \
                        '\n.FrqType "All"' + \
                        '\n.Type "Anisotropic"' + \
                        '\n.SetMaterialUnit "GHz", "um"' + \
                        '\n.MuX "1"' + \
                        '\n.MuY "1"' + \
                        '\n.MuZ "1"' + \
                        '\n.EpsilonX '+ '"' + str(Values["X"]) + '"' + \
                        '\n.EpsilonY '+ '"' + str(Values["Y"]) + '"' + \
                        '\n.EpsilonZ '+ '"' + str(Values["Z"]) + '"' + \
                        '\n.DispModelEps "None"' + \
                        '\n.AddDispEpsPole1stOrderX "0","0"' + \
                        '\n.Color "1","0", "0"' + \
                        '\n.Create' + \
                    '\nEnd With' + \
                '\nEnd Sub'

    return component    





def Material_Silicon(Name):
    """Add silicon to the Material Library

    Args:
        Name (str): Name of the Material

    Returns:
        str: String with the VBA code
    """
    # Add Meterial . Values is array with X, Y and Z Values for the Anisotropic Material Permittivity
    component = 'Sub Main () ' \
                    '\nWith Material' + \
                        '\n.Reset' + \
                        '\n.Name ' + '"' + Name + '"' + \
                        '\n.FrqType "All"' + \
                        '\n.Type "Normal"' + \
                        '\n.SetMaterialUnit "GHz", "um"' + \
                        '\n.Mu "1"' + \
                        '\n.Epsilon "11.9"' + \
                        '\n.Rho 2330' + \
                        '\n.ThermalConductivity "148"' + \
                        '\n.ThermalType "Normal"' + \
                        '\n.MechanicsType "Normal"' + \
                        '\n.SpecificHeat "700"' + \
                        '\n.YoungsModulus "112"' + \
                        '\n.PoissonsRatio "0.28"' + \
                        '\n.Color "0.5","0.5", "0.5"' + \
                        '\n.Create' + \
                    '\nEnd With' + \
                '\nEnd Sub'

    return component




def Material_Au(Name):
    """Add gold to the Material Library
    Args:
        Name (str): Name of the Material

    Returns:
        str: String with the VBA code
    """
    # Add Meterial Gold.
    component = 'Sub Main () ' \
                    '\nWith Material' + \
                        '\n.Reset' + \
                        '\n.Name ' + '"' + Name + '"' + \
                        '\n.FrqType "All"' + \
                        '\n.Type "lossy metal"' + \
                        '\n.SetMaterialUnit "GHz", "um"' + \
                        '\n.Mu "1"' + \
                        '\n.Sigma "4.561e+007"' + \
                        '\n.Rho "19320"' + \
                        '\n.ThermalConductivity "130"' + \
                        '\n.ThermalType "Normal"' + \
                        '\n.SpecificHeat "700"' + \
                        '\n.YoungsModulus "78"' + \
                        '\n.PoissonsRatio "42"' + \
                        '\n.Color "1","0.84", "0"' + \
                        '\n.Create' + \
                    '\nEnd With' + \
                '\nEnd Sub'

    return component




def Curve(CurveName, Points):

    CurveName  = CurveName
    # Extract the first points as starting Points
    Start_PointX = Points['X'][0]
    Start_PointY = Points['X'][0]
    tmp = []
    tmp.append('Sub Main () ' + '\nWith Polygon' + '\n.Reset' + '\n.Name ' + '"' + CurveName  + '"' + '\n.Curve ' + '"' + CurveName  + '"' + '\n.Point ' + '"' + str(Start_PointX) + '"' + ',' + '"' + str(Start_PointY) + '"')
    for i in range(1, len(Points['X'])):
        b = '\n.LineTo ' + '"' + str(Points['X'][i]) + '"' + ',' + '"' + str(Points['Y'][i]) + '"'
        tmp.append(b)
    tmp.append('\n.Create' + '\nEnd With' + '\nEnd Sub')
    CurveData = ''.join(tmp)

    return CurveData




def ToSolid(SolidName, CurveName = "Polygon", NameFolder = None, Material = None ):
    """This function transfer an function or polynom to solid object.

    Args:
        SolidName (str): Name of the solid object
        CurveName (str, optional): Name of the curve. Defaults to "Polygon".
        NameFolder (str, optional): Name of the folder. Defaults to Curve name.
        Material (str, optional): Material for of the Bond wire. For now you need to load the material 
                                  in your simulation and then use this function. Otherwise the material 
                                  will not be found. Defaults to Aluminum.

    Returns:
        str: String with VBA Code 
    """

    #Check Curve Name
    if CurveName != "Polygon":
        CurveName = CurveName
    else:
        pass

    # Folder Name
    if NameFolder != None:
        NameFolder = NameFolder
    else:
        NameFolder = CurveName

    # Material Def
    if Material == None:
        Material = "Aluminum"
    else:
        Material = Material



    SolidWire = 'Sub Main () ' \
                    '\nWith Wire' + \
                        '\n.Reset' + \
                        '\n.SolidName "component1:'+ SolidName + '"' + \
                        '\n.Name ' + '"' + CurveName + '"' + \
                        '\n.Folder ' + '"' + NameFolder +'"' + \
                        '\n.Material ' + '"' + Material + '"' + \
                        '\n.KeepWire "False"' + \
                        '\n.ConvertToSolidShape' + \
                    '\nEnd With' + \
                '\nEnd Sub'
    return SolidWire




def Poligon_2D(WGName, Points):
    """Create the 2D poligon for tRib waveguide

    Args:
        WGName (str): Name of the poligon
        Points (dict): Dictionary with the Points:
                        Points['X'] = []
                        Points['Y'] = []

    Returns:
        str: String with VBA Code 
    """

    WGName  = WGName
    # Extract the first points as starting Points
    Start_PointX = Points['X'][0]
    Start_PointY = Points['Y'][0]
    tmp = []
    tmp.append('Sub Main () ' + '\nWith Polygon' + '\n.Reset' + '\n.Name ' + '"' + WGName  + '"' + '\n.Curve ' + '"' + WGName  + '"' + '\n.Point ' + '"' + str(Start_PointX) + '"' + ',' + '"' + str(Start_PointY) + '"')
    for i in range(1, len(Points['X'])):
        b = '\n.LineTo ' + '"' + str(Points['X'][i]) + '"' + ',' + '"' + str(Points['Y'][i]) + '"'
        tmp.append(b)
    tmp.append('\n.Create' + '\nEnd With' + '\nEnd Sub')
    CurveData = ''.join(tmp)

    return CurveData




def Poligon_3D(WGName, Points):
    """Create the 2D poligon for tRib waveguide

    Args:
        WGName (str): Name of the poligon
        Points (dict): Dictionary with the Points:
                        Points['X'] = []
                        Points['Y'] = []
                        Points['Z'] = []   

    Returns:
        str: String with VBA Code 
    """

    WGName  = WGName
    # Extract the first points as starting Points
    Start_PointX = Points['X'][0]
    Start_PointY = Points['Y'][0]
    Start_PointZ = Points['Z'][0]

    tmp = []
    tmp.append('Sub Main () ' + '\nWith Polygon3D' + '\n.Reset' + '\n.Name ' + '"' + WGName  + '"' + '\n.Curve ' + '"' + WGName  + '"')
    for i in range(len(Points['X'])):
        b = '\n.Point ' + '"' + str(Points['X'][i]) + '"' + ',' + '"' + str(Points['Y'][i]) + '"'  + ',' + '"' + str(Points['Z'][i]) + '"'
        tmp.append(b)
    tmp.append('\n.Create' + '\nEnd With' + '\nEnd Sub')
    CurveData = ''.join(tmp)

    return CurveData








def RibWaveguide_ToSolid(SolidName, WaveguideName = "Rib_Waveguide", WG_Hight = None, Angle = None, NameFolder = None, Material = None, WGFolderName = None, WGName = None ):
    """This is ToSolid function that will allow the use to create the RibWaveguide. 
       TODO: Murrge the TOSolid and RibWaveguideToSolid functions later on. 

    Args:
        SolidName (str): Solid name
        WaveguideName (str, optional): Waveguide Name. Defaults to "Rib_Waveguide".
        WG_Hight (int/float, optional): Hight of the waveguide. Defaults to None.
        Angle (int/float, optional): Side angle of the waveguide. Defaults to None.
        NameFolder (str, optional): Folder name. Defaults to None.
        Material (str, optional): Material for of the Bond wire. For now you need to load the material 
                                  in your simulation and then use this function. Otherwise the material 
                                  will not be found. Defaults to None.
        WGFolderName (str, optional): Name of the folder where the poligon 3D or 2D is created. Defaults to None.
        WGName (str, optional): Name of the Poligon 3D or 2D. Defaults to None.

    Raises:
        ValueError: Error massage

    Returns:
        str: String with VBA Code 
    """

    #Check Curve Name
    if WaveguideName != "Rib_Waveguide":
        WaveguideName = WaveguideName
    else:
        pass

    # Folder Name
    if NameFolder != None:
        NameFolder = NameFolder
    else:
        NameFolder = WaveguideName

    # Material Def
    if Material == None:
        Material = "Copper (annealed)"
    else:
        Material = Material

    # Angle and Height. Pre define Hight = 2 and Angle = 30
    if WG_Hight == None:
        WG_Hight = 2
    else:
        WG_Hight = WG_Hight

    if Angle == None:
        Angle = 30
    else:
        Angle = Angle

    
     # Check if Curve Name and Folder Name are given
    if WGFolderName == None or WGName == None:
        raise ValueError("To create an Rib Waveguide you need to give the Waveguide Folder Name (WGFolderName) and the Rib Waveguide Name (WGName) that you already gave!")


    SolidWire = 'Sub Main () ' \
                    '\nWith ExtrudeCurve' + \
                        '\n.Reset' + \
                        '\n.Name ' + '"' + WaveguideName + '"' + \
                        '\n.Component ' + '"' + WaveguideName + '"' + \
                        '\n.Material ' + '"' + Material + '"' + \
                        '\n.Thickness ' + '"' + str(WG_Hight) +'"' + \
                        '\n.Twistangle "0"' + \
                        '\n.Taperangle ' + '"' + str(Angle) +'"' + \
                        '\n.Curve ' + '"' + WGFolderName + ':' + WGName + '"' + \
                        '\n.Create' + \
                    '\nEnd With' + \
                '\nEnd Sub'
    return SolidWire





def AddGlobalParameter(Parameters):
    """Add Global parameters to the working enviroment

    Args:
        Parameters (dict): Dictionary with the name and value of the parameter.
                           GlobalParam['Length']
    Returns:
        str: String with VBA Code 
    """

    # Parameters Values and Keys
    Param = list(Parameters.keys())
    Values = list(Parameters.values())

    tmp = []
    tmp.append('Sub Main () ')

    for i in range(len(Param)):
        b = '\nStoreParameter ' + '"' + str(Param[i]) + '"' + ',' + str(Values[i])
        tmp.append(b)
    tmp.append('\nEnd Sub')
    ParametersAdd = ''.join(tmp)
    return ParametersAdd




def DeleteGlobalParameter(Parameters):
    """Delete an Global parameter from the enviroment.

    Args:
        Parameters (dict): Dictionary with the name and value of the parameter.
                           GlobalParam['Length']

    Returns:
        str: String with VBA Code 
    """

    # Parameters Values and Keys
    Param = list(Parameters.keys())

    tmp = []
    tmp.append('Sub Main () ')

    for i in range(len(Param)):
        b = '\nDeleteParameter ' + '"' + str(Param[i]) + '"'
        tmp.append(b)
    tmp.append('\nEnd Sub')
    ParametersAdd = ''.join(tmp)
    return ParametersAdd





def BackgroundSet(Coordinates, Type = None):
    """Set the Simulation background.

    Args:
        Coordinates (dic): Dictionary with coordinates for the background. 
                            Coordinates["Xmin"] = int/float
                            Coordinates["Xmax"] = int/float
                            Coordinates["Ymin"] = int/float
                            Coordinates["Ymax"] = int/float
                            Coordinates["Zmin"] = int/float
                            Coordinates["Zmax"] = int/float
        Type (str, optional): Type of background material. Cehck CST for more information. Defaults to Normal.

    Returns:
        str: String with VBA Code 
    """

    #Check Background Type
    if Type == None:
        Type = "Normal"
    else:
        Type = Type

    BackgroundObj = 'Sub Main () '  \
                        '\nWith Background' + \
                            '\n.ResetBackground' + \
                            '\n.Type ' + '"' + Type + '"' + \
                            '\n.Epsilon "1.0"' + \
                            '\n.Mu "1.0"' + \
                            '\n.Rho "1.204"' + \
                            '\n.ThermalType "Normal"' + \
                            '\n.ThermalConductivity "0"' + \
                            '\n.SpecificHeat "1005"' + \
                            '\n.ApplyInAllDirections "False"' + \
                            '\n.XminSpace ' + '"' + str(Coordinates["Xmin"]) + '"' + \
                            '\n.XmaxSpace ' + '"' + str(Coordinates["Xmax"]) + '"' + \
                            '\n.YminSpace ' + '"' + str(Coordinates["Ymin"]) + '"' + \
                            '\n.YmaxSpace ' + '"' + str(Coordinates["Xmax"]) + '"' + \
                            '\n.ZminSpace ' + '"' + str(Coordinates["Zmin"]) + '"' + \
                            '\n.ZmaxSpace ' + '"' + str(Coordinates["Xmax"]) + '"' + \
                        '\nEnd With' + \
                        '\nEnd Sub'
    return BackgroundObj




def SetUnits(DictUnits):
    """Set the Simulation Global Units

    Args:
        DictUnits (dict): Dictionary with Units:
                          DictUnits['Dimensions'] - Length
                          DictUnits['Frequency'] - Frequency
                          DictUnits['Time'] - time
                          DictUnits['Temperature'] - Temperature
                          Not nessesery. Can be left out
                          DictUnits['Voltage'] - Voltage
                          DictUnits['Current'] - Current
                          DictUnits['Resistance'] - Resistance
                          DictUnits['Conductance'] - Conductance
                          DictUnits['Capacitance'] - Capacitance
                          DictUnits['Inductance'] - Inductance
    Returns:
        str: String with VBA Code 
    """

    #Check Values
    Keys = list(DictUnits.keys())
    tmp = []
    tmp.append('Sub Main () ' + '\nWith Units')
    
    if 'Dimensions' not in Keys:
        DictUnits['Dimensions'] = ("mm")
        tmp.append('\n.Geometry "mm"')
    else:
        tmp.append('\n.Geometry ' + '"' + str(DictUnits['Dimensions']) + '"')
    if 'Frequency' not in Keys:
        DictUnits['Frequency'] = ("hz")
        tmp.append('\n.Frequency "hz"')
    else:
        tmp.append('\n.Frequency ' + '"' + str(DictUnits['Frequency']) + '"')
    if 'Time' not in Keys:
        DictUnits['Time'] = ("s")
        tmp.append('\n.Time "s"')
    else:
        tmp.append('\n.Time ' + '"' + str(DictUnits['Time']) + '"')
    if 'Temperature' not in Keys:
        DictUnits['Temperature'] = 'degC'
        tmp.append('\n.Temperature "degC"')
    else:
        tmp.append('\n.Temperature ' + '"' + str(DictUnits['Temperature']) + '"')
    if 'Voltage' not in Keys:
        pass
    else:
        tmp.append('\n.Voltage ' + '("' + str(DictUnits['Voltage']) + '")')
    if 'Current' not in Keys:
        pass
        # DictUnits['Current'] = 'A'
    else:
        tmp.append('\n.Current ' + '("' + str(DictUnits['Current']) + '")')
    if 'Resistance' not in Keys:
        # DictUnits['Resistance'] = 'Ohm'
        pass
    else:
        tmp.append('\n.Resistance ' + '("' + str(DictUnits['Resistance']) + '")')
    if 'Conductance' not in Keys:
        # DictUnits['Conductance'] = 'Siemens'
        pass
    else:
        tmp.append('\n.Conductance ' + '("' + str(DictUnits['Conductance']) + '")')
    if 'Capacitance' not in Keys:
        # DictUnits['Capacitance'] =
        pass
    else:
        tmp.append('\n.Capacitance ' + '("' + str(DictUnits['Capacitance']) + '")')
    if 'Inductance' not in Keys:
        # DictUnits['Inductance'] = 'NanoH'
        pass
    else:
        tmp.append('\n.Inductance ' + '("' + str(DictUnits['Inductance']) + '")')
    # Units = 'Sub Main () ' \
    #             '\n With Units' + \
    #                 '\n.Geometry "mm"' + \
    #                 '\n.TemperatureUnit "Kelvin"' + \
    #                 '\n.Frequency "MHz"' + \
    #                 '\n.Time "ns"' + \
    #              '\n End With' + \
    #         '\n End Sub'


    # Units = 'Sub Main () ' \
    #             '\nWith Units' + \
    #                 '\n.Geometry ("m")'+ \
    #                 '\n.Frequency ("khz")'+ \
    #                 '\n.Time ("ns")'+ \
    #                 '\n.TemperatureUnit ("kelvin")'+ \
    #             '\nEnd With' + \
    #         '\nEnd Sub'
    tmp.append('\nEnd With' + '\nEnd Sub')
    Units = ''.join(tmp)
    return Units




def SetBoundary(Parameters):
    """Set Boundary box Parameters

    Args:
        Parameters (dict): Dictionary with Boundary Parameters 
                            Parameters["Xmin"] = str
                            Parameters["Xmax"] = str
                            Parameters["Ymin"] = str
                            Parameters["Ymax"] = str
                            Parameters["Zmin"] = str
                            Parameters["Zmax"] = str 
                            Parameters["Xsymmetry"] = str
                            Parameters["Ysymmetry"] = str
                            Parameters["Zsymmetry"] = str


    Returns:
        str: String with VBA Code 
    """
    # Boundary parameters
    Xmin = Parameters["Xmin"] 
    Xmax = Parameters["Xmax"] 
    Ymin = Parameters["Ymin"] 
    Ymax = Parameters["Ymax"] 
    Zmin = Parameters["Zmin"] 
    Zmax = Parameters["Zmax"] 
    Xsym = Parameters["Xsymmetry"] 
    Ysym = Parameters["Ysymmetry"] 
    Zsym = Parameters["Zsymmetry"] 


    BoundaryData = 'Sub Main () ' \
                    '\nWith Boundary' + \
                    '\n.Xmin ' + '"'+ str(Xmin) + '"' + \
                    '\n.Xmax ' + '"'+ str(Xmax) + '"' + \
                    '\n.Ymin ' + '"'+ str(Ymin) + '"' + \
                    '\n.Ymax ' + '"'+ str(Ymax) + '"' + \
                    '\n.Zmin ' + '"'+ str(Zmin) + '"' + \
                    '\n.Zmax ' + '"'+ str(Zmax) + '"' + \
                    '\n.Xsymmetry ' + '"'+ str(Xsym) + '"' + \
                    '\n.Ysymmetry ' + '"'+ str(Ysym) + '"' + \
                    '\n.Zsymmetry ' + '"'+ str(Zsym) + '"' + \
                    '\n.ApplyInAllDirections "False"' + \
                    '\n.XPeriodicShift "45.0"' + \
                    '\n.YPeriodicShift "0.0"' + \
                    '\n.ZPeriodicShift "0.0"' + \
                    '\n.PeriodicUseConstantAngles "True"' + \
                    '\n.SetPeriodicBoundaryAngles "30.0", "0.0"' + \
                    '\n.XminPotential ""' + \
                    '\n.XmaxPotential ""' + \
                    '\n.YminPotential ""' + \
                    '\n.YmaxPotential ""' + \
                    '\n.ZmaxPotential ""' + \
                '\nEnd With'  + \
           '\n End Sub'

    return BoundaryData



def WaveguidePort(Parameters):
    """Set the Waveguide Port for  

    Args:
        Parameters (dict): Dictionary with Port Parameters
                            Parameters["Orientation"] = Str with Port Orientation can be "Positive" or "Negative". For 
                                                        this function an 2 ports will be defined so please give an array with two Oriantations
                                                        like  Parameters["Orientation"] = ["Positive", "Positive"]
                            Parameters["Coordinates"] = Str witch Coordinates type, "Picks" is the best one!
                            Parameters["Span"] = Array with Array of port span [[Ymin, Ymax],[Zmin, Zmax]]
                            Parameters["Potential"] = Array with port Potential . For example [1,2]
                            Parameters["Port Number"] = Array with Port number [1,2]
                            Parameters["Polarity"] = Port Polarity can be be "Positive" or "Negative". For 
                                                    this function an 2 ports will be defined so please give an 
                                                    array with two Polaritys. For example Parameters["Polarity"] = ["Positive", "Positive"]
                            Parameters["Solid Name"] = Name of the Object on witch the Waveguide port will be created. For example "WG:WG1"
                            Parameters["Face ID"] = Array with the ID of the two picked faces. For example Parameters["Face ID"] = [2,4]

    Returns:
        str: String with VBA Code 
    """

    # Parameters to determin port position 
    Orientation = Parameters["Orientation"]
    Coordinates = Parameters["Coordinates"]
    # Choose coordinates
    Span11 = Parameters["Span"][0][0]
    Span12 = Parameters["Span"][0][1]
    Span21 = Parameters["Span"][1][0]
    Span22 = Parameters["Span"][1][1]
    Polarity = Parameters["Polarity"]
    SolidName = Parameters["Solid Name"]
    PickID = Parameters["Face ID"]
    #Return Ports in Dictionary for one Object 2 Ports can be define 
    Port = {}
    Port["1"] = None
    Port["2"] = None

    # Ports Numbers
    PortNum = Parameters["Port Number"]

    # loop true mode sets
    Potentials = Parameters["Potential"]


    #Start Loop for Ports
    for index, i in enumerate(Potentials,start = 1):
        data = 'Sub Main () ' \
                    '\nWith Port' + \
                    '\n.Reset' + \
                    '\n.PortNumber ' + '"' + str(PortNum[index-1]) + '"' + \
                    '\n.NumberOfModes "5"' + \
                    '\n.ReferencePlaneDistance "0"' + \
                    '\n.Coordinates ' + '"' + str(Coordinates) + '"' + \
                    '\n.Orientation ' +  '"' + str(Orientation[index-1]) + '"' + \
                    '\n.PortOnBound "False"' + \
                    '\n.ClipPickedPortToBound "False"' + \
                    '\n.YrangeAdd  ' + '"' + str(Span11[index-1]) + '"' + ',' + '"' + str(Span12[index-1]) + '"' + \
                    '\n.ZrangeAdd  ' + '"' + str(Span21[index-1]) + '"' + ',' + '"' + str(Span22[index-1]) + '"' + \
                    '\n.AdjustPolarization "True"' +\
                    '\n.PolarizationAngle "0"' +\
                    '\n.Create' + \
                '\nEnd With'  + \
                '\nEnd Sub'
        #  '\n.AddPotentialPicked  ' + '"' + str(index) + '"' + ',' + '"' + str(Polarity[index-1]) + '"' + ',' + '"' + str(SolidName) + '"'+ ',' + '"' + str(PickID[index-1]) + '"' + \
        Port[str(index)] = ''.join(data)

    return Port


def Pick(Parameters):
    """Pick function

    Args:
        Parameters (dict): Dictionary with Parameters   
                        PicParams["Option"] : Pcik option. For example "Face"
                        PicParams["Object"] : Name of the object on witch the pick will be executed.
                        PicParams["Face Number"] = int with the ID of the picked Face for example.
    Returns:
        str: String with VBA Code 
    """
    # Parameters for Picking an object
    FaceName = Parameters["Object"]
    Number = Parameters["Face Number"]

    if Parameters["Option"] == "Face":
        data = 'Sub Main () ' \
                    '\nPick.PickFaceFromId ' + '"' + str(FaceName) + '"' + "," +  '"' + str(Number)+ '"' \
                '\nEnd Sub'
    return data



def ClearAllPicks():
    """Clear all pciks

    Returns:
        str: String with VBA Code 
    """
    data = 'Sub Main () ' \
        '\nPick.ClearAllPicks ' +\
        '\nEnd Sub'
    return data




def SetSimFreqeuncy(Parameters):
    """Set Simulation Frequency

    Args:
        Parameters (dict): Dictionary with Parameters
                            Parameters["Min Frequency"] : Min Frequency of Simulation
                            Parameters["Max Frequency"] : Max Frequency of Simulation

    Returns:
        str: String with VBA Code 
    """
    FreqMin = Parameters["Min Frequency"]
    FreqMax = Parameters["Max Frequency"]
    
    data = 'Sub Main () ' \
        '\nWith Solver' + \
        '\n.Reset' + \
        '\n.FrequencyRange ' + '"' + str(FreqMin) + '"' + ',' + '"' + str(FreqMax) + '"' + \
        '\nEnd With' +\
        '\nEnd Sub'
    Port = ''.join(data)
    return Port



def SetSimWavelength(Parameters):
    """Set Simulation SetSimWavelength

    Args:
        Parameters (dict): Dictionary with Parameters
                            Parameters["Min Wavelength"] : Min Wavelength
                            Parameters["Max Wavelength"] : Max Wavelength
    Returns:
        str: String with VBA Code 
    """
    WavelengthMin = Parameters["Min Wavelength"]
    WavelengtMax = Parameters["Max Wavelength"]
    
    data = 'Sub Main () ' \
        '\nWith Solver' + \
        '\n.Reset' + \
        '\n.WavelengthRange ' + '"' + str(WavelengthMin) + '"' + ',' + '"' + str(WavelengtMax) + '"' + \
        '\nEnd With' +\
        '\nEnd Sub'
    Port = ''.join(data)
    return Port



def SetTimeSolver(Parameters):
    """Set Time solver Parameters

    Args:
        Parameters (dict): Dictionary with Parameters
                            Parameters["Accuracy"] : int from 10,15,20,25,30,35,40. Simulation accuracy. For example 20.    
                            Parameters["Caclculate Modes Only"] : Boolen. True if you want to calculate 
                            only the Ports modes. False to calculate the hole structure. 
                            Parameters["Auto Impedance"] : Booled. Set Port Impedance. True if you want to set manually
                            False otherwise.
                            Parameters["Impedance"] : int/float Port Impedance. 
                            Parameters["Source Port"] : Int. Set the Source Port. For example 1 or 2
    Returns:
        str: String with VBA Code 
    """
    Accuracy = Parameters["Accuracy"]
    ModesOnly = Parameters["Caclculate Modes Only"]
    AutoImpedance = Parameters["Auto Impedance"]
    Impedance = Parameters["Impedance"]
    Source = Parameters["Source Port"]

    data = 'Sub Main () ' \
        '\nWith Solver' + \
        '\n.Reset' + \
        '\n.Method "Hexahedral"' + \
        '\n.SteadyStateLimit ' + '"-' + str(Accuracy) + '"' + \
        '\n.CalculateModesOnly ' + '"' + str(ModesOnly) + '"' + \
        '\n.StimulationPort ' + '"' + str(Source) + '"' + \
        '\n.StimulationMode "All"' + \
        '\n.MeshAdaption "False"' + \
        '\n.SParaSymmetry "False"' +\
        '\n.AutoNormImpedance  ' + '"' + str(AutoImpedance) + '"' + \
        '\n.NormingImpedance ' + '"' + str(Impedance) + '"' + \
        '\n.SParaSymmetry "False"' +\
        '\nEnd With' + \
        '\nEnd Sub'
    Port = ''.join(data)
    return Port
        



def StartTimeSolver():
    """Start Time Simulation

    Returns:
        str: String with VBA Code 
    """
    data = 'Sub Main () ' \
        '\nWith Solver' + \
        '\n.Start' + \
        '\nEnd With' + \
        '\nEnd Sub'
    Port = ''.join(data)
    return Port




def CreateEfieldMonitor(Parameters):
    """Set Monitor

    Args:
        Parameters (dict): Dictionary with Parameters
                            Parameters["Wavelength"] : Int/float. This will set the wavelength. For this function
                            and only for this function you need to give the exact number. So 1.55 um will be Parameters["Wavelength"]  = 1.55e-6
                            Parameters["Monitor Type"] : str with monitor Type. Can be one of :
                                "Efield", "Hfield", "Surfacecurrent", "Powerflow", "Current",
                                "Powerloss", "Eenergy", "Elossdens", "Lossdens", "Henergy", 
                                "Farfield", "Fieldsource", "Spacecharge", "ParticleCurrentDensity", "Electrondensity" 
                            
    Raises:
        ValueError: Error massage

    Returns:
        str: String with VBA Code 
    """
    Wavelength = Parameters["Wavelength"]
    MonitorType = Parameters["Monitor Type"]
    Types = ["Efield", "Hfield", "Surfacecurrent", "Powerflow", "Current", "Powerloss", "Eenergy", "Elossdens", "Lossdens", "Henergy", "Farfield", "Fieldsource", "Spacecharge", "ParticleCurrentDensity", "Electrondensity" ]
    Speed_of_light = scipy.constants.c * 1e-6
    freq = Speed_of_light / Parameters["Wavelength"]

    if MonitorType in Types:
        Name = MonitorType + '_' + str(Wavelength)
        data = 'Sub Main ()' \
            '\nWith Monitor' + \
            '\n.Reset' + \
            '\n.Name ' + '"' + str(Name) + '"' + \
            '\n.Domain "Frequency"' + \
            '\n.FieldType ' + '"' + str(MonitorType) + '"' + \
            '\n.Frequency ' + '"' + str(freq) + '"' + \
            '\n.Create' + \
            '\nEnd With' + \
            '\nEnd Sub'
        Data = ''.join(data)
        return Data
    else:   
        raise ValueError("Parameters['Monitor Type'] is not in allowed types. Please choose one of the following types: ['Efield', 'Hfield', 'Surfacecurrent', 'Powerflow', 'Current', 'Powerloss', 'Eenergy', 'Elossdens', 'Lossdens', 'Henergy', 'Farfield', 'Fieldsource', 'Spacecharge', 'ParticleCurrentDensity', 'Electrondensity']")



