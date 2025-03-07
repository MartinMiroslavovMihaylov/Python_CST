import numpy as np



def Material(Name, Values):
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


def BondWire(NameWire, Coordinates, Height, Radius, BondwireType = "Spline", CenterPosition = 0.5, alpha = None, beta = None, Material = None, SolidWireModel = True, Termination = None, NameFolder = None):

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





def CurveWire(NameWire, Radius, Points = None, Material = None, SolidWireModel = True, Termination = None, NameFolder = None, CurveFolderName = None, CurveName = None):

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





def ToSolid(SolidName, CurveName = "Polygon", NameFolder = None, Material = None ):

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




def RibWG(WGName, Points):

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





def RibWaveguide_ToSolid(SolidName, WaveguideName = "Rib_Waveguide", WG_Hight = None, Angle = None, NameFolder = None, Material = None, WGFolderName = None, WGName = None ):

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




def Brick(BrickName, Coordinates, NameComponent = None, Material = None ):

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




def AddGlobalParameter(Parameters):

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




def SetBoundary():
    BoundaryData = 'Sub Main () ' \
                    '\nWith Boundary' + \
                    '\n.Xmin "electric"' + \
                    '\n.Xmax "electric"' + \
                    '\n.Ymin "electric"' + \
                    '\n.Ymax "electric"' + \
                    '\n.Zmin "electric"' + \
                    '\n.Zmax "electric"' + \
                    '\n.Xsymmetry "none"' + \
                    '\n.Ysymmetry "none"' + \
                    '\n.Zsymmetry "none"' + \
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



def WaveguidePort():
    Port = 'Sub Main () ' \
                '\nWith Port' + \
                '\n.Reset' + \
                '\n.PortNumber (1)' + \
                '\n.NumberOfModes (2)' + \
                '\n.AdjustPolarization (False)' + \
                '\n.PolarizationAngle (0.0)' + \
                '\n.ReferencePlaneDistance (0)' + \
                '\n.Coordinates ("Free")' + \
                '\n.Orientation ("xmax")' + \
                '\n.PortOnBound (False)' + \
                '\n.ClipPickedPortToBound (False)' + \
                '\n.Xrange (0, 0)' + \
                '\n.Yrange (-5, 5.4)' + \
                '\n.Zrange (-0.4, 0.4)' + \
                '\n.Create' + \
            '\nEnd With'  + \
            '\nEnd Sub'

    return Port




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

