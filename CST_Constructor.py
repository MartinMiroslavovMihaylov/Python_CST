import numpy as np 
import win32com.client
from matplotlib import pyplot as plt


def CstDefaultUnits(mws):
    units = mws.Units

    units.Geometry("mm")
    units.Frequency("kHz")
    units.TemperatureUnit("Kelvin")
    units.Time("ns")
    units.Voltage("V")
    units.Current("A")
    units.Resistance("Ohm")
    units.Conductance("Siemens")
    units.Capacitance("PikoF")
    units.Inductance("NanoH")

def CstDefineUnits(mws, Geometry, Frequency, Time, TemperatureUnit, Voltage, Current, Resistance, Conductance,
                   Capacitance, Inductance):
    units = mws.Units

    units.Geometry(Geometry)
    units.Frequency(Frequency)
    units.TemperatureUnit(TemperatureUnit)
    units.Time(Time)
    units.Voltage(Voltage)
    units.Current(Current)
    units.Resistance(Resistance)
    units.Conductance(Conductance)
    units.Capacitance(Capacitance)
    units.Inductance(Inductance)


def CstMeshInitiator(mws):
    FDSolver = mws.FDSolver
    mesh = mws.Mesh
    meshSetting = mws.MeshSettings
    meshAdaption3D = mws.MeshAdaption3D
    PostProcess1D = mws.PostProcess1D

    FDSolver.ExtrudeOpenBC('True')

    mesh.MergeThinPECLayerFixpoints('True')
    mesh.RatioLimit('20')
    mesh.AutomeshRefineAtPecLines('True', '6')
    mesh.FPBAAvoidNonRegUnite('True')
    mesh.ConsiderSpaceForLowerMeshLimit('False')
    mesh.MinimumStepNumber('5')
    mesh.AnisotropicCurvatureRefinement('True')
    mesh.AnisotropicCurvatureRefinementFSM('True')

    meshSetting.SetMeshType('Hex')
    meshSetting.Set('RatioLimitGeometry', '20')
    meshSetting.Set('EdgeRefinementOn', '1')
    meshSetting.Set('EdgeRefinementRatio', '6')

    meshSetting.SetMeshType('HexTLM')
    meshSetting.Set('RatioLimitGeometry', '20')

    meshSetting.SetMeshType('Tet')
    meshSetting.Set('VolMeshGradation', '1.5')
    meshSetting.Set('SrfMeshGradation', '1.5')

    meshAdaption3D.SetAdaptionStrategy('Energy')

    meshSetting.SetMeshType('Hex')
    meshSetting.Set('Version', '1%')

    mesh.MeshType('PBA')

    PostProcess1D.ActivateOperation('vswr', 'true')
    PostProcess1D.ActivateOperation('yz-matrices', 'true')



def CstDefineFrequencyRange(mws, frange1, frange2):
    Solver = mws.Solver
    mesh = mws.Mesh

    Solver.FrequencyRange(str(frange1), str(frange2))

    meshSettings = mws.MeshSettings
    meshSettings.SetMeshType('Hex')
    meshSettings.Set('Version', '1%')

    mesh.MeshType('PBA')



def CstDefineOpenBoundary(mws, minfrequency, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax):
    boundary = mws.Boundary
    plot = mws.Plot

    boundary.Xmin(Xmin)
    boundary.Xmax(Xmax)
    boundary.Ymin(Ymin)
    boundary.Ymax(Ymax)
    boundary.Zmin(Zmin)
    boundary.Zmax(Zmax)
    boundary.Xsymmetry('none')
    boundary.Ysymmetry('none')
    boundary.Zsymmetry('none')
    boundary.XminThermal('isothermal')
    boundary.XmaxThermal('isothermal')
    boundary.YminThermal('isothermal')
    boundary.YmaxThermal('isothermal')
    boundary.ZminThermal('isothermal')
    boundary.ZmaxThermal('isothermal')
    boundary.XsymmetryThermal('none')
    boundary.YsymmetryThermal('none')
    boundary.ZsymmetryThermal('none')
    boundary.ApplyInAllDirections('False')
    boundary.ApplyInAllDirectionsThermal('False')
    boundary.XminTemperature('')
    boundary.XminTemperatureType('None')
    boundary.XmaxTemperature('')
    boundary.XmaxTemperatureType('None')
    boundary.YminTemperature('')
    boundary.YminTemperatureType('None')
    boundary.YmaxTemperature('')
    boundary.YmaxTemperatureType('None')
    boundary.ZminTemperature('')
    boundary.ZminTemperatureType('None')
    boundary.ZmaxTemperature('')
    boundary.ZmaxTemperatureType('None')
    if Xmin == 'unit cell':
        boundary.XPeriodicShift('0.0')
        boundary.YPeriodicShift('0.0')
        boundary.ZPeriodicShift('0.0')
        boundary.PeriodicUseConstantAngles('False')
        boundary.SetPeriodicBoundaryAngles('0.0', '0.0')
        boundary.SetPeriodicBoundaryAnglesDirection('outward')
        boundary.UnitCellFitToBoundingBox('True')
        boundary.UnitCellDs1('0.0')
        boundary.UnitCellDs2('0.0')
        boundary.UnitCellAngle('90.0')
    if Xmin == 'expanded open':
        boundary.ReflectionLevel('0.0001')
        boundary.MinimumDistanceType('Fraction')
        boundary.MinimumDistancePerWavelengthNewMeshEngine('4')
        boundary.MinimumDistanceReferenceFrequencyType('CenterNMonitors')
        boundary.FrequencyForMinimumDistance(str(minfrequency))
        boundary.SetAbsoluteDistance('0.0')
        plot.DrawBox('True')




def Cstbrick(mws, Name, component, material, Xrange, Yrange, Zrange):
    brick = mws.Brick

    brick.Reset()
    brick.Name(Name)
    brick.component(component)
    brick.Material(material)
    brick.Xrange(str(Xrange[0]), str(Xrange[1]))
    brick.Yrange(str(Yrange[0]), str(Yrange[1]))
    brick.Zrange(str(Zrange[0]), str(Zrange[1]))
    brick.Create
    format(brick)




def Curve(mws, Name, Points):

    # Extract the first points as starting Points
    Start_PointX = Points['X'][0]
    Start_PointY = Points['X'][0]

    curve = mws.Polygon
    curve.Reset()
    curve.Name(Name)
    curve.Curve(Name)
    curve.Point(str(Start_PointX), str(Start_PointY))
    for i in range(1, len(Points['X'])):
        curve.LineTo(str(Points['X'][i]), str(Points['Y'][i]))
    curve.Create
    format(curve)


import sys
import win32com.client
import win32com.client.gencache
import os
import os.path
import datetime
sys.path.append(r"C:/Program Files (x86)/CST Studio Suite 2020/AMD64/python_cst_libraries/")


# win32com.client.Dispatch('CSTStudio.Application')
import win32com.client
import sys
sys.path.append(r"C:/Program Files (x86)/CST Studio Suite 2020/AMD64/python_cst_libraries")
import cst

cst =  win32com.client.Dispatch('CSTStudio.Application')

# from win32com.client import combrowse
# combrowse.main(modal=True)

screen = cst.GetObject("CSTStudio.Application").SelectedView.Control.Screen

# cst = win32com.client.dynamic.Dispatch("CSTStudio.Application")
cst.SetQuietMode(True)
new_mws = cst.NewMWS()
mws = cst.Active3D()


Geometry = 'm'
Frequency = 'GHz'
Time = 'ns'
TemperatureUnit = 'Kelvin'
Voltage = 'V'
Current = 'A'
Resistance = 'Ohm'
Conductance = 'S'
Capacitance = 'PikoF'
Inductance = 'NanoH'

CstDefineUnits(mws, Geometry, Frequency, Time, TemperatureUnit, Voltage, Current, Resistance, Conductance, Capacitance,
               Inductance)
CstMeshInitiator(mws)

CstDefineFrequencyRange(mws, 0.5, 4)


minfrequency = 0.5
x = 'expanded open'
CstDefineOpenBoundary(mws, minfrequency, x, x, x, x, x, x)


Name = 'Groundplane'
component = 'component1'
material = 'PEC'
Xrange = [-40, 40]
Yrange = [-40, 40]
Zrange = [0, 2]

Cstbrick(mws, Name, component, material, Xrange, Yrange, Zrange)


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

Name = 'Wire'
CurveName = 'CurveTest'
material = 'PEC'
Xrange = [-40, 40]
Yrange = [-40, 40]
Zrange = [0, 2]

Curve(mws, Name, Points)



units = mws.Units
units.Geometry("m")
units.Frequency("MHz")