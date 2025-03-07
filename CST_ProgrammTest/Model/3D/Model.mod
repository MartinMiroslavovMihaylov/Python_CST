'# MWS Version: Version 2025.0 - Aug 30 2024 - ACIS 34.0.1 -

'# length = mm
'# frequency = GHz
'# time = ns
'# frequency range: fmin = 1 fmax = 100
'# created = '[VERSION]2025.0|34.0.1|20240830[/VERSION]


'@ use template: Waveguide Coupler & Divider.cfg

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
'set the units
With Units
    .SetUnit "Length", "mm"
    .SetUnit "Frequency", "GHz"
    .SetUnit "Voltage", "V"
    .SetUnit "Resistance", "Ohm"
    .SetUnit "Inductance", "nH"
    .SetUnit "Temperature",  "degC"
    .SetUnit "Time", "ns"
    .SetUnit "Current", "A"
    .SetUnit "Conductance", "S"
    .SetUnit "Capacitance", "pF"
End With

ThermalSolver.AmbientTemperature "0"

'----------------------------------------------------------------------------

'set the frequency range
Solver.FrequencyRange "1", "100"

'----------------------------------------------------------------------------

With Background
     .Type "pec"
     .XminSpace "0.0"
     .XmaxSpace "0.0"
     .YminSpace "0.0"
     .YmaxSpace "0.0"
     .ZminSpace "0.0"
     .ZmaxSpace "0.0"
End With

' set boundary conditions to electric

With Boundary
     .Xmin "electric"
     .Xmax "electric"
     .Ymin "electric"
     .Ymax "electric"
     .Zmin "electric"
     .Zmax "electric"
     .Xsymmetry "none"
     .Ysymmetry "none"
     .Zsymmetry "none"
End With

Mesh.MinimumCurvatureRefinement "150"

With MeshSettings
     .SetMeshType "HexTLM"
     .Set "StepsPerWaveNear", "20"
     .Set "StepsPerBoxNear", "10"
     .Set "StepsPerWaveFar", "20"
     .Set "StepsPerBoxFar", "10"
     .Set "RatioLimitGeometry", "15"
End With

'----------------------------------------------------------------------------

Dim sDefineAt As String
sDefineAt = "1;50.5;100"
Dim sDefineAtName As String
sDefineAtName = "1;50.5;100"
Dim sDefineAtToken As String
sDefineAtToken = "f="
Dim aFreq() As String
aFreq = Split(sDefineAt, ";")
Dim aNames() As String
aNames = Split(sDefineAtName, ";")

Dim nIndex As Integer
For nIndex = LBound(aFreq) To UBound(aFreq)

Dim zz_val As String
zz_val = aFreq (nIndex)
Dim zz_name As String
zz_name = sDefineAtToken & aNames (nIndex)

' Define E-Field Monitors
With Monitor
    .Reset
    .Name "e-field ("& zz_name &")"
    .Dimension "Volume"
    .Domain "Frequency"
    .FieldType "Efield"
    .MonitorValue  zz_val
    .Create
End With

' Define H-Field Monitors
With Monitor
    .Reset
    .Name "h-field ("& zz_name &")"
    .Dimension "Volume"
    .Domain "Frequency"
    .FieldType "Hfield"
    .MonitorValue  zz_val
    .Create
End With

' Define Power flow Monitors
With Monitor
    .Reset
    .Name "power ("& zz_name &")"
    .Dimension "Volume"
    .Domain "Frequency"
    .FieldType "Powerflow"
    .MonitorValue  zz_val
    .Create
End With

Next

'----------------------------------------------------------------------------

With MeshSettings
     .SetMeshType "Hex"
     .Set "Version", 1%
End With

With Mesh
     .MeshType "PBA"
End With

'set the solver type
ChangeSolverType("HF Time Domain")

'----------------------------------------------------------------------------

'@ delete component: TestBrick

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Component.Delete "TestBrick"

'@ delete component: GND_Left_Electrode

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Component.Delete "GND_Left_Electrode"

'@ delete component: GND_right_Electrode

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Component.Delete "GND_right_Electrode"

'@ delete component: GND_Left_Electrode

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Component.Delete "GND_Left_Electrode"

'@ define brick: GND_Left_Electrode:solid1

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
With Brick
     .Reset 
     .Name "solid1" 
     .Component "GND_Left_Electrode" 
     .Material "Copper (annealed)" 
     .Xrange "-Length/2", "Length/2" 
     .Yrange "-Width_GND/2", "Width_GND/2" 
     .Zrange "0", "Hight" 
     .Create
End With

'@ delete component: component1

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Component.Delete "component1"

'@ delete component: default

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Component.Delete "default"

'@ delete component: GND_Left_Electrode

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Component.Delete "GND_Left_Electrode"

'@ delete component: Signal_Electrode

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Component.Delete "Signal_Electrode"

'@ define material: Folder1/LiNbO3

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
With Material 
     .Reset 
     .Name "LiNbO3"
     .Folder "Folder1"
     .Rho "0"
     .ThermalType "Normal"
     .ThermalConductivity "0"
     .SpecificHeat "0", "J/K/kg"
     .DynamicViscosity "0"
     .UseEmissivity "True"
     .Emissivity "0"
     .MetabolicRate "0.0"
     .VoxelConvection "0.0"
     .BloodFlow "0"
     .MechanicsType "Unused"
     .SolarRadiationAbsorptionType "Opaque"
     .Absorptance "0.0"
     .UseSemiTransparencyCalculator "False"
     .SolarRadTransmittance "0.0"
     .SolarRadReflectance "0.0"
     .SolarRadSpecimenThickness "0.0"
     .SolarRadRefractiveIndex "1.0"
     .SolarRadAbsorptionCoefficient "0.0"
     .IntrinsicCarrierDensityModel "none"
     .FrqType "all"
     .Type "Anisotropic"
     .MaterialUnit "Frequency", "GHz"
     .MaterialUnit "Geometry", "um"
     .MaterialUnit "Time", "ns"
     .MaterialUnit "Temperature", "degC"
     .EpsilonX "2.21"
     .EpsilonY "2.13"
     .EpsilonZ "2.21"
     .MuX "1"
     .MuY "1"
     .MuZ "1"
     .SigmaX "0"
     .SigmaY "0"
     .SigmaZ "0"
     .TanDX "0.0"
     .TanDY "0.0"
     .TanDZ "0.0"
     .TanDFreq "0.0"
     .TanDGiven "False"
     .TanDModel "ConstTanD"
     .SetConstTanDStrategyEps "AutomaticOrder"
     .ConstTanDModelOrderEpsX "3"
     .ConstTanDModelOrderEpsY "3"
     .ConstTanDModelOrderEpsZ "3"
     .DjordjevicSarkarUpperFreqEps "0"
     .SetElParametricConductivity "False"
     .ReferenceCoordSystem "Global"
     .CoordSystemType "Cartesian"
     .SigmaMX "0"
     .SigmaMY "0"
     .SigmaMZ "0"
     .TanDMX "0.0"
     .TanDMY "0.0"
     .TanDMZ "0.0"
     .TanDMFreq "0.0"
     .TanDMGiven "False"
     .TanDMModel "ConstTanD"
     .SetConstTanDStrategyMu "AutomaticOrder"
     .ConstTanDModelOrderMuX "3"
     .ConstTanDModelOrderMuY "3"
     .ConstTanDModelOrderMuZ "3"
     .DjordjevicSarkarUpperFreqMu "0"
     .SetMagParametricConductivity "False"
     .DispModelEps "None"
     .DispModelMu "None"
     .DispersiveFittingSchemeEps "Nth Order"
     .MaximalOrderNthModelFitEps "10"
     .ErrorLimitNthModelFitEps "0.1"
     .UseOnlyDataInSimFreqRangeNthModelEps "False"
     .DispersiveFittingSchemeMu "Nth Order"
     .MaximalOrderNthModelFitMu "10"
     .ErrorLimitNthModelFitMu "0.1"
     .UseOnlyDataInSimFreqRangeNthModelMu "False"
     .UseGeneralDispersionEps "False"
     .UseGeneralDispersionMu "False"
     .NLAnisotropy "False"
     .NLAStackingFactor "1"
     .NLADirectionX "1"
     .NLADirectionY "0"
     .NLADirectionZ "0"
     .Colour "0.839216", "0", "0" 
     .Wireframe "False" 
     .Reflection "False" 
     .Allowoutline "True" 
     .Transparentoutline "False" 
     .Transparency "0" 
     .Create
End With

'@ rename material: Folder1/LiNbO3 to: LiNbO3

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Material.Rename "Folder1/LiNbO3", "LiNbO3"

'@ delete material: Folder1

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Material.DeleteFolder "Folder1"

'@ define material colour: LiNbO3

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
With Material 
     .Name "LiNbO3"
     .Folder ""
     .Colour "0.839216", "0", "0" 
     .Wireframe "False" 
     .Reflection "False" 
     .Allowoutline "True" 
     .Transparentoutline "False" 
     .Transparency "0" 
     .ChangeColour 
End With

'@ delete curve: WG

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Curve.DeleteCurve "WG"

'@ delete curve: curve1

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Curve.DeleteCurve "curve1"

'@ delete curve: WG

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Curve.DeleteCurve "WG"

'@ delete wire folder: BondWire

'[VERSION]2025.0|34.0.1|20240830[/VERSION]
Wire.DeleteFolder "BondWire"

