import numpy as np 
import matplotlib.pyplot as plt
import os 
import sys
import cst
import cst.interface
import json
from scipy.integrate import quad
import scipy.integrate as integrate
from scipy.optimize import fsolve
from scipy.special import fresnel
plt.rcParams.update({"font.size":22})



class CST_Commands:

    def __init__(self):
        # Connect to existing CST Window or open new one 
        self.de = cst.interface.DesignEnvironment.connect_to_any_or_new()
        self.project_names_list = ["CS", "DS", "EMS", "FD3D", "MPS", "MWS", "PCBS", "PS"]
        self.prj = None
    


    def New_Project(self, name_project):
        # Start an new Project 
        if name_project in self.project_names_list:
            if name_project == "CS":
                # CST Cable Studio project
                self.prj = self.de.new_cs()
            elif name_project == "DS":
                # CST Design Studio project 
                self.prj = self.de.new_ds()
            elif name_project == "EMS":
                # CST EM Studio project
                self.prj = self.de.new_ems()
            elif name_project == "FD3D":
                # Filter Designer 3D project
                self.prj = self.de.new_fd3d()
            elif name_project == "MPS":
                # CST Mphysics Studio project
                self.prj = self.de.new_mps()
            elif name_project == "MWS":
                # CST Microwave Studio project
                self.prj = self.de.new_mws()
            elif name_project == "PCBS":
                # CST PCB Studio project
                self.prj = self.de.new_pcbs()
            elif name_project == "PS":
                # CST Particle Studio project
                self.prj = self.de.new_ps()
        else:
            raise ValueError(
                            """Possible new Project Names are:
                            'CS'   for CST Cable Studio project
                            'DS'   for CST Design Studio project
                            'EMS'  for CST EM Studio project
                            'FD3D' for Filter Designer 3D project
                            'MPS'  for CST Mphysics Studio project
                            'MWS'  for CST Microwave Studio project
                            'PCBS' for CST PCB Studio project
                            'PS'   for CST Particle Studio project"""
                            )
        
    

    def Save_Project(self, path, name, overwrite = False):
        self.prj.save(path + "/" + str(name) + ".cst", allow_overwrite=overwrite)



    def Open_Project(self, path):
        # open the project
        self.prj = self.de.open_project(path)



    ############################################################################
    # CST solver parameters
    ############################################################################
    

    def set_Units(self, Parameters):

        u_Lenght = Parameters["Unit Lenght"]
        u_Freq = Parameters["Unit Frequency"]
        u_Time = Parameters["Unit Time"]
        u_Temp = Parameters["Unit Temperature"]
        # define a String holding the VBA code for setting the units

        unit_vba = f"""
                        With Units
                            .SetUnit ("Length", "{u_Lenght}")
                            .SetUnit ("Frequency","{u_Freq}")
                            .SetUnit ("Time", "{u_Time}")
                            .SetUnit ("Temperature", "{u_Temp}")
                        End With
                        """
        
        # add the VBA code to the History List
        self.prj.model3d.add_to_history("set units", unit_vba)



    
    def add_json_material(self, path, name):
        # load materials from JSON file
        with open(os.path.join(path, str(name) + ".json")) as json_file:
            materials = json.load(json_file)



    def add_material(self, name):
        list_of_materials = ["Si", "LiNbO3", "SiO2", "Au", "Gold", "Al", "Glue"]

        if name in list_of_materials:
            if name == "Si":
                self.add_Si()
            elif name == "LiNbO3":
                Data = {}
                Data["X"] = 4.906
                Data["Y"] = 4.584095
                Data["Z"] = 4.906
                self.add_anisotropy_material("LiNbO3", Data)
            elif name == "SiO2":
                self.add_SiO2()
            elif name == "Au" or "Gold":
                self.add_Au()
            elif name == "Al":
                self.add_Al()
            elif name == "Glue":
                self.add_Glue()
            elif name == "ProbeMaterials":
                self.add_GGB_Probe_Material()
        else:
            raise ValueError(
                            """
                            For now only the following materials can be added with this python script:
                            "Si"
                            "LiNbO3"
                            "SiO2"
                            "Au or Gold"
                            "Al"
                            "Glue"
                            ProbeMaterials
                            """
                            )





    def add_anisotropy_material(self, name, Values):
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
        # Chix = round((Values["X"] - 1), 2)
        # Chiy = round((Values["Y"] - 1), 2)
        # Chiz = round((Values["Z"] - 1), 2 )
        Chix = Values["X"] 
        Chiy = Values["Y"] 
        Chiz = Values["Z"]
        Name = name


        vba_code = f"""
                        With Material
                        .Reset
                        .Name "{Name}"
                        .FrqType "All"
                        .Type "Anisotropic"
                        .SetMaterialUnit "GHz", "um"
                        .MuX "1"
                        .MuY "1"
                        .MuZ "1"
                        .EpsilonX "{Values["X"]}"
                        .EpsilonY "{Values["Y"]}"
                        .EpsilonZ "{Values["Z"]}"'
                        .DispModelEps "nonlinear2nd"
                        .DispCoeff0EpsX "{Chix}"
                        .DispCoeff0EpsY "{Chiy}"
                        .DispCoeff0EpsZ "{Chiz}"
                        .DispModelMu "None"
                        .AddDispEpsPole1stOrderX "0","0"
                        .Color "1","0", "0"
                        .Create
                    End With
                """
        self.prj.model3d.add_to_history("add LiNbO3 material", vba_code)


    def add_Si(self):
        """Add silicon to the Material Library

        Args:
            Name (str): Name of the Material

        Returns:
            str: String with the VBA code
        """
        # Add Silicon to materials
        vba_code = f"""
                        With Material
                        .Reset
                        .Name "Si"
                        .FrqType "All"
                        .Type "Normal"
                        .SetMaterialUnit "GHz", "um"
                        .Mu "1"
                        .Epsilon "11.9"
                        .Rho "2330"
                        .ThermalConductivity "148"
                        .ThermalType "Normal"
                        .MechanicsType "Normal"
                        .SpecificHeat "700"
                        .YoungsModulus "112"
                        .PoissonsRatio "0.28"
                        .Color "0.5","0.5", "0.5"
                        .Create
                        End With
                    """
        self.prj.model3d.add_to_history("add Si material", vba_code)
    

    def add_SiO2(self):
        """Add silicon dioxide to the Material Library

        Args:
            Name (str): Name of the Material

        Returns:
            str: String with the VBA code
        """
        # Add Meterial SiO2
        vba_code = f"""
                        With Material
                        .Reset
                        .Name "SiO2"
                        .FrqType "All"
                        .Type "Normal"
                        .SetMaterialUnit "GHz", "um"
                        .Mu "1"
                        .Epsilon "3.9"
                        .Color "0.529", "0.808", "0.922"
                        .Create
                        End With
                    """
        self.prj.model3d.add_to_history("add SiO2 material", vba_code)




    def add_Au(self):
        """Add gold to the Material Library
        Args:
            Name (str): Name of the Material

        Returns:
            str: String with the VBA code
        """
        # Add Meterial Gold.
        vba_code = f"""
                        With Material
                        .Reset
                        .Name "Au"
                        .FrqType "All"
                        .Type "lossy metal"
                        .SetMaterialUnit "GHz", "um"
                        .Mu "1"
                        .Sigma "4.561e+007"
                        .Rho "19320"
                        .ThermalConductivity "130"
                        .ThermalType "Normal"
                        .SpecificHeat "700"
                        .YoungsModulus "78"
                        .PoissonsRatio "42"
                        .Color "1","0.84", "0"
                        .Create
                        End With
                    """
        self.prj.model3d.add_to_history("add Au material", vba_code)

    def add_Al(self):
        """Add Aluminium to the Material Library
        Args:
            Name (str): Name of the Material

        Returns:
            str: String with the VBA code
        """
        # Add Meterial Gold.
        vba_code = f"""
                        With Material
                            .Reset
                            .Name "Al"
                            .Folder ""
                            .FrqType "static"
                            .Type "Normal"
                            .SetMaterialUnit "Hz", "mm"
                            .Epsilon "1"
                            .Mu "1.0"
                            .Kappa "3.56e+007"
                            .TanD "0.0"
                            .TanDFreq "0.0"
                            .TanDGiven "False"
                            .TanDModel "ConstTanD"
                            .KappaM "0"
                            .TanDM "0.0"
                            .TanDMFreq "0.0"
                            .TanDMGiven "False"
                            .TanDMModel "ConstTanD"
                            .DispModelEps "None"
                            .DispModelMu "None"
                            .DispersiveFittingSchemeEps "General 1st"
                            .DispersiveFittingSchemeMu "General 1st"
                            .UseGeneralDispersionEps "False"
                            .UseGeneralDispersionMu "False"
                            .FrqType "all"
                            .Type "Lossy metal"
                            .MaterialUnit "Frequency", "GHz"
                            .MaterialUnit "Geometry", "mm"
                            .MaterialUnit "Time", "s"
                            .MaterialUnit "Temperature", "Kelvin"
                            .Mu "1.0"
                            .Sigma "3.56e+007"
                            .Rho "2700.0"
                            .ThermalType "Normal"
                            .ThermalConductivity "237.0"
                            .SpecificHeat "900", "J/K/kg"
                            .MetabolicRate "0"
                            .BloodFlow "0"
                            .VoxelConvection "0"
                            .MechanicsType "Isotropic"
                            .YoungsModulus "69"
                            .PoissonsRatio "0.33"
                            .ThermalExpansionRate "23"
                            .ReferenceCoordSystem "Global"
                            .CoordSystemType "Cartesian"
                            .NLAnisotropy "False"
                            .NLAStackingFactor "1"
                            .NLADirectionX "1"
                            .NLADirectionY "0"
                            .NLADirectionZ "0"
                            .Colour "1", "1", "0"
                            .Wireframe "False"
                            .Reflection "False"
                            .Allowoutline "True"
                            .Transparentoutline "False"
                            .Transparency "0"
                            .Create
                        End With

                    """
        self.prj.model3d.add_to_history("add Al material", vba_code)


    def add_Glue(self):
        # Add Material Glue
        vba_code = f"""
                     With Material 
                    .Reset 
                    .Name "DAF_Glue"
                    .Folder ""
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
                    .Absorptance "0"
                    .MechanicsType "Unused"
                    .IntrinsicCarrierDensity "0"
                    .FrqType "all"
                    .Type "Normal"
                    .MaterialUnit "Frequency", "GHz"
                    .MaterialUnit "Geometry", "um"
                    .MaterialUnit "Time", "ns"
                    .MaterialUnit "Temperature", "degC"
                    .Epsilon "2.6"
                    .Mu "1"
                    .Sigma "0"
                    .TanD "0.0"
                    .TanDFreq "0.0"
                    .TanDGiven "False"
                    .TanDModel "ConstTanD"
                    .SetConstTanDStrategyEps "AutomaticOrder"
                    .ConstTanDModelOrderEps "3"
                    .DjordjevicSarkarUpperFreqEps "0"
                    .SetElParametricConductivity "False"
                    .ReferenceCoordSystem "Global"
                    .CoordSystemType "Cartesian"
                    .SigmaM "0"
                    .TanDM "0.0"
                    .TanDMFreq "0.0"
                    .TanDMGiven "False"
                    .TanDMModel "ConstTanD"
                    .SetConstTanDStrategyMu "AutomaticOrder"
                    .ConstTanDModelOrderMu "3"
                    .DjordjevicSarkarUpperFreqMu "0"
                    .SetMagParametricConductivity "False"
                    .DispModelEps  "None"
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
                    .Colour "0", "1", "1" 
                    .Wireframe "False" 
                    .Reflection "False" 
                    .Allowoutline "True" 
                    .Transparentoutline "False" 
                    .Transparency "0" 
                    .Create
                End With
                """
        self.prj.model3d.add_to_history("add Glue material", vba_code)
        
        
        
    def add_GGB_Probe_Material(self):
        # Define Probes Materials
        vba_code1 = f"""
                    With Material 
                        .Reset 
                        .Name "Teflon (PTFE) (loss free)" 
                        .Folder "" 
                        .Rho "2200.0"
                        .ThermalType "Normal"
                        .ThermalConductivity "0.2"
                        .SpecificHeat "1000", "J/K/kg"
                        .DynamicViscosity "0"
                        .UseEmissivity "True"
                        .Emissivity "0"
                        .MetabolicRate "0.0"
                        .VoxelConvection "0.0"
                        .BloodFlow "0"
                        .Absorptance "0"
                        .MechanicsType "Isotropic"
                        .YoungsModulus "0.5"
                        .PoissonsRatio "0.4"
                        .ThermalExpansionRate "140"
                        .IntrinsicCarrierDensity "0"
                        .FrqType "all"
                        .Type "Normal"
                        .MaterialUnit "Frequency", "GHz"
                        .MaterialUnit "Geometry", "mm"
                        .MaterialUnit "Time", "s"
                        .Epsilon "2.1"
                        .Mu "1.0"
                        .Sigma "0.0"
                        .TanD "0.0"
                        .TanDFreq "0.0"
                        .TanDGiven "False"
                        .TanDModel "ConstTanD"
                        .SetConstTanDStrategyEps "AutomaticOrder"
                        .ConstTanDModelOrderEps "1"
                        .DjordjevicSarkarUpperFreqEps "0"
                        .SetElParametricConductivity "False"
                        .ReferenceCoordSystem "Global"
                        .CoordSystemType "Cartesian"
                        .SigmaM "0.0"
                        .TanDM "0.0"
                        .TanDMFreq "0.0"
                        .TanDMGiven "False"
                        .TanDMModel "ConstSigma"
                        .SetConstTanDStrategyMu "AutomaticOrder"
                        .ConstTanDModelOrderMu "1"
                        .DjordjevicSarkarUpperFreqMu "0"
                        .SetMagParametricConductivity "False"
                        .DispModelEps  "None"
                        .DispModelMu "None"
                        .DispersiveFittingSchemeEps "1st Order"
                        .DispersiveFittingSchemeMu "1st Order"
                        .UseGeneralDispersionEps "False"
                        .UseGeneralDispersionMu "False"
                        .NLAnisotropy "False"
                        .NLAStackingFactor "1"
                        .NLADirectionX "1"
                        .NLADirectionY "0"
                        .NLADirectionZ "0"
                        .Colour "0.75", "0.95", "0.85" 
                        .Wireframe "False" 
                        .Reflection "False" 
                        .Allowoutline "True" 
                        .Transparentoutline "False" 
                        .Transparency "0" 
                        .Create
                    End With 

                    With Material 
                        .Reset 
                        .Name "epoxy_casting_CR110" 
                        .Folder "" 
                        .Rho "0"
                        .ThermalType "Normal"
                        .ThermalConductivity "0"
                        .SpecificHeat "0", "J/K/kg"
                        .DynamicViscosity "0"
                        .UseEmissivity "True"
                        .Emissivity "0"
                        .MetabolicRate "0"
                        .VoxelConvection "0"
                        .BloodFlow "0"
                        .Absorptance "0"
                        .MechanicsType "Unused"
                        .IntrinsicCarrierDensity "0"
                        .FrqType "all"
                        .Type "Normal"
                        .MaterialUnit "Frequency", "GHz"
                        .MaterialUnit "Geometry", "um"
                        .MaterialUnit "Time", "ns"
                        .MaterialUnit "Temperature", "K"
                        .Epsilon "2.5"
                        .Mu "1"
                        .Sigma "0"
                        .TanD "0.0"
                        .TanDFreq "0.0"
                        .TanDGiven "False"
                        .TanDModel "ConstTanD"
                        .SetConstTanDStrategyEps "AutomaticOrder"
                        .ConstTanDModelOrderEps "1"
                        .DjordjevicSarkarUpperFreqEps "0"
                        .SetElParametricConductivity "False"
                        .ReferenceCoordSystem "Global"
                        .CoordSystemType "Cartesian"
                        .SigmaM "0"
                        .TanDM "0.0"
                        .TanDMFreq "0.0"
                        .TanDMGiven "False"
                        .TanDMModel "ConstTanD"
                        .SetConstTanDStrategyMu "AutomaticOrder"
                        .ConstTanDModelOrderMu "1"
                        .DjordjevicSarkarUpperFreqMu "0"
                        .SetMagParametricConductivity "False"
                        .DispModelEps  "None"
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
                        .NonlinearMeasurementError "1e-1"
                        .NLAnisotropy "False"
                        .NLAStackingFactor "1"
                        .NLADirectionX "1"
                        .NLADirectionY "0"
                        .NLADirectionZ "0"
                        .Colour "0.501961", "0.501961", "0" 
                        .Wireframe "False" 
                        .Reflection "False" 
                        .Allowoutline "True" 
                        .Transparentoutline "False" 
                        .Transparency "0" 
                        .Create
                    End With 

        """
        self.prj.model3d.add_to_history(f"create Material Teflon (PTFE) (loss free) and epoxy_casting_CR110 for GGB Probe", vba_code1)



    def AddGlobalParameter(self, Parameters):
        """
        Add Global parameters to the working enviroment

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
        for i in range(len(Param)):
            tmp.append(f'StoreParameter "{Param[i]}", {Values[i]}')

        # join them with newlines
        line_block = "\n".join(tmp)
        
        vba_code = f"""
                    {line_block}
                    """
        self.prj.model3d.add_to_history("add global parameter", vba_code)


    def DeleteGlobalParameter(self, Parameters):
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
        for i in range(len(Param)):
            tmp.append(f'DeleteParameter "{Param[i]}"')

        # join them with newlines
        line_block = "\n".join(tmp)
        
        vba_code = f"""
                    {line_block}
                    """
        self.prj.model3d.add_to_history("delete global parameter", vba_code)

        
    def setBackground(self, Parameters):
        """Set the Simulation background.

        Args:
            Parameters (dic): Dictionary with coordinates for the background. 
                                Parameters["Type Background"] = str, Type of the Background. Defoult is "Normal"
                                Parameters["Xmin Background"] = int/float
                                Parameters["Xmin Background"] = int/float
                                Parameters["Xmin Background"] = int/float
                                Parameters["Xmin Background"] = int/float
                                Parameters["Xmin Background"] = int/float
                                Parameters["Xmin Background"] = int/float
            Type (str, optional): Type of background material. Cehck CST for more information. Defaults to Normal.

        Returns:
            str: String with VBA Code 
        """

        #Check Background Type
        if Parameters["Type Background"] == None:
            Type = "Normal"
        else:
            Type = Parameters["Type Background"]

        vba_code = f"""
                    With Background
                    .ResetBackground
                    .Type "{Type}"
                    .Epsilon "1.0"
                    .Mu "1.0"
                    .Rho "1.204"
                    .ThermalType "Normal"
                    .ThermalConductivity "0"
                    .SpecificHeat "1005"
                    .ApplyInAllDirections "False"
                    .XminSpace "{Parameters["Xmin Background"]}"
                    .XmaxSpace "{Parameters["Xmax Background"]}"
                    .YminSpace "{Parameters["Ymin Background"]}"
                    .YmaxSpace "{Parameters["Ymax Background"]}"
                    .ZminSpace "{Parameters["Zmin Background"]}"
                    .ZmaxSpace "{Parameters["Zmax Background"]}"
                    End With
                    """
        self.prj.model3d.add_to_history("set background", vba_code)       


    def setBoundary(self, Parameters):
        """Set Boundary box Parameters

        Args:
            Parameters (dict): Dictionary with Boundary Parameters 
                                Parameters["Xmin Boundary"] = str
                                Parameters["Xmax Boundary"] = str
                                Parameters["Ymin Boundary"] = str
                                Parameters["Ymax Boundary"] = str
                                Parameters["Zmin Boundary"] = str
                                Parameters["Zmax Boundary"] = str 
                                Parameters["Xsymmetry Boundary"] = str
                                Parameters["Ysymmetry Boundary"] = str
                                Parameters["Zsymmetry Boundary"] = str


        Returns:
            str: String with VBA Code 
        """
        # Boundary parameters
        Xmin = Parameters["Xmin Boundary"]
        Xmax = Parameters["Xmax Boundary"]
        Ymin = Parameters["Ymin Boundary"]
        Ymax = Parameters["Ymax Boundary"]
        Zmin = Parameters["Zmin Boundary"]
        Zmax = Parameters["Zmax Boundary"]
        Xsym = Parameters["Xsymmetry Boundary"]
        Ysym = Parameters["Ysymmetry Boundary"] 
        Zsym = Parameters["Zsymmetry Boundary"]

        vba_code = f"""
                    With Boundary
                    .Xmin "{Xmin}"
                    .Xmax "{Xmax}"
                    .Ymin "{Ymin}"
                    .Ymax "{Ymax}"
                    .Zmin "{Zmin}"
                    .Zmax "{Zmax}"
                    .Xsymmetry "{Xsym}"
                    .Ysymmetry "{Ysym}"
                    .Zsymmetry "{Zsym}"
                    .ApplyInAllDirections "False"
                    .XPeriodicShift "45.0"
                    .YPeriodicShift "0.0"
                    .ZPeriodicShift "0.0"
                    .PeriodicUseConstantAngles "True"
                    .SetPeriodicBoundaryAngles "30.0", "0.0"
                    .XminPotential ""
                    .XmaxPotential ""
                    .YminPotential ""
                    .YmaxPotential ""
                    .ZmaxPotential ""
                    End With
                    """
        self.prj.model3d.add_to_history("set boundary", vba_code)


    def setSimFreqeuncy(self, Parameters):
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
        
        vba_code = f"""
                    With Solver
                    .Reset
                    .FrequencyRange "{FreqMin}", "{FreqMax}"
                    End With
                    """
        self.prj.model3d.add_to_history("set sim frequency", vba_code)

    
    def setSimWavelength(self, Parameters):
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
        
        vba_code = f"""
                    With Solver
                    .Reset
                    .WavelengthRange "{WavelengthMin}", "{WavelengtMax}"
                    End With
                    """
        self.prj.model3d.add_to_history("set sim wavelength", vba_code)
        

   


############################################################################
# Building objects
############################################################################


    def Brick(self, Parameters):
        
        LenghtMin = Parameters["Brick Lenght Min"]
        LenghtMax = Parameters["Brick Lenght Max"]
        WidthMin = Parameters["Brick Width Min"]
        WidthMax = Parameters["Brick Width Max"]
        HightMin = Parameters["Brick Hight Min"]
        HightMax = Parameters["Brick Hight Max"]
        Name = Parameters["Brick Name"]
        Component_name = Parameters["Component Name"]

        Material = Parameters.get("Material", "PEC")

        vba_code = f"""
                    With Brick
                        .Reset
                        .Name "{Name}"
                        .Component "{Component_name}"
                        .Material "{Material}"
                        .Xrange "{LenghtMin/2}", "{LenghtMax/2}"
                        .Yrange "{WidthMin/2}", "{WidthMax/2}"
                        .Zrange "{HightMin/2}", "{HightMax/2}"
                        .Create
                    End With
                    """
        self.prj.model3d.add_to_history(f"create brick {Name}", vba_code)
        
        
    def Sphere(self, Parameters):
        Name = Parameters["Name"]
        ComponentName = Parameters["Component Name"]
        Material = Parameters["Material"]
        Axis = Parameters["Axis"]
        CentRadius = Parameters["Center Radius"]
        Centers = Parameters["Center Positions"]
        
        if "Top Radius" not in Parameters.keys():
            TopRadius = 0
        else:
            TopRadius = Parameters["Top Radius"]
        if "Bottom Radius" not in Parameters.keys():
            BotRadius = 0
        else:
            BotRadius = Parameters["Bottom Radius"]
            
                    
        vba_code = f"""
                    With Sphere 
                     .Reset 
                     .Name "{Name}" 
                     .Component "{ComponentName}" 
                     .Material "{Material}" 
                     .Axis "{Axis}" 
                     .CenterRadius "{CentRadius}" 
                     .TopRadius "{TopRadius}" 
                     .BottomRadius "{BotRadius}" 
                     .Center "{Centers["X"]}", "{Centers["X"]}", "{Centers["X"]}" 
                     .Segments "0" 
                     .Create 
                End With
                """
                
        self.prj.model3d.add_to_history(f"create sphere {Name}", vba_code)
        
        


        

    def Curve(self, CurveName, Points):

        CurveName  = CurveName
        # Extract the first points as starting Points
        Start_PointX = Points['X'][0]
        Start_PointY = Points['X'][0]
        lines = []

        for i in range(1, len(Points['X'])):
            lines.append(f'.LineTo "{Points["X"][i]}", "{Points["Y"][i]}"')

        # join them with newlines
        line_block = "\n".join(lines)

        vba_code = f"""
                    With Polygon 
                        .Reset 
                        .Name "{CurveName}" 
                        .Curve "{CurveName}"
                        .Point "{Start_PointX}", "{Start_PointY}" 
                        {line_block}
                        .Create
                    End With
                """
        self.prj.model3d.add_to_history("create curve", vba_code)
        
        
    
    def Elipse(self, Parameters):
        Name = Parameters["Name"]
        NameCurve = Parameters["Curve Name"]
        
        vba_code = f"""
                    With Ellipse
                     .Reset 
                     .Name "{Name}" 
                     .Curve "{NameCurve}" 
                     .XRadius "{Parameters["X Radius"]}"
                     .YRadius "{Parameters["Y Radius"]}" 
                     .Xcenter "{Parameters["X Center"]}" 
                     .Ycenter "{Parameters["Y Center"]}" 
                     .Segments "0" 
                     .Create
                End With
                
                        """
        self.prj.model3d.add_to_history(f"create elipse curve {Name}", vba_code)
        



    def Poligon_2D(self, WGName, Points):
        """
        Create the 2D poligon for tRib waveguide
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
        lines = []

        for i in range(1, len(Points['X'])):
            lines.append(f'.LineTo "{Points["X"][i]}", "{Points["Y"][i]}"')

        # join them with newlines
        line_block = "\n".join(lines)
        
        vba_code = f"""
                    With Polygon
                    .Reset
                    .Name "{WGName}"
                    .Curve "{WGName}"
                    .Point "{Start_PointX}", "{Start_PointY}"
                    {line_block}
                    .Create
                    End With
                    """

        self.prj.model3d.add_to_history("2d polygon", vba_code)  
            



    def Poligon_3D(self, Parameters):
        """
        Create the 2D poligon for tRib waveguide

        Args:
            WGName (str): Name of the poligon
            Points (dict): Dictionary with the Points:
                            Points['X'] = []
                            Points['Y'] = []
                            Points['Z'] = []   

        Returns:
            str: String with VBA Code 
        """

        Name = Parameters["Name"] 
        CurveName = Parameters["Curve Name"] 
        lines = []

        for i in range(0, len(Parameters["Point"]["X"])):
            lines.append(f'.Point "{Parameters["Point"]["X"][i]}", "{Parameters["Point"]["Y"][i]}", "{Parameters["Point"]["Z"][i]}"')

        # join them with newlines
        line_block = "\n".join(lines)

        vba_code = f"""
                    With Polygon3D
                    .Reset
                    .Name "{Name}"
                    .Curve "{CurveName}"
                    {line_block}
                    .Create
                    End With
                    """ 
        self.prj.model3d.add_to_history(f"3d polygon {CurveName}", vba_code) 




    
    # MZM Design
    def MZM(self, Parameters):
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
        

        self.add_anisotropy_material("LiNbO3", Data)
        self.add_SiO2()
        self.add_Au()


        # Calculate Top and Bottom Widths of the Waveguide
        x = abs(Height_WG / (np.cos((Angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - Height_WG ** 2)
        Width_WG_New = Width_WG + 2 * extention 

        # Global Parameters
        WidthObject = 2*Width_WG_New + Width_Signal + 2*Width_Electrodes + 4*Gap


        # Substrate Definition 
        Parameters = {}
        Parameters["Brick Lenght Min"] = -Length_WG/2
        Parameters["Brick Lenght Max"] = Length_WG/2
        Parameters["Brick Width Min"] = -WidthObject/2 
        Parameters["Brick Width Max"] = WidthObject/2
        Parameters["Brick Hight Min"] = -Height_Substrate/2
        Parameters["Brick Hight Max"] = Height_Substrate/2
        Parameters["Brick Name"] = "LNOI_Substrate"
        Parameters["Component Name"] = "LNOI_Substrate"
        Parameters["Material"] = "SiO2"
        self.Brick(Parameters)
     



        # Slab definition
        if Heigh_Slab == 0:
            HeightZ = Height_Substrate/2
        else:
            Parameters = {}
            Parameters["Brick Lenght Min"] = -Length_WG/2
            Parameters["Brick Lenght Max"] = Length_WG/2
            Parameters["Brick Width Min"] = -WidthObject/2 
            Parameters["Brick Width Max"] = WidthObject/2
            Parameters["Brick Hight Min"] = Height_Substrate/2
            Parameters["Brick Hight Max"] = Height_Substrate/2 + Heigh_Slab
            Parameters["Brick Name"] = "LNOI_Slab"
            Parameters["Component Name"] = "LNOI_Slab"
            Parameters["Material"] = "LiNbO3"
            self.Brick(Parameters)
            HeightZ = Height_Substrate/2 + Heigh_Slab

            
        # Calculate Top and Bottom Widths of the Waveguide
        x = abs(Height_WG / (np.cos((Angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - Height_WG ** 2)
        Width_WG_New = Width_WG + 2 * extention 

        #  Electrodes 
        NamesElectrodes = "Electrode_Left", "Electrode_Right", "Signal"

        # Plase Electrodes
        Parameters = {}
        Parameters["Brick Lenght Min"] = -Length_WG/2
        Parameters["Brick Lenght Max"] = Length_WG/2
        Parameters["Brick Width Min"] = -WidthObject/2 + Width_Electrodes
        Parameters["Brick Width Max"] = -WidthObject/2
        Parameters["Brick Hight Min"] = HeightZ
        Parameters["Brick Hight Max"] = HeightZ + Height_Electrodes
        Parameters["Brick Name"] = NamesElectrodes[0]
        Parameters["Component Name"] = NamesElectrodes[0]
        Parameters["Material"] = "Au"
        self.Brick(Parameters)
        


        Parameters = {}
        Parameters["Brick Lenght Min"] = -Length_WG/2
        Parameters["Brick Lenght Max"] = Length_WG/2
        Parameters["Brick Width Min"] = -WidthObject/2 + Width_Electrodes + Width_WG_New + Gap*2
        Parameters["Brick Width Max"] = -WidthObject/2 + Width_Electrodes + Width_WG_New + Gap*2 + Width_Signal
        Parameters["Brick Hight Min"] = HeightZ
        Parameters["Brick Hight Max"] = HeightZ + Height_Electrodes
        Parameters["Brick Name"] = NamesElectrodes[2]
        Parameters["Component Name"] = NamesElectrodes[2]
        Parameters["Material"] = "Au"
        self.Brick(Parameters)



        Parameters = {}
        Parameters["Brick Lenght Min"] = -Length_WG/2
        Parameters["Brick Lenght Max"] = Length_WG/2
        Parameters["Brick Width Min"] = WidthObject/2 
        Parameters["Brick Width Max"] = WidthObject/2 - Width_Electrodes 
        Parameters["Brick Hight Min"] = HeightZ
        Parameters["Brick Hight Max"] = HeightZ + Height_Electrodes
        Parameters["Brick Name"] = NamesElectrodes[1]
        Parameters["Component Name"] = NamesElectrodes[1]
        Parameters["Material"] = "Au"
        self.Brick(Parameters)

        

        # Plase Waveguides
        GND_Left_Corner = (-WidthObject/2 + Width_Electrodes) + Gap + Width_WG_New/2
        GND_Right_Corner = (WidthObject/2 - Width_Electrodes) - Gap - Width_WG_New/2

        PosLeft = [round((GND_Left_Corner  + Width_WG_New/2),2), round((GND_Left_Corner  - Width_WG_New/2),2)] 
        PosRight = [round((GND_Right_Corner - Width_WG_New/2),2), round((GND_Right_Corner + Width_WG_New/2),2)] 

        # New Waveguide Placement with 2D Polygon and Translation 
        PointsLeft = {}
        PointsLeft["X"] = [Length_WG/2, Length_WG/2, Length_WG/2 - Height_WG, Length_WG/2 - Height_WG, Length_WG/2]
        PointsLeft["Y"] = [PosLeft[0]/2, PosLeft[1]/2, round((PosLeft[1] + extention),2)/2 , round((PosLeft[0] - extention),2)/2, PosLeft[0]/2]

        # Waveguide and Waveguide to 
        # Translation Parameters
        Trans = {}
        Trans["Translate Type"] = "Translate"
        Trans["Angle X"] = 0
        Trans["Angle Y"] = 90
        Trans["Angle Z"] = 0
        Trans["Position X"] = Height_WG/2
        Trans["Position Y"] = 0 
        Trans["Position Z"] = HeightZ/2 + Height_WG/2

        self.Poligon_2D(WGName = "Waveguide_Left", Points = PointsLeft)
        Trans["Name Object"] = "Waveguide_Left"
        self.RotationTranslation(Trans)
        self.Translation(Trans)
        self.RibWaveguide_ToSolid("Waveguide_Left", WaveguideName = "Waveguide_Left", WG_Hight = Length_WG, Angle = 0, WGFolderName = "Waveguide_Left", WGName = "Waveguide_Left", Material="LiNbO3")
       



        PointsRight = {}
        PointsRight["X"] = [Length_WG/2, Length_WG/2, Length_WG/2 - Height_WG, Length_WG/2 - Height_WG, Length_WG/2]
        PointsRight["Y"] = [PosRight[0]/2, PosRight[1]/2, round((PosRight[1] - extention),2)/2 , round((PosRight[0] + extention),2)/2, PosRight[0]/2]


        # Waveguide and Waveguide to solid
        self.Poligon_2D(WGName = "Waveguide_Right", Points = PointsRight)
        Trans["Name Object"] = "Waveguide_Right"
        self.RotationTranslation(Trans)
        self.Translation(Trans)
        self.RibWaveguide_ToSolid("Waveguide_Right", WaveguideName = "Waveguide_Right", WG_Hight = -Length_WG, Angle = 0, WGFolderName = "Waveguide_Right", WGName = "Waveguide_Right", Material="LiNbO3")





    # Phase Modulator Design
    def PhaseModulator(self, Parameters):
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

        self.add_anisotropy_material("LiNbO3", Data)
        self.add_SiO2()
        self.add_Au()


        # Calculate Top and Bottom Widths of the Waveguide
        x = abs(Height_WG / (np.cos((Angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - Height_WG ** 2)
        Width_WG_New = Width_WG + 2 * extention 

        # Global Parameters
        WidthObject = 2*Width_WG_New + Width_Signal + Width_Electrodes + 2*Gap


        # Substrate Definition 
        Parameters = {}
        Parameters["Brick Lenght Min"] = -Length_WG/2
        Parameters["Brick Lenght Max"] = Length_WG/2
        Parameters["Brick Width Min"] = -(WidthObject+2)/2
        Parameters["Brick Width Max"] = (WidthObject+2)/2
        Parameters["Brick Hight Min"] = -Height_Substrate/2
        Parameters["Brick Hight Max"] = Height_Substrate/2
        Parameters["Brick Name"] = "LNOI_Substrate"
        Parameters["Component Name"] = "LNOI_Substrate"
        Parameters["Material"] = "SiO2"
        self.Brick(Parameters)




        # Slab definition
        if Heigh_Slab == 0:
            HeightZ = Height_Substrate/2
        else:
            Parameters = {}
            Parameters["Brick Lenght Min"] = -Length_WG/2
            Parameters["Brick Lenght Max"] = Length_WG/2
            Parameters["Brick Width Min"] = -(WidthObject+2)/2
            Parameters["Brick Width Max"] = (WidthObject+2)/2
            Parameters["Brick Hight Min"] = Height_Substrate/2
            Parameters["Brick Hight Max"] = Height_Substrate/2 + Heigh_Slab
            Parameters["Brick Name"] = "LNOI_Slab"
            Parameters["Component Name"] = "LNOI_Slab"
            Parameters["Material"] = "LiNbO3"
            self.Brick(Parameters)
            HeightZ = Height_Substrate/2 + Heigh_Slab

            
        # Calculate Top and Bottom Widths of the Waveguide
        x = abs(Height_WG / (np.cos((Angle) * np.pi / 180)))  # in Radians
        extention = np.sqrt(x ** 2 - Height_WG ** 2)
        Width_WG_New = Width_WG + 2 * extention 

        #  Electrodes 
        NamesElectrodes = "Electrode", "Signal"

        # Plase Electrodes
        Parameters = {}
        Parameters["Brick Lenght Min"] = -Length_MZM/2
        Parameters["Brick Lenght Max"] = Length_MZM/2
        Parameters["Brick Width Min"] = -WidthObject/2 
        Parameters["Brick Width Max"] = -WidthObject/2 + Width_Electrodes
        Parameters["Brick Hight Min"] = HeightZ
        Parameters["Brick Hight Max"] = HeightZ + Height_Electrodes
        Parameters["Brick Name"] = NamesElectrodes[0]
        Parameters["Component Name"] = NamesElectrodes[0]
        Parameters["Material"] = "Au"
        self.Brick(Parameters)


        Parameters = {}
        Parameters["Brick Lenght Min"] = -Length_MZM/2
        Parameters["Brick Lenght Max"] = Length_MZM/2
        Parameters["Brick Width Min"] = -WidthObject/2 + Width_Electrodes + Width_WG_New + Gap*2
        Parameters["Brick Width Max"] = -WidthObject/2 + Width_Electrodes + Width_WG_New + Gap*2 + Width_Signal
        Parameters["Brick Hight Min"] = HeightZ
        Parameters["Brick Hight Max"] = HeightZ + Height_Electrodes
        Parameters["Brick Name"] = NamesElectrodes[1]
        Parameters["Component Name"] = NamesElectrodes[1]
        Parameters["Material"] = "Au"
        self.Brick(Parameters)


        # Plase Waveguides
        GND_Left_Corner = (-WidthObject/2 + Width_Electrodes) + Gap + Width_WG_New/2
        
        PosLeft = [round((GND_Left_Corner  + Width_WG_New/2),2), round((GND_Left_Corner  - Width_WG_New/2),2)] 
    

        # New Waveguide Placement with 2D Polygon and Translation 
        PointsLeft = {}
        PointsLeft["X"] = [Length_WG/2, Length_WG/2, Length_WG/2 - Height_WG, Length_WG/2 - Height_WG, Length_WG/2]
        PointsLeft["Y"] = [PosLeft[0]/2, PosLeft[1]/2, round((PosLeft[1] + extention),2)/2 , round((PosLeft[0] - extention),2)/2, PosLeft[0]/2]

        # Waveguide and Waveguide to 
        # Translation Parameters
        Trans = {}
        Trans["Translate Type"] = "Translate"
        Trans["Angle X"] = 0
        Trans["Angle Y"] = 90
        Trans["Angle Z"] = 0
        Trans["Position X"] = Height_WG/2
        Trans["Position Y"] = 0 
        Trans["Position Z"] = HeightZ/2 + Height_WG/2

        self.Poligon_2D(WGName = "Waveguide_Left", Points = PointsLeft)
        Trans["Name Object"] = "Waveguide_Left"
        self.RotationTranslation(Trans)
        self.Translation(Trans)
        self.RibWaveguide_ToSolid("Waveguide_Left", WaveguideName = "Waveguide_Left", WG_Hight = Length_WG, Angle = 0, WGFolderName = "Waveguide_Left", WGName = "Waveguide_Left", Material="LiNbO3")




    def Squere_Waveguide(self, Parameters):
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

        self.add_anisotropy_material("LiNbO3", Data)
        self.add_SiO2()
        self.add_Au()

        

        # Global Parameters for the WG
        WidthObject = 4*Width_WG 
    

        # Substrate Definition 
        Parameters = {}
        Parameters["Brick Lenght Min"] = -Length_WG/2
        Parameters["Brick Lenght Max"] = Length_WG/2
        Parameters["Brick Width Min"] = -WidthObject/2
        Parameters["Brick Width Max"] = WidthObject/2 
        Parameters["Brick Hight Min"] = -Height_Substrate/2
        Parameters["Brick Hight Max"] = Height_Substrate/2
        Parameters["Brick Name"] = "LNOI_Substrate"
        Parameters["Component Name"] = "LNOI_Substrate"
        Parameters["Material"] = "SiO2"
        self.Brick(Parameters)



        # Slab definition
        if Heigh_Slab == 0:
            HeightZ = Height_Substrate/2
        else:
            Parameters = {}
            Parameters["Brick Lenght Min"] = -Length_WG/2
            Parameters["Brick Lenght Max"] = Length_WG/2
            Parameters["Brick Width Min"] = -WidthObject/2
            Parameters["Brick Width Max"] = WidthObject/2 
            Parameters["Brick Hight Min"] = Height_Substrate/2
            Parameters["Brick Hight Max"] = Height_Substrate/2 + Heigh_Slab
            Parameters["Brick Name"] = "LNOI_Slab"
            Parameters["Component Name"] = "LNOI_Slab"
            Parameters["Material"] = "LiNbO3"
            self.Brick(Parameters)
            HeightZ = Height_Substrate/2 + Heigh_Slab
        
        # Create squere Waveguide
        Parameters = {}
        Parameters["Brick Lenght Min"] = -Length_WG/2
        Parameters["Brick Lenght Max"] = Length_WG/2
        Parameters["Brick Width Min"] = -Width_WG/2 
        Parameters["Brick Width Max"] = Width_WG/2 
        Parameters["Brick Hight Min"] = HeightZ
        Parameters["Brick Hight Max"] = HeightZ + Hight_WG
        Parameters["Brick Name"] = "WG"
        Parameters["Component Name"] = "WG"
        Parameters["Material"] = "LiNbO3"
        self.Brick(Parameters)





    def BondWire(self, NameWire, Coordinates, Height, Radius, BondwireType = "Spline", CenterPosition = 0.5, alpha = None, beta = None, Material = None, SolidWireModel = True, Termination = None, NameFolder = None):
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
            vba_code = f"""
                        With Wire
                        .Reset
                        .Folder "{NameFolder}"
                        .Type "Bondwire"
                        .Name "{NameWire}"
                        .BondWireType "{BondwireType}"
                        .Point1 "{Coordinates["X1"]}", "{Coordinates["Y1"]}", "{Coordinates["Z1"]}", "False"
                        .Point2 "{Coordinates["X2"]}", "{Coordinates["Y2"]}", "{Coordinates["Z2"]}", "False"
                        .Height "{Height}"
                        .Radius "{Radius}"
                        .RelativeCenterPosition "{CenterPosition}"
                        .Material "{Material}"
                        .SolidWireModel "{SolidWireModel}"
                        .Termination "{Termination}"
                        .Add
                        End With
                        """
            self.prj.model3d.add_to_history(f"create bondwire {NameWire}", vba_code)
            

        else:
            vba_code = f"""
                    With Wire
                    .Reset
                        .Folder "{NameFolder}"
                        .Type "Bondwire"
                        .Name "{NameWire}"
                        .BondWireType "{BondwireType}"
                        .Point1 "{Coordinates["X1"]}", "{Coordinates["Y1"]}", "{Coordinates["Z1"]}", "False"
                        .Point2 "{Coordinates["X2"]}", "{Coordinates["Y2"]}", "{Coordinates["Z2"]}", "False"
                        .Height "{Height}"
                        .Radius "{Radius}"
                        .alpha "{alpha}"
                        .beta "{beta}"
                        .Material "{Material}"
                        .SolidWireModel "{SolidWireModel}"
                        .Termination "{Termination}"
                        .Add
                        End With
                        """
            self.prj.model3d.add_to_history(f"create bondwire {NameWire}", vba_code)




    def CurveWire(self, NameWire, Radius, Points = None, Material = None, SolidWireModel = True, Termination = None, NameFolder = None, CurveFolderName = None, CurveName = None):
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


        vba_code = f"""
                    With Wire
                    .Reset
                    .Folder "{NameFolder}"
                    .Type "Curvewire"
                    .Name "{NameWire}"
                    .Curve "{CurveFolderName}:{CurveName}"
                    .Radius "{Radius}"
                    .SolidWireModel "{SolidWireModel}"
                    .Material "{Material}"
                    .Termination "{Termination}"
                    .Add
                    End With
                    """
        self.prj.model3d.add_to_history("create curve wire", vba_code)



    def Cylinder(self, Parameters):

        Name = Parameters["Cylinder Name"]
        Component_name = Parameters["Component Name"]
        Material = Parameters["Material"]
        Outer_Radius = Parameters["Outer Radius"]
        Inner_Radius = Parameters["Inner Radius"]
        Orentation_Axis = Parameters["Orentation Axis"]
        
        
        vba_code = f"""
                    With Cylinder 
                    .Reset 
                    .Name "{Name}" 
                    .Component "{Component_name}" 
                    .Material "{Material}" 
                    .OuterRadius "{Outer_Radius}" 
                    .InnerRadius "{Inner_Radius}" 
                    .Axis "{Orentation_Axis}" 
                    """
        if Orentation_Axis == "X":
            X_min = Parameters["X min"]
            X_max = Parameters["X max"]
            Zcenter = Parameters["Z center"]
            Ycenter = Parameters["Y center"]
            vba_code1 = f"""
                    .Xrange "{X_min}", "{X_max}" 
                    .Zcenter "{Zcenter}" 
                    .Ycenter "{Ycenter}" 
                    .Segments "0" 
                    .Create 
                End With
                """
            full_vba_code = vba_code + vba_code1

        elif Orentation_Axis == "Y":
            Y_min = Parameters["Y min"]
            Y_max = Parameters["Y max"]
            Zcenter = Parameters["Z center"]
            Xcenter = Parameters["X center"]
            vba_code1 = f"""
                    .Yrange "{Y_min}", "{Y_max}" 
                    .Zcenter "{Zcenter}" 
                    .Xcenter "{Xcenter}" 
                    .Segments "0" 
                    .Create 
                End With
                """
            full_vba_code = vba_code + vba_code1
        elif Orentation_Axis == "Z":
            Z_min = Parameters["Z min"]
            Z_max = Parameters["Z max"]
            Xcenter = Parameters["X center"]
            Ycenter = Parameters["Y center"]
            vba_code1 = f"""
                    .Zrange "{Z_min}", "{Z_max}" 
                    .Xcenter "{Xcenter}" 
                    .Ycenter "{Ycenter}" 
                    .Segments "0" 
                    .Create 
                End With
                """
            full_vba_code = vba_code + vba_code1
        
        else:
            raise ValueError("Orentation Axis can be only X, Y and Z!")
        
        self.prj.model3d.add_to_history("create cylinder", full_vba_code)




    def Conus(self, Parameters):
        
        Name = Parameters["Cylinder Name"]
        Component_name = Parameters["Component Name"]
        Material = Parameters["Material"]
        Outer_Radius = Parameters["Top Radius"]
        Inner_Radius = Parameters["Bottom Radius"]
        Orentation_Axis = Parameters["Orentation Axis"]
        vba_code = f"""
                    With Cone 
                     .Reset 
                     .Name "{Name}" 
                     .Component "{Component_name}" 
                     .Material "{Material}" 
                     .TopRadius "{Outer_Radius}" 
                     .BottomRadius "{Inner_Radius}" 
                     .Axis "{Orentation_Axis}"
                     """
        if Orentation_Axis == "X":
            X_min = Parameters["X min"]
            X_max = Parameters["X max"]
            Zcenter = Parameters["Z center"]
            Ycenter = Parameters["Y center"]
            vba_code1 = f"""
                    
                    .Xrange "{X_min}", "{X_max}" 
                    .Zcenter "{Zcenter}" 
                    .Ycenter "{Ycenter}" 
                    .Segments "0" 
                    .Create 
                End With
                """
            full_vba_code = vba_code + vba_code1

        elif Orentation_Axis == "Y":
            Y_min = Parameters["Y min"]
            Y_max = Parameters["Y max"]
            Zcenter = Parameters["Z center"]
            Xcenter = Parameters["X center"]
            vba_code1 = f"""
                    .Yrange "{Y_min}", "{Y_max}" 
                    .Zcenter "{Zcenter}" 
                    .Xcenter "{Xcenter}" 
                    .Segments "0" 
                    .Create 
                End With
                """
            full_vba_code = vba_code + vba_code1
            
        elif Orentation_Axis == "Z":
            Z_min = Parameters["Z min"]
            Z_max = Parameters["Z max"]
            Xcenter = Parameters["X center"]
            Ycenter = Parameters["Y center"]
            vba_code1 = f"""
                    .Zrange "{Z_min}", "{Z_max}" 
                    .Xcenter "{Xcenter}" 
                    .Ycenter "{Ycenter}" 
                    .Segments "0" 
                    .Create 
                End With
                """
            full_vba_code = vba_code + vba_code1
            
        self.prj.model3d.add_to_history("create cylinder", full_vba_code)




    def GGB_Probe(self, Parameters):
        
        
        # Global Parameters
        ComponentName = Parameters["Component Name"]
        angleOrientation = Parameters["Orientation Angle"]
        self.add_GGB_Probe_Material()
       
        # Create GGB Probes Body Cylinder
        DictComponent = {}
        DictComponent["Cylinder Name"] = "Body_Probe"
        DictComponent["Component Name"] = ComponentName
        DictComponent["Material"] = "PEC"
        DictComponent["Outer Radius"] =  290
        DictComponent["Inner Radius"] =  170
        DictComponent["Z min"] = 0
        DictComponent["Z max"] = 900
        DictComponent["X center"] = 0   
        DictComponent["Y center"] = 0
        DictComponent["Orentation Axis"] = "Z"
        
        self.Cylinder(DictComponent)
        
        
        # Create GGB Probes Signal Cylinder
        DictComponent = {}
        DictComponent["Cylinder Name"] = "Signal_Probe"
        DictComponent["Component Name"] = ComponentName
        DictComponent["Material"] = "PEC"
        DictComponent["Outer Radius"] =  50
        DictComponent["Inner Radius"] =  0
        DictComponent["Z min"] = 0
        DictComponent["Z max"] = 900
        DictComponent["X center"] = 0   
        DictComponent["Y center"] = 0
        DictComponent["Orentation Axis"] = "Z"
        
        self.Cylinder(DictComponent)
        
        
        # Create GGB Probes Isolation Cylinder
        DictComponent = {}
        DictComponent["Cylinder Name"] = "Isolation_Probe"
        DictComponent["Component Name"] = ComponentName
        DictComponent["Material"] = "Teflon (PTFE) (loss free)"
        DictComponent["Outer Radius"] =  170
        DictComponent["Inner Radius"] = 50
        DictComponent["Z min"] = 0
        DictComponent["Z max"] = 900
        DictComponent["X center"] = 0   
        DictComponent["Y center"] = 0
        DictComponent["Orentation Axis"] = "Z"
        
        self.Cylinder(DictComponent)
        
        
        
        
        # Create Signal Conus
        DictComponent = {}
        DictComponent["Cylinder Name"] = "Conus_Probe_Sig"
        DictComponent["Component Name"] = ComponentName
        DictComponent["Material"] = "PEC"
        DictComponent["Top Radius"] = 50
        DictComponent["Bottom Radius"] = 51.91/2
        DictComponent["Z min"] = -197.30
        DictComponent["Z max"] = 0
        DictComponent["X center"] = 0   
        DictComponent["Y center"] = 0
        DictComponent["Orentation Axis"] = "Z"
        
        self.Conus(DictComponent)
        
        
        # Create Signal Conus Tip
        DictComponent = {}
        DictComponent["Cylinder Name"] = "Conus_Probe_Sig_Tip"
        DictComponent["Component Name"] = ComponentName
        DictComponent["Material"] = "PEC"
        DictComponent["Top Radius"] = 51.91/2
        DictComponent["Bottom Radius"] = 10.31/2
        DictComponent["Z min"] = -197.30 - 45.10
        DictComponent["Z max"] = -197.30
        DictComponent["X center"] = 0   
        DictComponent["Y center"] = 0
        DictComponent["Orentation Axis"] = "Z"
        
        self.Conus(DictComponent)
        
        
        # Create Ground forms Left
        DictComponent = {}
        DictComponent["Name"] = "GND_Curve_Pin"
        DictComponent["Curve Name"] = "GND_Curve_Pin"
        Points = {}
        Points['X'] = [0, 0, 0, 0, 0, 0]
        Points['Y'] = [105, 130, 85, -5, -105, 105]
        Points['Z'] = [0, -185, -240, -200, 0, 0 ]   
        DictComponent["Point"] = Points
        
        self.Poligon_3D(DictComponent)
        
        
        
        DictComponent = {}
        DictComponent["Name"] = ComponentName
        DictComponent["Component Name"] = "Porbe_GND_L_Tip"
        DictComponent["Material"] = "PEC"
        DictComponent["Thickness"] = 25
        DictComponent["Curve Folder"] = "GND_Curve_Pin"
        DictComponent["Curve Name"] = "GND_Curve_Pin"
        
        
        self.CurveToSolid(DictComponent)
        
        
        # Translate to right 
        DictComponent = {}
        DictComponent["Translate Type"] = "Translate"
        DictComponent["Name Object"] = f"{ComponentName}:Porbe_GND_L_Tip"
        DictComponent["Position X"] = 25/2
        DictComponent["Position Y"] = -186.1
        DictComponent["Position Z"] = 0
        DictComponent["Structure Type"] = "Shape"
        
        
        self.Translation(DictComponent)
        
        
        
        
        
        #  Create Ground forms Right
        DictComponent = {}
        DictComponent["Name"] = "GND_Curve_Pin"
        DictComponent["Curve Name"] = "GND_Curve_Pin"
        Points = {}
        Points['X'] = [0, 0, 0, 0, 0, 0]
        Points['Y'] = [-105, -130, -85, 5, 105, -105]
        Points['Z'] = [0, -185, -240, -200, 0, 0 ]   
        DictComponent["Point"] = Points
        
        self.Poligon_3D(DictComponent)
        
        
        
        DictComponent = {}
        DictComponent["Name"] = ComponentName
        DictComponent["Component Name"] = "Porbe_GND_R_Tip"
        DictComponent["Material"] = "PEC"
        DictComponent["Thickness"] = 25
        DictComponent["Curve Folder"] = "GND_Curve_Pin"
        DictComponent["Curve Name"] = "GND_Curve_Pin"
        
        self.CurveToSolid(DictComponent)
        
        
        # Translate to right 180 degree
        DictComponent = {}
        DictComponent["Translate Type"] = "Translate"
        DictComponent["Name Object"] = f"{ComponentName}:Porbe_GND_R_Tip"
        DictComponent["Position X"] = -25/2
        DictComponent["Position Y"] = 186.1
        DictComponent["Position Z"] = 0
        DictComponent["Structure Type"] = "Shape"
        
        self.Translation(DictComponent)
        
        
        # Create Sphere
        DictComponent = {}
        DictComponent["Name"] = "Elipse"
        DictComponent["Component Name"] = ComponentName
        DictComponent["Material"] =  "epoxy_casting_CR110"
        DictComponent["Axis"] = "Z"
        DictComponent["Center Radius"] = 100
        Pos = {}
        Pos["X"] = 0
        Pos["Y"] = 0
        Pos["Z"] = 0
        DictComponent["Center Positions"] = Pos
        
        self.Sphere(DictComponent)
        
        
        #Transform Sphere to Elipse 
        DictComponent = {}
        DictComponent["Translate Type"] = "Scale" 
        DictComponent["Name Object"] = f"{ComponentName}:Elipse" 
        Pos = {}
        Pos["X"] = 0
        Pos["Y"] = 0
        Pos["Z"] = 0
        DictComponent["Center Positions"] = Pos
        Scale = {}
        Scale["X"] = 1.5
        Scale["Y"] = 2
        Scale["Z"] = 1
        DictComponent["Scale Factors"] = Scale
        
        self.Translation(DictComponent)
        
        
        # Subtract
        DictComponent = {}
        DictComponent["Name Cut Structure"] = f"{ComponentName}:Elipse"
        Structures = [f"{ComponentName}:Conus_Probe_Sig" , f"{ComponentName}:Isolation_Probe", f"{ComponentName}:Porbe_GND_L_Tip", f"{ComponentName}:Porbe_GND_R_Tip", f"{ComponentName}:Body_Probe", f"{ComponentName}:Signal_Probe"]
        
        
        for i in range(len(Structures)):
            DictComponent["Name Structure to Cut"] = Structures[i]    
            self.Cut_structures(DictComponent)
            
            
        # Rebuild the missing parts after the sphere 
        # Create GGB Probes Body Cylinder    
        DictComponent = {}
        DictComponent["Cylinder Name"] = "Body_Probe"
        DictComponent["Component Name"] = ComponentName
        DictComponent["Material"] = "PEC"
        DictComponent["Outer Radius"] =  290
        DictComponent["Inner Radius"] =  170
        DictComponent["Z min"] = 0
        DictComponent["Z max"] = 900
        DictComponent["X center"] = 0   
        DictComponent["Y center"] = 0
        DictComponent["Orentation Axis"] = "Z"
        
        self.Cylinder(DictComponent)
        
        
        
        
        # Create GGB Probes Signal Cylinder
        DictComponent = {}
        DictComponent["Cylinder Name"] = "Signal_Probe"
        DictComponent["Component Name"] = ComponentName
        DictComponent["Material"] = "PEC"
        DictComponent["Outer Radius"] =  50
        DictComponent["Inner Radius"] =  0
        DictComponent["Z min"] = 0
        DictComponent["Z max"] = 900
        DictComponent["X center"] = 0   
        DictComponent["Y center"] = 0
        DictComponent["Orentation Axis"] = "Z"
        
        self.Cylinder(DictComponent)

        
        # Create GGB Probes Isolation Cylinder
        DictComponent = {}
        DictComponent["Cylinder Name"] = "Isolation_Probe"
        DictComponent["Component Name"] = ComponentName
        DictComponent["Material"] = "Teflon (PTFE) (loss free)"
        DictComponent["Outer Radius"] =  170
        DictComponent["Inner Radius"] = 50
        DictComponent["Z min"] = 0
        DictComponent["Z max"] = 900
        DictComponent["X center"] = 0   
        DictComponent["Y center"] = 0
        DictComponent["Orentation Axis"] = "Z"
        
        self.Cylinder(DictComponent)
        
   
        # Create Signal Conus
        DictComponent = {}
        DictComponent["Cylinder Name"] = "Conus_Probe_Sig"
        DictComponent["Component Name"] = ComponentName
        DictComponent["Material"] = "PEC"
        DictComponent["Top Radius"] = 50
        DictComponent["Bottom Radius"] = 51.91/2
        DictComponent["Z min"] = -197.30
        DictComponent["Z max"] = 0
        DictComponent["X center"] = 0   
        DictComponent["Y center"] = 0
        DictComponent["Orentation Axis"] = "Z"
        
        self.Conus(DictComponent)
        

        # Create Ground forms Left
        DictComponent = {}
        DictComponent["Name"] = "GND_Curve_Pin"
        DictComponent["Curve Name"] = "GND_Curve_Pin"  
        Points['X'] = [0, 0, 0, 0, 0, 0]
        Points['Y'] = [105, 130, 85, -5, -105, 105]
        Points['Z'] = [0, -185, -240, -200, 0, 0 ]   
        DictComponent["Point"] = Points
        
        self.Poligon_3D(DictComponent)
        
        
        
        DictComponent = {}
        DictComponent["Name"] = ComponentName
        DictComponent["Component Name"] = "Porbe_GND_L_Tip"
        DictComponent["Material"] = "PEC"
        DictComponent["Thickness"] = 25
        DictComponent["Curve Folder"] = "GND_Curve_Pin"
        DictComponent["Curve Name"] = "GND_Curve_Pin"
        
        self.CurveToSolid( DictComponent)
        
        
        # Translate to right 
        DictComponent = {}
        DictComponent["Translate Type"] = "Translate"
        DictComponent["Name Object"] = f"{ComponentName}:Porbe_GND_L_Tip"
        DictComponent["Position X"] = 25/2
        DictComponent["Position Y"] = -186.1
        DictComponent["Position Z"] = 0
        DictComponent["Structure Type"] = "Shape"
        
        self.Translation(DictComponent)
        
        
        #  Create Ground forms Right
        DictComponent = {}
        
        DictComponent["Name"] = "GND_Curve_Pin"
        DictComponent["Curve Name"] = "GND_Curve_Pin"
        
        Points['X'] = [0, 0, 0, 0, 0, 0]
        Points['Y'] = [-105, -130, -85, 5, 105, -105]
        Points['Z'] = [0, -185, -240, -200, 0, 0 ]   
        DictComponent["Point"] = Points
        
        self.Poligon_3D(DictComponent)
        
        
        DictComponent = {}
        DictComponent["Name"] = ComponentName
        DictComponent["Component Name"] = "Porbe_GND_R_Tip"
        DictComponent["Material"] = "PEC"
        DictComponent["Thickness"] = 25
        DictComponent["Curve Folder"] = "GND_Curve_Pin"
        DictComponent["Curve Name"] = "GND_Curve_Pin"
        
        
        self.CurveToSolid(DictComponent)
        
        
        # Translate to right 180 degree
        DictComponent = {}
        DictComponent["Translate Type"] = "Translate"
        DictComponent["Name Object"] = f"{ComponentName}:Porbe_GND_R_Tip"
        DictComponent["Position X"] = -25/2
        DictComponent["Position Y"] = 186.1
        DictComponent["Position Z"] = 0
        DictComponent["Structure Type"] = "Shape"
        
        self.Translation(DictComponent)
        
        
        # Translate probe to -30 degree
        DictComponent = {}
        DictComponent["Translate Type"] = "Rotate"
        DictComponent["Name Object"] = ComponentName
        Pos["X"] = 0
        Pos["Y"] = 0
        Pos["Z"] = 0
        DictComponent["Center Positions"] = Pos
        angle = {}
        angle["X"] = 0
        angle["Y"] = angleOrientation
        angle["Z"] = 0
        DictComponent["Angle Values"] = angle
        
        self.Translation(DictComponent)
        
        # Subtract
        ObjectsToCut = [f"{ComponentName}:Conus_Probe_Sig_Tip", f"{ComponentName}:Porbe_GND_L_Tip", f"{ComponentName}:Porbe_GND_R_Tip" , f"{ComponentName}:Elipse",  f"{ComponentName}:Body_Probe", f"{ComponentName}:Isolation_Probe"]
        if angleOrientation<0:
            angleCutPlate = 30
        else:
            angleCutPlate = -30
            
        for i in range(len(ObjectsToCut)):
                
            # 30 Degree Plate for the cut 
            DictComponent["Brick Lenght Max"] = 1000 
            DictComponent["Brick Lenght Min"] = -1000
            DictComponent["Brick Width Max"] = 750 
            DictComponent["Brick Width Min"] = -750 
            DictComponent["Brick Hight Max"] = 0
            DictComponent["Brick Hight Min"] = -300
            DictComponent["Brick Name"] = "Cut_Plate"
            DictComponent["Component Name"] = ComponentName
            DictComponent["Material"] = "PEC"
            self.Brick(DictComponent)
            
            
            
            DictComponent = {}
            DictComponent["Translate Type"] = "Rotate"
            DictComponent["Name Object"] = f"{ComponentName}:Cut_Plate"
            Pos["X"] = 0
            Pos["Y"] = 0
            Pos["Z"] = 0
            DictComponent["Center Positions"] = Pos
            angle = {}
            angle["X"] = 0
            angle["Y"] = angleCutPlate
            angle["Z"] = 0
            DictComponent["Angle Values"] = angle
            
            
            self.Translation(DictComponent)
            

            DictComponent = {}
            DictComponent["Translate Type"] = "Translate"
            DictComponent["Name Object"] = f"{ComponentName}:Cut_Plate"
            DictComponent["Position X"] = 0
            DictComponent["Position Y"] = 0
            DictComponent["Position Z"] = -144
            DictComponent["Structure Type"] = "Shape"
            
            self.Translation(DictComponent)
            
            DictComponent = {}
            DictComponent["Name Structure to Cut"]  = f"{ComponentName}:Cut_Plate"
            DictComponent["Name Cut Structure"] = ObjectsToCut[i]
            self.Cut_structures(DictComponent)
            
            
        # Tranlate Probe to Tips Z = 0 Positiopn 
        DictComponent = {}
        DictComponent["Translate Type"] = "Rotate"
        DictComponent["Name Object"] = ComponentName
        Pos["X"] = 0
        Pos["Y"] = 0
        Pos["Z"] = 0
        DictComponent["Center Positions"] = Pos
        angle = {}
        angle["X"] = 0
        angle["Y"] = angleOrientation
        angle["Z"] = 0
        DictComponent["Angle Values"] = angle
        
        self.Translation(DictComponent)
        
        DictComponent = {}
        DictComponent["Translate Type"] = "Translate"
        DictComponent["Name Object"] = ComponentName
        DictComponent["Position X"] = 0
        DictComponent["Position Y"] = 0
        DictComponent["Position Z"] = 124.5
        
        self.Translation(DictComponent)
        
       




############################################################################
# Transform to Solid
############################################################################

    def ToSolid(self, SolidName, CurveName = "Polygon", NameFolder = None, Material = None ):
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


        vba_code = f"""
                    With Wire
                    .Reset
                    .SolidName "{SolidName}:{SolidName}"
                    .Name "{CurveName}"
                    .Folder "{NameFolder}"
                    .Material "{Material}"
                    .KeepWire "False"
                    .ConvertToSolidShape
                    End With
                    """
        self.prj.model3d.add_to_history(f"to solid {CurveName}", vba_code)

    
    def CurveToSolid(self,  Parameters):
        
        Name = Parameters["Name"] 
        NameComponent = Parameters["Component Name"]
        Material = Parameters["Material"]
        Thickness = Parameters["Thickness"]
        CurveFolder = Parameters["Curve Folder"]
        CurveName = Parameters["Curve Name"]
        
        if "angle" not in Parameters.keys():
            Angle = 0
        else:
            Angle = Parameters["angle"]
        
        vba_code = f"""
                    With ExtrudeCurve
                    .Reset
                    .Name "{NameComponent}"
                    .Component "{Name}"
                    .Material "{Material}"
                    .Thickness "{Thickness}"
                    .Twistangle "0"
                    .Taperangle "{Angle}"
                    .Curve "{CurveFolder}:{CurveName}"
                    .Create
                    End With
                    """
        self.prj.model3d.add_to_history(f"Curve to solid {NameComponent}", vba_code)
                
        
        





    def RibWaveguide_ToSolid(self, SolidName, WaveguideName = "Rib_Waveguide", WG_Hight = None, Angle = None, NameFolder = None, Material = None, WGFolderName = None, WGName = None ):
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

        vba_code = f"""
                    With ExtrudeCurve
                    .Reset
                    .Name "{WaveguideName}"
                    .Component "{WaveguideName}"
                    .Material "{Material}"
                    .Thickness "{WG_Hight}"
                    .Twistangle "0"
                    .Taperangle "{Angle}"
                    .Curve "{WGFolderName}:{WGName}"
                    .Create
                    End With
                    """
        self.prj.model3d.add_to_history("create rib waveguide", vba_code) 



############################################################################
# Ports
############################################################################


    def WaveguidePort(self, Parameters):
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
            vba_code = f"""
                        With Port
                        .Reset
                        .PortNumber "{PortNum[index-1]}"
                        .NumberOfModes "5"\
                        .ReferencePlaneDistance "0"
                        .Coordinates "{Coordinates}"
                        .Orientation "{Orientation[index-1]}"
                        .PortOnBound "False"
                        .ClipPickedPortToBound "False"
                        .YrangeAdd  "{Span11[index-1]}", "{Span12[index-1]}"
                        .ZrangeAdd  "{Span21[index-1]}", "{Span22[index-1]}"
                        .AdjustPolarization "True"
                        .PolarizationAngle "0"
                        .Create
                        End With
                        """
            self.prj.model3d.add_to_history("create waveguide port", vba_code) 




    def WaveguidePortWithPins(self, Parameters):
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
            PicParams (dict): Dictionary with Pick Faces Parameters
                                PickParams["Face Number"] = Integer with the number of the picked face of the structure

        Returns:
            str: String with VBA Code 
        """

        # Parameters to determin port position 
        Orientation = Parameters["Orientation"]
        Coordinates = Parameters["Coordinates"]
        # Choose coordinates
        PortNumber = Parameters["Port Number"]
        Span = Parameters["Span"]
        set_number = Parameters["Picked Port Number"] 
        Number_of_picks = Parameters["Number of picks"]
        
        
        lines = []

        lines = []
        if Number_of_picks > 2:
            potential = Parameters["Picked Port Polarity"]
            marked_face = Parameters["Picked Component Name"]
            facet_Id = Parameters["Face Number"]
            for i in range(len(marked_face)):
                lines.append(f'.AddPotentialPicked "{set_number}", "{potential[i]}", "{marked_face[i]}", "{facet_Id[i]}"')

            # join them with newlines
            line_block = "\n".join(lines)
            

            vba_code = f"""
                        With Port
                            .Reset
                            .PortNumber "{PortNumber}"
                            .Label ""
                            .Folder ""
                            .NumberOfModes "1"
                            .AdjustPolarization "False"
                            .PolarizationAngle "0.0"
                            .ReferencePlaneDistance "0.0"
                            .TextSize "50"
                            .TextMaxLimit "0"
                            .Coordinates "{Coordinates}"
                            .Orientation "{Orientation}"
                            .PortOnBound "True"
                            .ClipPickedPortToBound "False"
                            .Xrange "0", "0"
                            .Yrange "0", "0"
                            .Zrange "0", "0"
                            .XrangeAdd "{Span[0]}", "{Span[0]}"
                            .YrangeAdd "{Span[1]}", "{Span[1]}"
                            .ZrangeAdd "{Span[2]}", "{Span[2]}" 
                            .SingleEnded "False"
                            {line_block}
                            .WaveguideMonitor "False"
                            .Create
                        End With
                        """
            self.prj.model3d.add_to_history(f"create waveguide port {PortNumber}", vba_code)


        else:
            vba_code = f"""
                        With Port
                        .Reset
                        .PortNumber "{PortNumber}"
                        .Label ""
                        .Folder ""
                        .NumberOfModes "1"
                        .AdjustPolarization "False"
                        .PolarizationAngle "0"
                        .ReferencePlaneDistance "0"
                        .TextSize "50"
                        .TextMaxLimit "1"
                        .Coordinates "{Coordinates}"
                        .Orientation "{Orientation} "
                        .PortOnBound "True"
                        .ClipPickedPortToBound "False"
                        .Xrange "0", "100"
                        .Yrange "0", "100"
                        .Zrange "0", "0"
                        .XrangeAdd "{Span[0]}", "{Span[0]}"
                        .YrangeAdd "{Span[1]}", "{Span[1]}"
                        .ZrangeAdd "{Span[2]}", "{Span[2]}"
                        .SingleEnded "False"
                        .WaveguideMonitor "False"
                        .Create
                        End With
                        """
            self.prj.model3d.add_to_history("create waveguide port", vba_code)

        
    def MoveWaveguidePorts(self, Parameters):
        Distance = Parameters["Distance"]
        Span11 = Parameters["Span"][0][0]
        Span12 = Parameters["Span"][0][1]
        Span21 = Parameters["Span"][1][0]
        Span22 = Parameters["Span"][1][1]

        # Ports Numbers
        PortNum = Parameters["Port Number"]

        Port = {}
        Port["1"] = None
        Port["2"] = None


        # Xposition of thw teo ends 
        Xpos = [-(0.2 + Distance/2), (0.2 + Distance/2)]
        #Start Loop for Ports
        for index, i in enumerate(PortNum,start = 1):
            vba_code = f"""
                        With Port
                        .Reset
                        LoadContentForModify "{PortNum[index-1]}"
                        .Coordinates "Free"
                        .ClipPickedPortToBound "False"
                        .Yrange  "{Span11[index-1]}", "{Span12[index-1]}"
                        .Zrange  "{Span21[index-1]}", "{Span22[index-1]}"
                        .Xrange "{Xpos[index - 1]}", "{Xpos[index - 1]}"
                        .NumberOfModes "5"
                        .AdjustPolarization "True"
                        .PolarizationAngle "0"
                        .Modify
                        End With
                        """
            self.prj.model3d.add_to_history("move waveguide port", vba_code)

    def SetDiscretePort(self, Parameters):

        PortNumber = Parameters["Discrete Port Number"]
        if Parameters["Discrete Port Type"] in ["Voltage", "S-Parameters", "Current"]:
            PortType = Parameters["Discrete Port Type"]
        else:
            raise ValueError("Port Type can be on of 'Voltage', 'S-Parameters' or 'Current'")
        PortImpedance = Parameters["Port Impedance"]
        Voltage = Parameters["Port Voltage"]
        Current = Parameters["Port Current"]
        PortRadius = Parameters["Port Radius"]


        # Port = {}
        # Port["1"] = None
        # Port["2"] = None

        if "Discrete Port Coordinates" in list(Parameters.keys()):
            Coordinates = Parameters["Discrete Port Coordinates"]
            vba_code = f"""
                    With DiscretePort
                    .Reset
                    .PortNumber "{PortNumber}"
                    .Type "{PortType}"
                        .Label ""
                        .Folder ""
                        .Impedance "{PortImpedance}"
                        .Voltage "{Voltage}"
                        .Current "{Current}"
                        .Monitor "True"
                        .Radius "{PortRadius}"
                        .SetP1 "False", "{Coordinates["X1"]}", "{Coordinates["Y1"]}", "{Coordinates["Z1"]}"
                        .SetP2 "False", "{Coordinates["X2"]}", "{Coordinates["Y2"]}", "{Coordinates["Z2"]}"
                        .InvertDirection "False"
                        .LocalCoordinates "False"
                        .Wire ""
                        .Position "end1"
                        .Create
                        End With
                        """

            self.prj.model3d.add_to_history("create discrete waveguide port", vba_code)

        else:
            vba_code = f"""
                With DiscretePort
                .Reset
                .PortNumber "{PortNumber}"
                .Type "{PortType}"
                .Label ""
                .Folder ""
                .Impedance "{PortImpedance}"
                .Impedance "{PortImpedance}"
                .Voltage "{Voltage}"
                .Current "{Current}"
                .Monitor "True"
                .Radius "{PortRadius}"
                .SetP1 "True", "-23", "-2.9618802153517", "1.7"
                .SetP2 "True", "-23", "3", "1.7"
                .InvertDirection "False"
                .LocalCoordinates "False"
                .Wire ""
                .Position "end1"
                .Create 
                End With
                """
            self.prj.model3d.add_to_history("move waveguide port", vba_code)
            



# def WaveguidePorts_on_Electrodes_MZM(self, Parameters, Obj):
#     for i in range(len(Parameters["Port Number"])):
#         for j in range(len(Parameters["Electrodes_Names"])):
#             PicParams = {}
#             PicParams["Option"] = "Face"
#             PicParams["Object"] = Parameters["Electrodes_Names"][j]
#             PicParams["Face Number"] = Parameters["Face ID"][i]
#             PickFace = Pick(PicParams)
#             Parameters["Solid Name"] = Parameters["Electrodes_Names"][j]
#             Obj.schematic.execute_vba_code(PickFace, timeout=None)
        
#         Port = WaveguidePortWithPins(Parameters, PicParams)
#         Obj.schematic.execute_vba_code(Port[str(i+1)], timeout=None)





# def Optical_WaveguidePorts_MZM(self, Parameters, Obj):
#     for i in range(len(Parameters["Port Number"])):
#         for j in range(len(Parameters["Electrodes_Names"])):
#             PicParams = {}
#             PicParams["Option"] = "Face"
#             PicParams["Object"] = Parameters["Electrodes_Names"][j]
#             PicParams["Face Number"] = Parameters["Face ID"][i]
#             PickFace = Pick(PicParams)
#             Parameters["Solid Name"] = Parameters["Electrodes_Names"][j]
#             Obj.schematic.execute_vba_code(PickFace, timeout=None)
        
#         Port = WaveguidePort(Parameters)
#         Obj.schematic.execute_vba_code(Port[str(i+1)], timeout=None)






# def Discrete_Port(self, Parameters, Obj):
#     for i in range(len(Parameters["Port Number"])):
#         Parameters["Discrete Port Number"] = Parameters["Port Number"][i]
#         for j in range(len(Parameters["Port Number"])):
#             PicParams = {}
#             PicParams["Option"] = "Centerpoint"
#             PicParams["Object"] = Parameters["Electrodes_Names"][j]
#             PicParams["Face Number"] = Parameters["Face ID"][i]
#             PickFace = Pick(PicParams)
#             Parameters["Solid Name"] = Parameters["Electrodes_Names"][j]
#             Obj.schematic.execute_vba_code(PickFace, timeout=None)

#         DiscretePort = SetDiscretePort(Parameters)
#         Obj.schematic.execute_vba_code(DiscretePort, timeout=None)
    
############################################################################
# Picks
############################################################################
                

    def Pick(self, Parameters):
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
            vba_code = f"""
                        Pick.PickFaceFromId "{FaceName}","{Number}"
                    """
            self.prj.model3d.add_to_history("set pick", vba_code)
        elif Parameters["Option"] == "Centerpoint":
            vba_code = f"""
                        Pick.PickCenterpointFromId "{FaceName}", "{Number}"
                    """
            self.prj.model3d.add_to_history("set pick", vba_code)
        else:
            raise ValueError("Wrong Object, can be choose from 'Face' and 'Centerpoint' for the moment!")

        
    def ClearAllPicks(self):
        """Clear all pciks

        Returns:
            str: String with VBA Code 
        """
        vba_code = f"""
                Pick.ClearAllPicks
                """ 
        self.prj.model3d.add_to_history("delete pick", vba_code)    





                
                    


############################################################################
# Solvers
############################################################################

    def setTimeSolver(self, Parameters):
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
        MeshType = Parameters["Solver Mesh Type"]

        if MeshType == "FIT":
            vba_code = f"""
                        With Solver
                        .Reset
                        .Method "Hexahedral"
                        .SteadyStateLimit "-{Accuracy}"
                        .CalculateModesOnly "{ModesOnly}"
                        .StimulationPort "{Source}"
                        .StimulationMode "All"
                        .MeshAdaption "False"
                        .SParaSymmetry "True"
                        .AutoNormImpedance  "{AutoImpedance}"
                        .NormingImpedance "{Impedance}"
                        End With
                        """
            self.prj.model3d.add_to_history("set time solver", vba_code)  

        else:
            vba_code = f"""
                        'Mesh.SetCreator "High Frequency"
                        'With Solver
                        '.Reset
                        '.Method "Hexahedral TLM"
                        '.SteadyStateLimit "-{Accuracy}"
                        '.StimulationPort "All"
                        '.StimulationMode "All"
                        '.MeshAdaption "False"
                        '.SParaSymmetry "True"
                        '.AutoNormImpedance "True"
                        '.NormingImpedance ' + '"{Impedance}"
                        '.StoreTDResultsInCache  "False"
                        '.RunDiscretizerOnly "False"
                        '.SuperimposePLWExcitation "False"
                        'End With
                        """
            self.prj.model3d.add_to_history("set time solver", vba_code) 






    def setFreqSolver(self, Parameters):
        """Set Frequency solver Parameters

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
        # Accuracy = Parameters["Accuracy"]
        # ModesOnly = Parameters["Caclculate Modes Only"]
        # AutoImpedance = Parameters["Auto Impedance"]
        # Impedance = Parameters["Impedance"]
        # Source = Parameters["Source Port"]

        vba_code = f"""
            Mesh.SetCreator "High Frequency" 

            With FDSolver
                .Reset 
                .SetMethod "Tetrahedral", "General purpose" 
                .OrderTet "Second" 
                .OrderSrf "First" 
                .Stimulation "All", "All" 
                .ResetExcitationList 
                .AutoNormImpedance "True" 
                .NormingImpedance "50" 
                .ModesOnly "False" 
                .ConsiderPortLossesTet "True" 
                .SetShieldAllPorts "False" 
                .AccuracyHex "1e-6" 
                .AccuracyTet "1e-4" 
                .AccuracySrf "1e-3" 
                .LimitIterations "False" 
                .MaxIterations "0" 
                .SetCalcBlockExcitationsInParallel "True", "True", "" 
                .StoreAllResults "False" 
                .StoreResultsInCache "False" 
                .UseHelmholtzEquation "True" 
                .LowFrequencyStabilization "True" 
                .Type "Auto" 
                .MeshAdaptionHex "False" 
                .MeshAdaptionTet "True" 
                .AcceleratedRestart "True" 
                .FreqDistAdaptMode "Distributed" 
                .NewIterativeSolver "True" 
                .TDCompatibleMaterials "False" 
                .ExtrudeOpenBC "False" 
                .SetOpenBCTypeHex "Default" 
                .SetOpenBCTypeTet "Default" 
                .AddMonitorSamples "True" 
                .CalcPowerLoss "True" 
                .CalcPowerLossPerComponent "False" 
                .SetKeepSolutionCoefficients "MonitorsAndMeshAdaptation" 
                .UseDoublePrecision "False" 
                .UseDoublePrecision_ML "True" 
                .MixedOrderSrf "False" 
                .MixedOrderTet "False" 
                .PreconditionerAccuracyIntEq "0.15" 
                .MLFMMAccuracy "Default" 
                .MinMLFMMBoxSize "0.3" 
                .UseCFIEForCPECIntEq "True" 
                .UseEnhancedCFIE2 "True" 
                .UseFastRCSSweepIntEq "false" 
                .UseSensitivityAnalysis "False" 
                .UseEnhancedNFSImprint "True" 
                .UseFastDirectFFCalc "True" 
                .RemoveAllStopCriteria "Hex"
                .AddStopCriterion "All S-Parameters", "0.01", "2", "Hex", "True"
                .AddStopCriterion "Reflection S-Parameters", "0.01", "2", "Hex", "False"
                .AddStopCriterion "Transmission S-Parameters", "0.01", "2", "Hex", "False"
                .RemoveAllStopCriteria "Tet"
                .AddStopCriterion "All S-Parameters", "0.01", "2", "Tet", "True"
                .AddStopCriterion "Reflection S-Parameters", "0.01", "2", "Tet", "False"
                .AddStopCriterion "Transmission S-Parameters", "0.01", "2", "Tet", "False"
                .AddStopCriterion "All Probes", "0.05", "2", "Tet", "True"
                .RemoveAllStopCriteria "Srf"
                .AddStopCriterion "All S-Parameters", "0.01", "2", "Srf", "True"
                .AddStopCriterion "Reflection S-Parameters", "0.01", "2", "Srf", "False"
                .AddStopCriterion "Transmission S-Parameters", "0.01", "2", "Srf", "False"
                .SweepMinimumSamples "3" 
                .SetNumberOfResultDataSamples "1001" 
                .SetResultDataSamplingMode "Automatic" 
                .SweepWeightEvanescent "1.0" 
                .AccuracyROM "1e-4" 
                .AddSampleInterval "", "", "1", "Automatic", "True" 
                .AddSampleInterval "", "", "", "Automatic", "False" 
                .MPIParallelization "False"
                .UseDistributedComputing "False"
                .NetworkComputingStrategy "RunRemote"
                .NetworkComputingJobCount "3"
                .UseParallelization "True"
                .MaxCPUs "1024"
                .MaximumNumberOfCPUDevices "2"
                .HardwareAcceleration "False"
                .MaximumNumberOfGPUs "1"
            End With

            With IESolver
                .Reset 
                .UseFastFrequencySweep "True" 
                .UseIEGroundPlane "False" 
                .SetRealGroundMaterialName "" 
                .CalcFarFieldInRealGround "False" 
                .RealGroundModelType "Auto" 
                .PreconditionerType "Auto" 
                .ExtendThinWireModelByWireNubs "False" 
                .ExtraPreconditioning "False" 
            End With

            With IESolver
                .SetFMMFFCalcStopLevel "0" 
                .SetFMMFFCalcNumInterpPoints "6" 
                .UseFMMFarfieldCalc "True" 
                .SetCFIEAlpha "0.500000" 
                .LowFrequencyStabilization "False" 
                .LowFrequencyStabilizationML "True" 
                .Multilayer "False" 
                .SetiMoMACC_I "0.0001" 
                .SetiMoMACC_M "0.0001" 
                .DeembedExternalPorts "True" 
                .SetOpenBC_XY "True" 
                .OldRCSSweepDefintion "False" 
                .SetRCSOptimizationProperties "True", "100", "0.00001" 
                .SetAccuracySetting "Medium" 
                .CalculateSParaforFieldsources "True" 
                .ModeTrackingCMA "True" 
                .NumberOfModesCMA "3" 
                .StartFrequencyCMA "-1.0" 
                .SetAccuracySettingCMA "Default" 
                .FrequencySamplesCMA "0" 
                .SetMemSettingCMA "Auto" 
                .CalculateModalWeightingCoefficientsCMA "True" 
                .DetectThinDielectrics "True" 
                .UseLegacyRadiatedPowerCalc "False" 
            End With

        """

        self.prj.model3d.add_to_history("set frequency solver", vba_code) 

                        
    def StartTimeSolver(self):
        """Start Time Simulation

        Returns:
            str: String with VBA Code 
        """
        vba_code = f"""
                    With Solver
                    .Start
                    End With
                    """
        self.prj.model3d.add_to_history("start time solver", vba_code)



    # Set Optical Simulation Domain 
    def setOpticalSimulationProperties(self, Parameters):
        # Units Properties
        Length = Parameters['Dimensions'] 
        Frequency = Parameters['Frequency'] 
        Time = Parameters['Time'] 
        Temperature = Parameters['Temperature'] 

        # Set Background
        Type = Parameters["Type Background"]
        Xmin = Parameters["Xmin Background"] 
        Xmax = Parameters["Xmax Background"]
        Ymin = Parameters["Ymin Background"] 
        Ymax = Parameters["Ymax Background"]
        Zmin = Parameters["Zmin Background"] 
        Zmax = Parameters["Zmax Background"]

        # Operational Wavelength
        WavelengthMin = Parameters["Min Wavelength"]
        WavelengtMax = Parameters["Max Wavelength"]

        # Set Boundary
        XminBound = Parameters["Xmin Boundary"] 
        XmaxBound = Parameters["Xmax Boundary"] 
        YminBound = Parameters["Ymin Boundary"] 
        YmaxBound = Parameters["Ymax Boundary"] 
        ZminBound = Parameters["Zmin Boundary"] 
        ZmaxBound = Parameters["Zmax Boundary"] 
        XSimetryBound = Parameters["Xsymmetry Boundary"] 
        YSimetryBound = Parameters["Ysymmetry Boundary"] 
        ZSimetryBound = Parameters["Zsymmetry Boundary"] 

        # Mesh Settings
        MeshType = Parameters["Mesh Type"] 
        CellNear = Parameters["Mesh Cells Near Object"] 
        CellFar = Parameters["Mesh Cells far Object"] 


        # VBA Code
        vba_code = f"""
                    With Units
                        .SetUnit "Length", "{Length}"
                        .SetUnit "Frequency", "{Frequency}"
                        .SetUnit "Voltage", "V"
                        .SetUnit "Resistance", "Ohm"
                        .SetUnit "Inductance", "nH"
                        .SetUnit "Temperature", "{Temperature}"
                        .SetUnit "Time", "{Time}"
                        .SetUnit "Current", "A"
                        .SetUnit "Conductance", "S"
                        .SetUnit "Capacitance", "pF"
                    End With
                ThermalSolver.AmbientTemperature "0"
                Plot.DrawBox "True"
                    With Background
                        .Type "{Type}"
                        .Epsilon "1.0"
                        .Mu "1.0"
                        .Rho "1.204"
                        .ThermalType "Normal"
                        .ThermalConductivity "0.026"
                        .SpecificHeat "1005", "J/K/kg"
                        .XminSpace "{Xmin}"
                        .XmaxSpace "{Xmax}"
                        .YminSpace "{Ymin}"
                        .YmaxSpace "{Ymax}"
                        .ZminSpace "{Zmin}"
                        .ZmaxSpace "{Zmax}"
                    End With
                    With Boundary
                        .Xmin "{XminBound}"
                        .Xmax "{XmaxBound}"
                        .Ymin "{YminBound}"
                        .Ymax "{YmaxBound}"
                        .Zmin "{ZminBound}"
                        .Zmax "{ZmaxBound}"
                        .Xsymmetry "{XSimetryBound}"
                        .Ysymmetry "{YSimetryBound}"
                        .Zsymmetry "{ZSimetryBound}"
                    End With
                    With Solver
                        .Reset
                        .WavelengthRange "{WavelengthMin}", "{WavelengtMax}"
                    End With
                    With Mesh
                        .MergeThinPECLayerFixpoints "True"
                        .RatioLimit "20"
                        .AutomeshRefineAtPecLines "True", "6"
                        .FPBAAvoidNonRegUnite "True"
                        .ConsiderSpaceForLowerMeshLimit "False"
                        .MinimumStepNumber "5"
                        .AnisotropicCurvatureRefinement "True"
                        .AnisotropicCurvatureRefinementFSM "True"
                    End With
                    With MeshSettings
                        .SetMeshType "Hex"
                        .Set "RatioLimitGeometry", "20"
                        .Set "EdgeRefinementOn", "1"
                        .Set "EdgeRefinementRatio", "6"
                    End With
                    With MeshSettings
                        .SetMeshType "Tet"
                        .Set "VolMeshGradation", "1.5"
                        .Set "SrfMeshGradation", "1.5"
                    End With
                    With MeshSettings
                        .SetMeshType "HexTLM"
                        .Set "RatioLimitGeometry", "20"
                    End With
                MeshAdaption3D.SetAdaptionStrategy "Energy"
                ChangeProblemType "Optical"
                    With MeshSettings
                        .SetMeshType "Hex"
                        .Set "Version", 1%
                        .Set "StepsPerWaveNear", "{CellNear}"
                        .Set "StepsPerWaveFar", "{CellFar}"
                        .Set "WavelengthRefinementSameAsNear", "0" 
                        .Set "StepsPerBoxNear", "{CellNear}"
                        .Set "StepsPerBoxFar", "{CellFar}"
                    End With
                    With Mesh
                        .MeshType "{MeshType}"
                    End With
                ChangeSolverType("HF Time Domain")
                """
        self.prj.model3d.add_to_history("set optical simulation", vba_code)





      
    # Set Optical Simulation Domain 
    def setElectricalSimulationProperties(self, Parameters):
        # Units Properties
        Length = Parameters['Dimensions'] 
        Frequency = Parameters['Frequency'] 
        Time = Parameters['Time'] 
        Temperature = Parameters['Temperature'] 

        # Set Background
        Type = Parameters["Type Background"]
        Xmin = Parameters["Xmin Background"] 
        Xmax = Parameters["Xmax Background"]
        Ymin = Parameters["Ymin Background"] 
        Ymax = Parameters["Ymax Background"]
        Zmin = Parameters["Zmin Background"] 
        Zmax = Parameters["Zmax Background"]


        # Operational Frequency 
        FreqMin = Parameters["Min Frequency"]
        FreqMax = Parameters["Max Frequency"]


        # Set Boundary
        XminBound = Parameters["Xmin Boundary"] 
        XmaxBound = Parameters["Xmax Boundary"] 
        YminBound = Parameters["Ymin Boundary"] 
        YmaxBound = Parameters["Ymax Boundary"] 
        ZminBound = Parameters["Zmin Boundary"] 
        ZmaxBound = Parameters["Zmax Boundary"] 
        XSimetryBound = Parameters["Xsymmetry Boundary"] 
        YSimetryBound = Parameters["Ysymmetry Boundary"] 
        ZSimetryBound = Parameters["Zsymmetry Boundary"] 

        # Mesh Settings
        MeshType = Parameters["Mesh Type"] 
        CellNear = Parameters["Mesh Cells Near Object"] 
        CellFar = Parameters["Mesh Cells far Object"] 




        # VBA Code
        vba_code = f"""
                    With Units
                        .SetUnit "Length", "{Length}"
                        .SetUnit "Frequency", "{Frequency}"
                        .SetUnit "Voltage", "V"
                        .SetUnit "Resistance", "Ohm"
                        .SetUnit "Inductance", "nH"
                        .SetUnit "Temperature", "{Temperature}"
                        .SetUnit "Time", "{Time}"
                        .SetUnit "Current", "A"
                        .SetUnit "Conductance", "S"
                        .SetUnit "Capacitance", "pF"
                    End With
                ThermalSolver.AmbientTemperature "0"
                Plot.DrawBox "True"
                    With Background
                        .Type "{Type}"
                        .Epsilon "1.0"
                        .Mu "1.0"
                        .Rho "1.204"
                        .ThermalType "Normal"
                        .ThermalConductivity "0.026"
                        .SpecificHeat "1005", "J/K/kg"
                        .XminSpace "{Xmin}"
                        .XmaxSpace "{Xmax}"
                        .YminSpace "{Ymin}"
                        .YmaxSpace "{Ymax}"
                        .ZminSpace "{Zmin}"
                        .ZmaxSpace "{Zmax}"
                    End With
                    With Boundary
                        .Xmin "{XminBound}"
                        .Xmax "{XmaxBound}"
                        .Ymin "{minBound}"
                        .Ymax "{YmaxBound}"
                        .Zmin "{ZminBound}"
                        .Zmax "{ZmaxBound}"
                        .Xsymmetry "{XSimetryBound}"
                        .Ysymmetry "{YSimetryBound}"
                        .Zsymmetry "{ZSimetryBound}"
                    End With
                    With Solver
                        .Reset
                        .FrequencyRange "{FreqMin}", "{FreqMax}"
                    End With
                        Mesh.MinimumCurvatureRefinement "150"
                        With MeshSettings
                            .SetMeshType "HexTLM"
                            .Set "StepsPerWaveNear", "{CellNear}"
                        .Set "StepsPerWaveFar", "{CellFar}"
                            .Set "StepsPerBoxNear", "{CellNear}"
                        .Set "StepsPerBoxFar", "{CellFar}"
                            .Set "RatioLimitGeometry", "15"
                        End With
                        With MeshSettings
                            .SetMeshType "Hex"
                            .Set "Version", 1%
                        End With
                        With Mesh
                            .MeshType "PBA"
                        End With
                        ChangeSolverType("HF Time Domain")
                """
        self.prj.model3d.add_to_history("set electrical simulation", vba_code)


    def ChangeSolverType(self, Type):
        time = ["TIME", "time", "Time", "t", "T"]
        if Type in time:
            vba_code = f"""
                        ChangeSolverType "HF Time Domain"
                        """
            self.prj.model3d.add_to_history("change solver", vba_code)
        
        else: 
            pass
    



    def setDomainSolverType(self, Domain):

        domain_list = ["Time", "Freq", "EigenMode", "Integral", "Asymptotic", "Multilayer"]
        if Domain in domain_list:
            if Domain == "Time":
                domain_type = "HF Time Domain"
            elif Domain == "Freq":
                domain_type = "HF Frequency Domain"
            elif Domain == "EigenMode":
                domain_type = "HF Eigenmode"
            elif Domain == "Integral":
                domain_type =  "HF IntegralEq"
            elif Domain == "Asymptotic":
                domain_type = "HF Asymptotic"
            elif Domain == "Multilayer":
                domain_type = "HF Multilayer"

        vba_code = f"""
                    ChangeSolverType "{domain_type}"
                    """
        self.prj.model3d.add_to_history("set freq solver", vba_code)


############################################################################
# Monitora
############################################################################


    def CreateEfieldMonitor(self, Parameters):
        """Set Monitor

        Args:
            Parameters (dict): Dictionary with Parameters
                                Parameters["Monitor Wavelength"] : Int/float. This will set the wavelength. For this function
                                and only for this function you need to give the exact number. So 1.55 um will be Parameters["Wavelength"]  = 1.55e-6
                                Parameters["Monitor Frequency"] : Set Frequency of the Monitor
                                Parameters["Domain"] : Str can be "Frequency" or "Wavelength"
                                Parameters["Monitor Type"] : str with monitor Type. Can be one of :
                                    "Efield", "Hfield", "Surfacecurrent", "Powerflow", "Current",
                                    "Powerloss", "Eenergy", "Elossdens", "Lossdens", "Henergy", 
                                    "Farfield", "Fieldsource", "Spacecharge", "ParticleCurrentDensity", "Electrondensity" 
                                
        Raises:
            ValueError: Error massage

        Returns:
            str: String with VBA Code 
        """
        MonitorType = Parameters["Monitor Type"]
        Domain = Parameters["Domain"]
        Types = ["Efield", "Hfield", "Surfacecurrent", "Powerflow", "Current", "Powerloss", "Eenergy", "Elossdens", "Lossdens", "Henergy", "Farfield", "Fieldsource", "Spacecharge", "ParticleCurrentDensity", "Electrondensity" ]
    
        

        if MonitorType in Types:
            if "Monitor Frequency" in Parameters.keys():
            # if Parameters["Monitor Frequency"] != None:
                Freq = Parameters["Monitor Frequency"]
                Name = MonitorType + ' (wl=' + str(Freq) + ')'
                vba_code = f"""
                    With Monitor
                    .Reset
                    .Name "{Name}"
                    .Domain "{Domain}"
                    .FieldType "{MonitorType}"
                    .MonitorValue "{Freq}"
                    .UseSubvolume "False"
                    .Coordinates "Structure"
                    .SetSubvolume "-5.5", "2.5", "-3.5", "3.5", "-2.5", "3.8"
                    .SetSubvolumeOffset "0.0", "0.0", "0.0", "0.0", "0.0", "0.0"
                    .SetSubvolumeInflateWithOffset "False"
                    .Create
                    End With
                    """
                self.prj.model3d.add_to_history("create Efield monitor", vba_code)

            else:
                Wavelength = Parameters["Monitor Wavelength"]
                Name = MonitorType + ' (wl=' + str(Wavelength) + ')'
                vba_code = f"""
                    With Monitor
                    .Reset
                    .Name "{Name}"
                    .Domain "{Domain}"
                    .FieldType "{MonitorType}"
                    .MonitorValue "{Wavelength}"
                    .UseSubvolume "False"
                    .Coordinates "Structure"
                    .SetSubvolume "-5.5", "2.5", "-3.5", "3.5", "-2.5", "3.8"
                    .SetSubvolumeOffset "0.0", "0.0", "0.0", "0.0", "0.0", "0.0"
                    .SetSubvolumeInflateWithOffset "False"
                    .Create
                    End With
                    """
        else:   
            raise ValueError("Parameters['Monitor Type'] is not in allowed types. Please choose one of the following types: ['Efield', 'Hfield', 'Surfacecurrent', 'Powerflow', 'Current', 'Powerloss', 'Eenergy', 'Elossdens', 'Lossdens', 'Henergy', 'Farfield', 'Fieldsource', 'Spacecharge', 'ParticleCurrentDensity', 'Electrondensity']")



############################################################################
# Mesh
############################################################################

    def SetMesh(self, Parameters):
        """Sets the type of the mesh. The user can define the mesh cells per wavelength near and far
            from the simulations object.
        Args:
            Parameters (dict): Dictionarty with Parameters
                                Parameters["Mesh Type"] : str Type of the Mesh. You can choose between:  
                                                        PBA - Hexahedral mesh with Perfect Boundary Approximation
                                                        HexahedralTLM
                                                        CFD
                                                        CFDNew
                                                        Staircase - Hexahedral mesh with staircase cells
                                                        Tetrahedral - Tetrahedral mesh
                                                        Surface - Surface mesh
                                                        SurfaceMLS - urface multi layer mesh
                                                        Planar - Planar 2D mesh
                                Parameters["Mesh Cells Near Object"] : int Cells per Wavelength near the simulation object
                                Parameters["Mesh Cells far Object"] : int Cells per Wavelength far from the simulation object
        Returns:
            str: String with VBA Code 
        """
        

        MeshType = Parameters["Mesh Type"]
        CellNear = Parameters["Mesh Cells Near Object"]
        CellFar = Parameters["Mesh Cells far Object"]


        vba_code = f"""
            With Mesh
            .MeshType "{MeshType}"
            .SetCreator "High Frequency"
            End With
            nWith MeshSettings
            .SetMeshType "Hex"
            .Set "Version", 1%
            .Set "StepsPerWaveNear", "{CellNear}"
            .Set "StepsPerWaveFar", "{CellFar}"
            .Set "WavelengthRefinementSameAsNear", "0" 
            .Set "StepsPerBoxNear", "{CellNear}"
            .Set "StepsPerBoxFar", "{CellFar}"
            .Set "MaxStepNear", "0"
            .Set "MaxStepFar", "0"
            .Set "ModelBoxDescrNear", "maxedge"
            .Set "ModelBoxDescrFar", "maxedge"
            .Set "UseMaxStepAbsolute", "0"
            .Set "GeometryRefinementSameAsNear", "0"
            End With
            """
        self.prj.model3d.add_to_history("set mesh", vba_code)





############################################################################
# Translations of objects
############################################################################

    

    def RotationTranslation(self, Parameters):


        Name = Parameters["Name Object"]
        AngleX = Parameters["Angle X"]
        AngleY = Parameters["Angle Y"]
        AngleZ = Parameters["Angle Z"]

        vba_code = f"""
                    With Transform
                    .Reset
                    .Name "{Name}"
                    .Origin "ShapeCenter"
                    .Center "0", "0", "0"
                    .Angle "{AngleX}", "{AngleY}", "{AngleZ}"
                    .MultipleObjects "False"
                    .GroupObjects "False"
                    .Repetitions "1"
                    .MultipleSelection "False"
                    .Transform "Curve", "Rotate"
                    End With
                    """
        self.prj.model3d.add_to_history("rotational translation", vba_code)




    def Translation(self, Parameters):
        Type = Parameters["Translate Type"]
        Type_list = ["Translate", "Scale", "Rotate"]
        
        
        

        if Type in Type_list:
            if Type == "Translate":
                Name = Parameters["Name Object"]
                PosX = Parameters["Position X"]
                PosY = Parameters["Position Y"]
                PosZ = Parameters["Position Z"]
                vba_code = f"""
                            With Transform
                            .Reset
                            .Name "{Name}"
                            .Vector "{PosX}", "{PosY}", "{PosZ}"
                            .UsePickedPoints "False"
                            .InvertPickedPoints "False"
                            .MultipleObjects "False"
                            .GroupObjects "False"
                            .Repetitions "1"
                            .MultipleSelection "False"
                            .Transform "Shape", "Translate
                            End With
                            """
                self.prj.model3d.add_to_history("translation", vba_code)
                
            elif Type == "Scale":  
                Name = Parameters["Name Object"]
                CenterPos = Parameters["Center Positions"]
                ScaleFactor = Parameters["Scale Factors"]
        
                vba_code = f"""                                        
                            With Transform 
                                 .Reset 
                                 .Name "{Name}" 
                                 .Origin "Free" 
                                 .Center "{CenterPos["X"]}", "{CenterPos["Y"]}", "{CenterPos["Z"]}" 
                                 .ScaleFactor "{ScaleFactor["X"]}", "{ScaleFactor["Y"]}", "{ScaleFactor["Z"]}" 
                                 .MultipleObjects "False" 
                                 .GroupObjects "False" 
                                 .Repetitions "1" 
                                 .MultipleSelection "False" 
                                 .AutoDestination "True" 
                                 .Transform "Shape", "Scale" 
                            End With
                         """
                self.prj.model3d.add_to_history("translation scale {Name}", vba_code)
            elif Type == "Rotate" :
                Name = Parameters["Name Object"]
                CenterPos = Parameters["Center Positions"]
                Angle = Parameters["Angle Values"]

                vba_code = f"""                                        
                            With Transform 
                                 .Reset 
                                 .Name "{Name}" 
                                 .Origin "Free" 
                                 .Center "{CenterPos["X"]}", "{CenterPos["Y"]}", "{CenterPos["Z"]}" 
                                 .Angle "{Angle["X"]}", "{Angle["Y"]}", "{Angle["Z"]}" 
                                 .MultipleObjects "False" 
                                 .GroupObjects "False" 
                                 .Repetitions "1" 
                                 .MultipleSelection "False" 
                                 .AutoDestination "True" 
                                 .Transform "Shape", "Rotate" 
                            End With

                         """
                self.prj.model3d.add_to_history("translation scale {Name}", vba_code)
                        
        
    
    def Translation_Scale(self, Parameters):

        Name = Parameters["Name Object"]
        CenterPos = Parameters["Center Positions"]
        ScaleFactor = Parameters["Scale Factors"]

        vba_code = f"""                                        
                    With Transform 
                         .Reset 
                         .Name "{Name}" 
                         .Origin "Free" 
                         .Center "{CenterPos["X"]}", "{CenterPos["Y"]}", "{CenterPos["Z"]}" 
                         .ScaleFactor "{ScaleFactor["X"]}", "{ScaleFactor["Y"]}", "{ScaleFactor["Z"]}" 
                         .MultipleObjects "False" 
                         .GroupObjects "False" 
                         .Repetitions "1" 
                         .MultipleSelection "False" 
                         .AutoDestination "True" 
                         .Transform "Shape", "Scale" 
                    End With
                 """
        self.prj.model3d.add_to_history("translation scale {Name}", vba_code)
        
        
        
        
    def Translation_Rotation(self, Parameters):

        Name = Parameters["Name Object"]
        CenterPos = Parameters["Center Positions"]
        Angle = Parameters["Angle Values"]

        vba_code = f"""                                        
                    With Transform 
                         .Reset 
                         .Name "{Name}" 
                         .Origin "Free" 
                         .Center "{CenterPos["X"]}", "{CenterPos["Y"]}", "{CenterPos["Z"]}" 
                         .Angle "{Angle["X"]}", "{Angle["Y"]}", "{Angle["Z"]}" 
                         .MultipleObjects "False" 
                         .GroupObjects "False" 
                         .Repetitions "1" 
                         .MultipleSelection "False" 
                         .AutoDestination "True" 
                         .Transform "Shape", "Rotate" 
                    End With

                 """
        self.prj.model3d.add_to_history("translation scale {Name}", vba_code)
    



    def Cut_structures(self, Parameters):
        CutElement = Parameters["Name Cut Structure"]
        CutterElement = Parameters["Name Structure to Cut"]
        vba_code = f"""
                    Solid.Subtract "{CutElement}", "{CutterElement}"
        """
        self.prj.model3d.add_to_history("Cur Structure {CutElement}", vba_code)
        
        




#################################################################################
# Curves Class
#################################################################################

class Curves:
    
    def __init__(self, Span_X, Span_Y, Length):
        self.Span_X = Span_X
        self.Span_Y = Span_Y
        self.Length = Length
        self.t = np.arange(0, 1, 1/self.Length)
        self.x_Points = []
        self.y_Points = []
        
        
    def Bezier_Curve(self):
        """
        The Cubic Bezier Polynoms are taken from Wikipedia!

        Returns
        -------
        curve : 2 dimentional np array
            curve[:,0] - X Param of the Bezier Curve
            curve[:,1] - Y Param of the Bezier Curve

        """
        for i in range(len(self.t)):
            self.x_Points.append(float((1-self.t[i])**3*0 + 3*(1-self.t[i])**2*self.t[i]*self.Span_X/2 + 3*(1-self.t[i])*self.t[i]**2*self.Span_X/2 + self.t[i]**3*self.Span_X))
            self.y_Points.append(float(((1-self.t[i]**3)*0 + 3*(1-self.t[i])**2*self.t[i]*0 + 3*(1-self.t[i])*self.t[i]**2*self.Span_Y + self.t[i]**3*self.Span_Y)))
            
        x_Points = np.array(self.x_Points,dtype=float)
        y_Points = np.array(self.y_Points,dtype=float)  
        curve = np.vstack((x_Points, y_Points)).T
        
        return curve
    
    
    
    def Cosinus_Curve(self):
        '''
        The cosinus function is SpanY*(cos((pi/(2*SpanX))*t)^2)  ---->> t in the range of 0 to SpanY

        Returns
        -------
        curve : 2 dimentional np array
            curve[:,0] - X Param of the Bezier Curve
            curve[:,1] - Y Param of the Bezier Curve

        '''
        P = self.Span_Y
        L = self.Span_X
        stepSize = L/len(self.t)
        t = np.arange(0,L, stepSize) 
        Func = P*(np.cos((np.pi/(2*L))*t)**2)  
        curve = np.vstack((t, Func[::-1])).T
        
        return curve
    
    
    
    
    def Euler_Curve(self):
        """ 
        Returns a smooth, continuous-curvature Euler S-bend from (0, 0) to (Span_X, Span_Y)
        using Fresnel integrals and linear curvature.
        """
        Span_X = self.Span_X
        Span_Y = self.Span_Y
        num_pts = self.Length
        
        # Arc length parameter
        # s = np.arange(-1, 1, 2/num_pts)
        s = np.linspace(-1, 1, num_pts)
        
        # Fresnel integrals (Euler spiral)
        S, C = fresnel(s)

        # Normalize to range [-1, 1]
        C = C - C[0]
        S = S - S[0]
        C = C / (C[-1] - C[0])
        S = S / (S[-1] - S[0])

        # Euler spiral from -90° to +90° → makes an S-bend
        x = S * Span_X
        y = C * Span_Y
        

        curve = np.vstack((x, y)).T
    
        return curve

