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
        """Connect to existing CST Window or open new working window.
           Possible CST Projects:
                                - CST Cable Studio project
                                - CST Design Studio project
                                - CST EM Studio project
                                - Filter Designer 3D project
                                - CST Mphysics Studio project
                                - CST Microwave Studio project
                                - CST PCB Studio project
                                - CST Particle Studio project

        """
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
        """
        Save CST Project to path with an given name

        Parameters
        ----------
        path : str
            Path to where the CST file should be saved.
        name : str
            Name of the file to save
        overwrite : boolen
            If the save file exist you can set "overwrite=True" to overwrite the existing save project. The default is False.

        Returns
        -------
        None.

        """


        self.prj.save(path + "/" + str(name) + ".cst", allow_overwrite=overwrite)



    def Open_Project(self, path):
        """
        Open existing project.

        Parameters
        ----------
        path : str
            Path and name to the project that you want ot open. 
            For example : path = "C:/Test_folder/Test_cst_Project.cst"

        Returns
        -------
        None.

        """
        
        # open the project
        self.prj = self.de.open_project(path)


    def removeObject(self, Parameters):
        """
        Delete objects in project

        Parameters
        ----------
        Parameters : dict 
            Dictionary with all needed parameters for this function.
        Parameters["Type"] : str
            Type of the Deleted objet. It can be:  
                                                                'Folder'
                                                                'Material'
                                                                'Component'
                                                                'Port'
                                                                'Curve' 
        Parameters["Name"] : str
            Name of the object to delete. When 'Port' choosen you only need 
            to give the number of the port , like Parameters["Name"] = "1"
            
            For example Parameters["Type"] = 'Component'
            Parameters["Name"] = 'Box'
            This will delte an component called box. 
        
                                                                    
        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        

        type = Parameters["Type"]
        name = Parameters["Name"]

        if type == "Folder":
            vba_code = f"""
                         Wire.DeleteFolder "{name}"
                        """
            self.prj.model3d.add_to_history("set units", vba_code)

        elif type == "Material":
            vba_code = f"""
                         Material.Delete "{name}"

                        """
            self.prj.model3d.add_to_history("set units", vba_code)

        elif type == "Component":
            vba_code = f"""
                         Component.Delete "{name}"

                        """
            self.prj.model3d.add_to_history("set units", vba_code)

        elif type == "Port":
            vba_code = f"""
                         Port.Delete "{name}"

                        """
            self.prj.model3d.add_to_history("set units", vba_code)

        elif type == "Curve":
            vba_code = f"""
                         Curve.DeleteCurve "{name}"

                        """
            self.prj.model3d.add_to_history("set units", vba_code)
        
        else:
            raise ValueError("Not supportet Type. You can choose between, 'Component', 'Folder', Material, 'Port', or 'Curve' !")
  




    ############################################################################
    # CST solver parameters
    ############################################################################
    

    def set_Units(self, Parameters):
        """
        Set CST enviroment global units.

        Parameters
        ----------
        Parameters : dict 
            Dictionary with all needed parameters for this function.
        Parameters["Unit Lenght"] : str
            Measurement unit for lenght. For example "um"
        Parameters["Unit Frequency"] : str
            Measurement unit for frequency. For example "GHz"
        Parameters["Unit Time"] : str
            Measurement unit for time. For example "ns"
        Parameters["Unit Temperature"] : str
            Measurement unit for temperatur. For example "K"

        Returns
        -------
        None.

        """
        

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
        """
        Add a JSON material file if one is provided by the manufacturer or another source.

        Parameters
        ----------
        path : str
            Path to file.
        name : str
            Name of the file..

        Returns
        -------
        None.

        """
        
        # load materials from JSON file
        with open(os.path.join(path, str(name) + ".json")) as json_file:
            materials = json.load(json_file)



    def add_material(self, name):
        """
        Add an pre-set of Materials that are for now only available. An extra materials need to be added here. 
           

        Parameters
        ----------
        name : str
            Materials that are added to this CST_Constructor:
                - "Si" Silicon
                - "LiNbO3" Lithium niobate x-cut
                - "SiO2" silicon dioxide
                - "Au" Gold
                - "Al" Aluminium
                - "Glue" Special Glue for the Bondwire floating shield according to KIT Model

        Raises
        ------
        ValueError
            Give Error back if you choose some matherial that is not in the list of materials.

        Returns
        -------
        None.

        """
        
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
        """
        Add Anisotropic material to Material Libs

        Parameters
        ----------
        name : str
            Name of the Material
        Values : dict
            Dictionary with the material Values:
                Values["X"] : X Epsilon Value
                Values["Y"] : Y Epsilon Value
                Values["Z"] : Z Epsilon Value

        Returns
        -------
        None.

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
        """
        Add Silicon to Material 

        Returns
        -------
        None.

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
                        .Colour "0.749", "0.749", "0.560" 
                        .Create
                        End With
                    """
        self.prj.model3d.add_to_history("add Si material", vba_code)
    

    def add_SiO2(self):
        """
        Add Silicon dioxide to Materials

        Returns
        -------
        None.

        """
        
        # Add Meterial SiO2
        vba_code = f"""
                        With Material
                        .Reset 
                        .Name "SiO2"
                        .Rho "2270.0"
                        .ThermalType "Normal"
                        .ThermalConductivity "0.32"
                        .SpecificHeat "1000", "J/K/kg"
                        .DynamicViscosity "0"
                        .UseEmissivity "True"
                        .Emissivity "0"
                        .MetabolicRate "0.0"
                        .VoxelConvection "0.0"
                        .BloodFlow "0"
                        .MechanicsType "Isotropic"
                        .YoungsModulus "66"
                        .PoissonsRatio "0.17"
                        .ThermalExpansionRate "0.56"
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
                        .Type "Normal"
                        .MaterialUnit "Frequency", "GHz"
                        .MaterialUnit "Geometry", "um"
                        .MaterialUnit "Time", "s"
                        .Epsilon "3.8"
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
                        .Colour "0.529", "0.808", "0.922" 
                        .Wireframe "False" 
                        .Reflection "False" 
                        .Allowoutline "True" 
                        .Transparentoutline "False" 
                        .Transparency "0" 
                        .Create
                        End With
                    """
        self.prj.model3d.add_to_history("add SiO2 material", vba_code)




    def add_Au(self):
        """
        Add gold to Material

        Returns
        -------
        None.

        """
        
        # Add Meterial Gold.
        vba_code = f"""
                        With Material
                        .Reset
                        .Name "Au"
                        .FrqType "static"
                        .Type "Normal"
                        .SetMaterialUnit "Hz", "mm"
                        .Epsilon "1"
                        .Mu "1.0"
                        .Kappa "4.561e+007"
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
                        .Sigma "4.561e+007"
                        .Rho "19320.0"
                        .ThermalType "Normal"
                        .ThermalConductivity "314.0"
                        .SpecificHeat "130", "J/K/kg"
                        .MetabolicRate "0"
                        .BloodFlow "0"
                        .VoxelConvection "0"
                        .MechanicsType "Isotropic"
                        .YoungsModulus "78"
                        .PoissonsRatio "0.42"
                        .ThermalExpansionRate "14"
                        .ReferenceCoordSystem "Global"
                        .CoordSystemType "Cartesian"
                        .NLAnisotropy "False"
                        .NLAStackingFactor "1"
                        .NLADirectionX "1"
                        .NLADirectionY "0"
                        .NLADirectionZ "0"
                        .Color "1","0.84", "0"
                        .Wireframe "False"
                        .Reflection "False"
                        .Allowoutline "True"
                        .Transparentoutline "False"
                        .Transparency "0"
                        .Create
                        End With
                    """
        self.prj.model3d.add_to_history("add Au material", vba_code)



    def add_Al(self):
        """
        Add Aluminium to Material

        Returns
        -------
        None.

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
                            .Colour "0.643", "0.666", "0.635" 
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
        """
        Add special glue for the floating shield bond wires according to the KIT model.

        Returns
        -------
        None.

        """
        
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
                    .SolarRadiationAbsorptionType "Opaque"
                    .Absorptance "0"
                    .UseSemiTransparencyCalculator "False"
                    .SolarRadTransmittance "0.0"
                    .SolarRadReflectance "0.0"
                    .SolarRadSpecimenThickness "0.0"
                    .SolarRadRefractiveIndex "1.0"
                    .SolarRadAbsorptionCoefficient "0.0"
                    .IntrinsicCarrierDensity "0"
                    .IntrinsicCarrierDensityModel "none"
                    .FrqType "all"
                    .Type "Normal"
                    .MaterialUnit "Frequency", "GHz"
                    .MaterialUnit "Geometry", "um"
                    .MaterialUnit "Time", "ns"
                    .MaterialUnit "Temperature", "degC"
                    .Epsilon "2.6"
                    .Mu "1"
                    .Sigma "0"
                    .TanD "0.043"
                    .TanDFreq "100"
                    .TanDGiven "True"
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
        """
        Add an materials according to GGB Probes. 

        Returns
        -------
        None.

        """
        
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

        Parameters
        ----------
        Parameters : dict
            Dictionary with the name and value of the parameter. 
                For Example:
                    Parameters = {}
                    Parameters['Length'] = 20
                    This will add Length = 20 to your CSt global enviroments. 

        Returns
        -------
        None.

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
        """
        Delete an Global parameter from the enviroment.

        Parameters
        ----------
        Parameters : str
            Name of parameter to be removed.

        Returns
        -------
        None.

        """
        

        # Parameters Values and Keys

        # Param = list(Parameters.keys())
        Param = Parameters
        # tmp = []
        # for i in range(len(Param)):
        #     tmp.append(f'DeleteParameter "{Param[i]}"')
        tmp = f'DeleteParameter "{Param}"'

        # # join them with newlines
        # line_block = "\n".join(tmp)
        
        # vba_code = f"""
        #             {line_block}
        #             """
        vba_code = f"""
                        {tmp}
                    """

        self.prj.model3d.add_to_history("delete global parameter", vba_code)


        
    def setBackground(self, Parameters):
        
        """
        Set the Simulation background.

        Parameters
        ----------
        Parameters : dict 
            Dictionary with all needed parameters for this function.
            Parameters["Type Background"] : str
                Type of the Background. Defoult is "Normal". It can be set to:
                    "PEC"
                    "Normal"
                    "Anisotropic"
                    "Anisotropic"     
            Parameters["Xmax Background"] : int/float
                Background X-max value
            Parameters["Xmin Background"] : int/float
                Background X-min value
            Parameters["Ymax Background"] : int/float
                Background Y-max value
            Parameters["Ymin Background"] : int/float
                Background Y-min value
            Parameters["Zmax Background"] : int/float
                Background Z-max value
            Parameters["Zmin Background"] : int/float
                Background Z-max value
        
        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        

        #Check Background Type
        Types = ["PEC", "Normal", "Anisotropic", "Anisotropic"]
        
        if Parameters["Type Background"] not in Types:
            raise ValueError("Background Type can be set to: 'PEC', 'Normal', 'Anisotropic' or 'Lossy metal'")
        else:
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
        """
        Set Boundary box Parameters

        Parameters
        ----------
        Parameters["Xmin Boundary"] : str
            Type of the X-min Boundary
        Parameters["Xmax Boundary"] : str
            Type of the X-max Boundary
        Parameters["Ymin Boundary"] : str
            Type of the Y-min Boundary
        Parameters["Ymax Boundary"] : str
            Type of the Y-max Boundary
        Parameters["Zmin Boundary"] : str
            Type of the Z-min Boundary
        Parameters["Zmax Boundary"] : str 
            Type of the Z-max Boundary
        Parameters["Xsymmetry Boundary"] : str
            Type of the X Boundary Symmetry 
        Parameters["Ysymmetry Boundary"] : str
            Type of the Y Boundary Symmetry 
        Parameters["Zsymmetry Boundary"] : str
            Type of the Z Boundary Symmetry 

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        
        Types = ["electric", "magnetic", "open", "expanded open", "periodic", "conducting wall", "unit cell"]
        Symmetry = ["none", "electric", "magnetic"]
        list_values = list(Parameters.values())
        
        for v in list_values:
            if v not in Types + Symmetry:  # combine both lists
                raise ValueError(f"""Invalid value '{v}'! 
        You can choose from the Boundary types: {', '.join(Types)}.
        For Symmetry you can choose between: {', '.join(Symmetry)}.""")
    
               
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
        """
        Set Simulation Frequency

        Parameters
        ----------
        Parameters["Min Frequency"] : int/float
            Min Frequency of Simulation
        Parameters["Max Frequency"] : int/float
            Max Frequency of Simulation

        Returns
        -------
        None.

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
        """
        Set Simulation Wavelength

        Parameters
        ----------
        Parameters["Min Wavelength"] : int/float
            Min Wavelength
        Parameters["Max Wavelength"] : int/float
            Max Wavelength
            

        Returns
        -------
        None.

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
        """
        Create an Brick object.

        Parameters
        ----------
        Parameters["Brick Lenght Min"] : int/float
            Min Lenght. It will be center so the lenght that you give will be set to Lenght/2
        Parameters["Brick Lenght Max"] : int/float
            Max Lenght. It will be center so the lenght that you give will be set to Lenght/2
        Parameters["Brick Width Min"] : int/float
            Min Width. It will be center so the width that you give will be set to Width/2
        Parameters["Brick Width Max"] : int/float
            Max Width. It will be center so the width that you give will be set to Width/2
        Parameters["Brick Hight Min"] : int/float
            Min Hight. It will be center so the hight that you give will be set to Hight/2
        Parameters["Brick Hight Max"] : int/float
            Max Hight. It will be center so the hight that you give will be set to Hight/2
        Parameters["Material"] : str
            Set Material for the Brick object
        Parameters["Brick Name"] : str
            Set Brick Name
        Parameters["Component Name"] : str
            Set Component Name. The strcuture will be made as "Component Name:Brick Name"

        Returns
        -------
        None.

        """
        
        
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
        """
        Create Sphere object.

        Parameters
        ----------
        Parameters["Axis"] : str
            Set an oriantation Axis. Can be "X", "Y" or "Z".
        Parameters["Center Radius"] : int/float
            Radius of the sphere
        Parameters["Top Radius"] : int/float
            Set sphere top radius
        Parameters["Bottom Radius"] : int/float
            Set sphere bottom radius
        Parameters["Center Positions"] : dict
            Dictionary with XCenter, YCenter and ZCenter. 
                For example:
                    Positions = {}
                    Positions["X"] = 2
                    Positions["Y"] = 2
                    Positions["Z"] = 2
                    Parameters["Center Positions"] = Positions
        Parameters["Material"] : str
            Set Materials for the Sphere object
        Parameters["Name"] : str
            Set Sphere Name.
        Parameters["Component Name"] : str
            Set Component Name. The strcuture will be made as "Component Name:Name".
        

        Returns
        -------
        None.

        """
        
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
        
        


        

    def Curve(self, Parameters):
        """
        Create 2D Curve in the CST enviroment.

        Parameters
        ----------
        Parameters["Curve Name"] : str
            Set name of the Curves
        Parameters["Points"] : dict
            Dictionary of Curve Points. 
                For example. 
                    from CST_Constructor import Curves
                     # Define Curves Parameters and Data
                    Lenght = 100
                    Offset = 40
                    points = 100
                    # Generate the Bezier and Cos points
                    ObjCurves = Curves(Lenght, Offset, points)
                    CosinusCurve = ObjCurves.Cosinus_Curve()
                    Parameters["Points"] = CosinusCurve
        
        Returns
        -------
        None.

        """
        
        # Set Parameters
        CurveName = Parameters["Curve Name"]
        Points = Parameters["Points"]

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
        """
        Create Elipse object, in Z-orientation in the moment, in CST enviroment. 

        Parameters
        ----------
        Parameters["Name"] : str
            Set Name of the Elipse object.
        Parameters["Curve Name"] : str
            Set Component Name. The strcuture will be made as "Curve Name:Name".
        Parameters["X Radius"] : int/float
            Set radius in X-Direction of the elipse
        Parameters["Y Radius"] : int/float
            Set radius in Y-Direction of the elipse
        Parameters["X Center"] : int/float
            Set position of X-Center
        Parameters["Y Center"] : int/float
            Set position of Y-Center

        Returns
        -------
        None.

        """
        
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
        



    def Poligon_2D(self, Parameters):
        """
        Create the 2D poligon for tRib waveguide

        Parameters
        ----------
        Parameters["Waveguide Name"] : str
            Set name of the poligon.
        Parameters["Points X"] : list of int/float
            Set poligon dictionary wiht X Points
         Parameters["Points Y"] : list of int/float
             Set poligon dictionary wiht X Points

        Returns
        -------
        None.

        """
        

        WGName  = Parameters["Waveguide Name"]
        Points_X = Parameters["Points X"]
        Points_Y = Parameters["Points Y"]

    
        # Extract the first points as starting Points
        Start_PointX = Points_X[0]
        Start_PointY = Points_Y[0]
        lines = []

        for i in range(1, len(Points_X)):
            lines.append(f'.LineTo "{Points_X[i]}", "{Points_Y[i]}"')

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
        Create the 3D poligon for tRib waveguide

        Parameters
        ----------
        Parameters["Name"] : str
            Set name of the poligon.
        Parameters["Curve Name"]  : str
            Set Component Name. The strcuture will be made as "Curve Name:Name".
        Parameters["Point"] : dict
            Dictionary with 'X' and 'Y' points for the poligon. 

        Returns
        -------
        None.

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
        """
        Create an MZM Modulator. Materials used:
                                                Gold - For the electrodes
                                                LiNbO3 - For Optical Waveguides
                                                SiO2 - For Substrate

        Parameters
        ----------
        Parameters["Lenght_Electrodes"] : int/float
            Set the length of the Electrodes. The Waveguides will be 2 (Units) longer then the electrodes.
        Parameters["Width GND"] : int/float
            Set the width of the GND electrodes.
        Parameters["Width Signal"] : int/float
            Set the width of the Signal Electrode.
        Parameters["Width WG"] : int/float
            Set the top width of the optical waveguide. It is an Rib waveguide.
        Parameters["Gap"] : int/float
            Set the gap between Signal and optical Waveguide.
        Parameters["angle"] : int/float
            Set the angle of the side wall of the optical waveguide.
        Parameters["High Electrodes"] : int/float
            Set the hight of the Electodes.
        Parameters["High WG"] : int/float
            Set the hight of the optical Waveguide
        Parameters["High Slab"] : int/float
            Set the hight of the Slab. When choosen "0" no Slab will be implemented.
        Parameters["High Substrate"] : int/float
            Set the hight of the substrate.

        Returns
        -------
        None.

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

        # self.Poligon_2D(WGName = "Waveguide_Left", Points = PointsLeft)
        Parameters_2D = {}
        Parameters_2D["Waveguide Name"] = "Waveguide_Left"
        Parameters_2D["Points X"] = PointsLeft["X"]
        Parameters_2D["Points Y"] = PointsLeft["Y"]
        Parameters_2D["Points"] = PointsLeft
        self.Poligon_2D(Parameters_2D)
        Trans = {}
        Trans["Translate Type"] = "Rotate"
        Trans["Name Object"] = "Waveguide_Left"
        Trans["Center Positions"] = [0,0,0]
        Trans["Angle Values"] = [0, 90, 0]
        Trans["Object to rotate"] = "Curve"
        Trans['Origin'] = "ShapeCenter"
        self.Translation(Trans)
        Trans = {}
        Trans["Translate Type"] = "Translate"
        Trans["Name Object"] = "Waveguide_Left"
        Trans["Position"] = [Height_WG/2, 0, HeightZ/2 + Height_WG/2]
        Trans["Object to rotate"] = "Curve"
        self.Translation(Trans)
        Parameters_RibWaveguide_toSolid = {}
        Parameters_RibWaveguide_toSolid["Waveguide Name"] = "Waveguide_Left"
        Parameters_RibWaveguide_toSolid["Waveguide Hight"] = Length_WG
        Parameters_RibWaveguide_toSolid["Angle"] = 0
        Parameters_RibWaveguide_toSolid["Name Folder"] = "Waveguide_Left"
        Parameters_RibWaveguide_toSolid["Material"]  ="LiNbO3"
        Parameters_RibWaveguide_toSolid["Waveguide Folder Name"] = "Waveguide_Left"
        Parameters_RibWaveguide_toSolid["Waveguide Name"] = "Waveguide_Left"
        self.RibWaveguide_ToSolid(Parameters_RibWaveguide_toSolid)
        


        PointsRight = {}
        PointsRight["X"] = [Length_WG/2, Length_WG/2, Length_WG/2 - Height_WG, Length_WG/2 - Height_WG, Length_WG/2]
        PointsRight["Y"] = [PosRight[0]/2, PosRight[1]/2, round((PosRight[1] - extention),2)/2 , round((PosRight[0] + extention),2)/2, PosRight[0]/2]
        

        # Waveguide and Waveguide to solid
        Parameters_2D = {}
        Parameters_2D["Waveguide Name"] = "Waveguide_Right"
        Parameters_2D["Points X"] = PointsRight["X"]
        Parameters_2D["Points Y"] = PointsRight["Y"]
        self.Poligon_2D(Parameters_2D)
        Trans = {}
        Trans["Translate Type"] = "Rotate"
        Trans["Name Object"] = "Waveguide_Right"
        Trans["Center Positions"] = [0,0,0]
        Trans["Angle Values"] = [0, 90, 0]
        Trans["Object to rotate"] = "Curve"
        Trans['Origin'] = "ShapeCenter"
        self.Translation(Trans)
        Trans = {}
        Trans["Translate Type"] = "Translate"
        Trans["Name Object"] = "Waveguide_Right"
        Trans["Position"] = [Height_WG/2, 0, HeightZ/2 + Height_WG/2]
        Trans["Object to rotate"] = "Curve"
        self.Translation(Trans)
        Parameters_RibWaveguide_toSolid = {}
        Parameters_RibWaveguide_toSolid["Waveguide Name"] = "Waveguide_Right"
        Parameters_RibWaveguide_toSolid["Waveguide Hight"] = -Length_WG
        Parameters_RibWaveguide_toSolid["Angle"] = 0
        Parameters_RibWaveguide_toSolid["Name Folder"] = "Waveguide_Right"
        Parameters_RibWaveguide_toSolid["Material"]  ="LiNbO3"
        Parameters_RibWaveguide_toSolid["Waveguide Folder Name"] = "Waveguide_Right"
        Parameters_RibWaveguide_toSolid["Waveguide Name"] = "Waveguide_Right"
        self.RibWaveguide_ToSolid(Parameters_RibWaveguide_toSolid)
        # self.RibWaveguide_ToSolid("Waveguide_Right", WaveguideName = "Waveguide_Right", WG_Hight = -Length_WG, Angle = 0, WGFolderName = "Waveguide_Right", WGName = "Waveguide_Right", Material="LiNbO3")





    # Phase Modulator Design
    def PhaseModulator(self, Parameters):
        """
        Create an Phase Modulator. Materials used:
                                                Gold - For the electrodes
                                                LiNbO3 - For Optical Waveguides
                                                SiO2 - For Substrate

        Parameters
        ----------
        Parameters["Lenght_Electrodes"] : int/float
            Set the length of the Electrodes. The Waveguides will be 2 (Units) longer then the electrodes.
        Parameters["Width GND"] : int/float
            Set the width of the GND electrodes.
        Parameters["Width Signal"] : int/float
            Set the width of the Signal Electrode.
        Parameters["Width WG"] : int/float
            Set the top wWidth of the optical waveguide. It is an Rib waveguide.
        Parameters["Gap"] : int/float
            Set the gap between Signal and optical Waveguide.
        Parameters["angle"] : int/float
            Set the angle of the side wall of the optical waveguide.
        Parameters["High Electrodes"] : int/float
            Set the hight of the Electodes.
        Parameters["High WG"] : int/float
            Set the hight of the optical Waveguide.
        Parameters["High Slab"] : int/float
            Set the hight of the Slab. When choosen "0" no Slab will be implemented.
        Parameters["High Substrate"] : int/float
            Set the hight of the substrate

        Returns
        -------
        None.

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

        Parameters_2D = {}
        Parameters_2D["Waveguide Name"] = "Waveguide_Left"
        Parameters_2D["Points X"] = PointsLeft["X"]
        Parameters_2D["Points Y"] = PointsLeft["Y"]
        self.Poligon_2D(Parameters_2D)

        Trans["Translate Type"] = "Rotate"
        Trans["Name Object"] = "Waveguide_Left"
        Trans["Center Positions"] = [0,0,0]
        Trans["Angle Values"] = [0, 90, 0]
        Trans["Object to rotate"] = "Curve"
        Trans['Origin'] = "ShapeCenter"
        self.Translation(Trans)
        Trans = {}
        Trans["Translate Type"] = "Translate"
        Trans["Name Object"] = "Waveguide_Left"
        Trans["Position"] = [Height_WG/2, 0, HeightZ/2 + Height_WG/2]
        Trans["Object to rotate"] = "Curve"
        self.Translation(Trans)
        
        Parameters_RibWaveguide_toSolid = {}
        Parameters_RibWaveguide_toSolid["Waveguide Name"] = "Waveguide_Left"
        Parameters_RibWaveguide_toSolid["Waveguide Hight"] = Length_WG
        Parameters_RibWaveguide_toSolid["Angle"] = 0
        Parameters_RibWaveguide_toSolid["Name Folder"] = "Waveguide_Left"
        Parameters_RibWaveguide_toSolid["Material"]  ="LiNbO3"
        Parameters_RibWaveguide_toSolid["Waveguide Folder Name"] = "Waveguide_Left"
        Parameters_RibWaveguide_toSolid["Waveguide Name"] = "Waveguide_Left"
        self.RibWaveguide_ToSolid(Parameters_RibWaveguide_toSolid)
     


    def Squere_Waveguide(self, Parameters):
        """
        This function generate and simple straight waveguide. Materials used:
                                                Gold - For the electrodes
                                                LiNbO3 - For Optical Waveguides
                                                SiO2 - For Substrate

        Parameters
        ----------
        Parameters["Lenght WG"] : int/float
            Set the length of the Waveguide.
        Parameters["High_WG"] : int/float
            Set the hight of the optical Waveguide.
        Parameters["Width WG"] : int/float
            Set the top width of the optical waveguide. It is an Rib waveguide.
        Parameters["Substrate Height"] : int/float
            Set the hight of the substrate..
        Parameters["Slab Heigh"] : int/float
            Set the hight of the slab. When choosen "0" no Slab will be implemented. .

        Returns
        -------
        None.

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



    def GSG_Bondwire_ChipToChip_connection(self, Parameters):
        """
        This function generate and two Chips connected by bondwires setup
        The Setup is taken from IHP the Bondpads are from Aluminiam the 
        bondwires are also from Aluminium. The Glue and floating shield over the 
        bondwires was taken from KIT Paper setup. This setup is for GSG Pads.
        This Function have the Option to import the GGP GSg probes. 
        Materials:
                    Gold - Floating Shield
                    DAF Glue - For the glue between the bondwires
                    Al - For bondwires and bondpads
                    SiO2 - Layer under the bondpads on the boths chips
                    Si - Substrate Layer below the SiO2

        Parameters
        ----------
        Parameters["PAD Width GND"] : int/float
            Set the GND Bondpads width.
        Parameters["PAD Width Signal"] : int/float
            Set the signal bondpads width.
        Parameters["PAD Length"] : int/float
            Set the signal and GND bondpads lenght.
        Parameters["PAD Thickness"] : int/float
            Set the signal and GND bondpads thickness.
        Parameters["PADs Distance"] : int/float
            Set the distance from first chip bondpads to the secound chip bondpads.
        Parameters["Bonwire height"] : int/float
            Set the bondwire hight in the middle point.
        Parameters["Bonwire radius"] : int/float
            Set the bondwire ridius.
        Parameters["Floating Shield"] : boolen
            Set an floating shield with glue as dielectricum to improve the Bandwidth of the Bondwires connections.
            If set True you need to set Parameters["Floating Shield Thickness"] and Parameters["Glue Thickness"] too!
        Parameters["Glue Thickness"] : int/float
            Set the DAF glue thickness: measured from the SiO₂ layer of the first chip upward to the SiO₂ layer of the second chip.
        Parameters["Floating Shield Thickness"] : int/float
            Set the floating shield gold metal thickness.
        Parameters["Accuracy"] : int/float
            FDTD Solver Accuracy. Can be :
                                    80 dB
                                    60 dB
                                    50 dB
                                    40 dB
                                    35 dB
                                    30 dB
                                    25 dB
                                    20 dB
                                    'no check'

        Returns
        -------
        None.

        """
        

        # Define Curves Parameters and Data
        Lenght = 100
        Offset = 40
        points = 100
        # Generate the Bezier and Cos points
        ObjCurves = Curves(Lenght, Offset, points)
        CosinusCurve = ObjCurves.Cosinus_Curve()
 

        # Import Materials into CST enviroment
        self.add_Al()
        self.add_Au()
        self.add_Glue()
        self.add_Si()
        self.add_SiO2()

        PAD_Width = Parameters["PAD Width GND"] 
        SIG_PAD_Width = Parameters["PAD Width Signal"] 
        PAD_Length = Parameters["PAD Length"]
        PAD_Thickness = Parameters["PAD Thickness"]
        Component_Name = ["_input", "_output"]
        Names_Chips = ["Chip_L", "Chip_R"]
        PAD_Dist_given = Parameters["PADs Distance"]
        PAD_Dist = [ -PAD_Dist_given - PAD_Length, PAD_Dist_given + PAD_Length]
        Bondwire_Height = Parameters["Bonwire height"]
        Bondwire_Radius = Parameters["Bonwire radius"]
        SpanY = Parameters["Port Y Span"]
        SpanZ = Parameters["Port Z Span"]
        Facet_Number = [4, 6]
        Solver_Accuracy = Parameters["Accuracy"] 
        Solver_Impedance_Status = True
        Solver_Impedance_Value = 50
        Solver_Source_Port = "All"
        Solver_Mesh_Type = "TLM"
        Solver_Calc_Modes = False
        Probes = Parameters["Probes"]

    
        # Set Chip Distance from Pad to pad and the spacing
        # this will set the Chip width too
        Dist = np.arange(-200, 400, 200)


        # Create Pads Loop
        Brick_Parameters = {}
        Material_Pads = "Al"
        Name = ["GND_L", "Sig", "GND_R"]

        for j in range(len(Component_Name)):
            for i in range(len(Name)):
                if Name[i].split("_")[0] == "GND":
                    # Create squere Electrodes 
                    Brick_Parameters["Brick Lenght Max"] = PAD_Dist[j] + PAD_Length
                    Brick_Parameters["Brick Lenght Min"] = PAD_Dist[j] - PAD_Length
                    Brick_Parameters["Brick Width Max"] = Dist[i] + PAD_Width
                    Brick_Parameters["Brick Width Min"] = Dist[i] - PAD_Width
                    Brick_Parameters["Brick Hight Max"] = PAD_Thickness * 2
                    Brick_Parameters["Brick Hight Min"] = 0 
                    Brick_Parameters["Brick Name"] = Name[i] + Component_Name[j]
                    Brick_Parameters["Component Name"] = Name[i] 
                    Brick_Parameters["Material"] = Material_Pads
                    self.Brick(Brick_Parameters)
                else:
                    # Create squere Electrodes 
                    Brick_Parameters["Brick Lenght Max"] = PAD_Dist[j] + PAD_Length
                    Brick_Parameters["Brick Lenght Min"] = PAD_Dist[j] - PAD_Length
                    Brick_Parameters["Brick Width Max"] = Dist[i] + SIG_PAD_Width
                    Brick_Parameters["Brick Width Min"] = Dist[i] - SIG_PAD_Width
                    Brick_Parameters["Brick Hight Max"] = PAD_Thickness * 2
                    Brick_Parameters["Brick Hight Min"] = 0 
                    Brick_Parameters["Brick Name"] = Name[i] + Component_Name[j]
                    Brick_Parameters["Component Name"] = Name[i] 
                    Brick_Parameters["Material"] = Material_Pads
                    self.Brick(Brick_Parameters)
        

        # Create two Chips Left and Right Loop 
        Material_Substrate = "Si"
        Material_iso = "SiO2"
        Brick_Parameters = {}

        for i in range(len(Names_Chips)):
            # Create SiO2 Layer
            # Parameters["Brick Lenght Max"] = PAD_Dist[1] + PAD_Length
            # Parameters["Brick Lenght Min"] = PAD_Dist[0] - PAD_Length
            if PAD_Dist[i] > 0:
                Brick_Parameters["Brick Lenght Max"] = PAD_Dist[i] + PAD_Length 
                Brick_Parameters["Brick Lenght Min"] = PAD_Dist[i] - PAD_Length - 50
            else:
                Brick_Parameters["Brick Lenght Max"] = PAD_Dist[i] - PAD_Length 
                Brick_Parameters["Brick Lenght Min"] = PAD_Dist[i] + PAD_Length + 50
            Brick_Parameters["Brick Width Max"] = max(Dist) + PAD_Width*3
            Brick_Parameters["Brick Width Min"] = min(Dist) - PAD_Width*3
            Brick_Parameters["Brick Hight Max"] = 0
            Brick_Parameters["Brick Hight Min"] = -19.92
            Brick_Parameters["Brick Name"] = "SiO2_Layer" + Component_Name[i]
            Brick_Parameters["Component Name"] = Names_Chips[i]
            Brick_Parameters["Material"] = Material_iso
            self.Brick(Brick_Parameters)

            # Create Substrate
            # Parameters["Brick Lenght Max"] = PAD_Dist[1] + PAD_Length
            # Parameters["Brick Lenght Min"] = PAD_Dist[0] - PAD_Length
            if PAD_Dist[i] > 0:
                Brick_Parameters["Brick Lenght Max"] = PAD_Dist[i] + PAD_Length 
                Brick_Parameters["Brick Lenght Min"] = PAD_Dist[i] - PAD_Length - 50
            else:
                Brick_Parameters["Brick Lenght Max"] = PAD_Dist[i] - PAD_Length 
                Brick_Parameters["Brick Lenght Min"] = PAD_Dist[i] + PAD_Length + 50
            Brick_Parameters["Brick Width Max"] = max(Dist) + PAD_Width*3
            Brick_Parameters["Brick Width Min"] = min(Dist) - PAD_Width*3
            Brick_Parameters["Brick Hight Max"] = -19.92
            Brick_Parameters["Brick Hight Min"] = -298.86
            Brick_Parameters["Brick Name"] = "Substrate_Chip" + Component_Name[i]
            Brick_Parameters["Component Name"] = Names_Chips[i]
            Brick_Parameters["Material"] = Material_Substrate
            self.Brick(Brick_Parameters)

        # Create GND PADS Bondwires
        Names_Wires = ["Wire_L_GND_1", "Wire_L_GND_2", "Wire_R_GND_1", "Wire_R_GND_2"]
        Dist2 = [-125, -75, 75, 125]

        for k in range(len(Dist2)):

            Parameters_Bondwire = {}
            Parameters_Bondwire['X1'] = PAD_Dist[0]/2 
            Parameters_Bondwire['Y1'] = Dist2[k]
            Parameters_Bondwire['Z1'] = PAD_Thickness
            Parameters_Bondwire['X2'] = PAD_Dist[1]/2 
            Parameters_Bondwire['Y2'] = Dist2[k]
            Parameters_Bondwire['Z2'] = PAD_Thickness

            ParametersWire = {}
            ParametersWire["Name Wire"] = Names_Wires[k]
            ParametersWire["Coordinates"] = Parameters_Bondwire
            ParametersWire["Height"] = Bondwire_Height
            ParametersWire["Radius"] = Bondwire_Radius
            ParametersWire["Bondwire Type"] = "Spline"
            ParametersWire["Center Position"] = 0.5
            ParametersWire["Material"] = "Al"
            ParametersWire["SolidWireModel"] = True
            ParametersWire["Termination"] = "extended"
            ParametersWire["NameFolder"] = Names_Wires[k] + "_BondWire"
            self.BondWire(ParametersWire)
     

            Parameters_toSolid = {}
            Parameters_toSolid["Solid Name"] = Names_Wires[k]
            Parameters_toSolid["Curve Name"] = Names_Wires[k]
            Parameters_toSolid["Name Folder"] =  Names_Wires[k] + "_BondWire"
            Parameters_toSolid["Material"] =  "Al"
            self.ToSolid(Parameters_toSolid)
   
        
        # Create Signal PADS Bondwires
        Names_Wires = ["Wire_Sig"]
        Dist2 = [0]


        for k in range(len(Dist2)):
            Parameters_Bondwire = {}
            Parameters_Bondwire['X1'] = PAD_Dist[0]/2 
            Parameters_Bondwire['Y1'] = Dist2[k]
            Parameters_Bondwire['Z1'] = PAD_Thickness
            Parameters_Bondwire['X2'] = PAD_Dist[1]/2 
            Parameters_Bondwire['Y2'] = Dist2[k]
            Parameters_Bondwire['Z2'] = PAD_Thickness

            ParametersWire = {}
            ParametersWire["Name Wire"] = Names_Wires[k]
            ParametersWire["Coordinates"] = Parameters_Bondwire
            ParametersWire["Height"] = Bondwire_Height
            ParametersWire["Radius"] = Bondwire_Radius
            ParametersWire["Bondwire Type"] = "Spline"
            ParametersWire["Center Position"] = 0.5
            ParametersWire["Material"] = "Al"
            ParametersWire["SolidWireModel"] = True
            ParametersWire["Termination"] = "extended"
            ParametersWire["NameFolder"] = Names_Wires[k] + "_BondWire"
            self.BondWire(ParametersWire)
     

            Parameters_toSolid = {}
            Parameters_toSolid["Solid Name"] = Names_Wires[k]
            Parameters_toSolid["Curve Name"] = Names_Wires[k]
            Parameters_toSolid["Name Folder"] =  Names_Wires[k] + "_BondWire"
            Parameters_toSolid["Material"] =  "Al"
            self.ToSolid(Parameters_toSolid)

        # Create Top Plate for Bond Wires and cladding
         # Create Top Plate for Bond Wires and cladding
        if Parameters["Floating Shield"] == True:

            # Set glue thickness and floating shield thickness
            Glue_Thickness = Parameters["Glue Thickness"]
            FloatingShieldThickness = Parameters["Floating Shield Thickness"]
            # PAD_Thickness = 2.8
            Name_Shield = "Floating_Shield"
            Component_Name = "Floating_Shield"
            Material = "Au"
            Material_Clad = "DAF_Glue"


            # Create cladding
            Brick_Parameters = {}
            Brick_Parameters["Brick Lenght Max"] = PAD_Dist[0]/2 - PAD_Length
            Brick_Parameters["Brick Lenght Min"] = PAD_Dist[1]/2 + PAD_Length
            Brick_Parameters["Brick Width Max"] = max(Dist) + PAD_Width*2
            Brick_Parameters["Brick Width Min"] = min(Dist) - PAD_Width*2
            Brick_Parameters["Brick Hight Max"] = 2 * (PAD_Thickness + Bondwire_Height + Bondwire_Radius*2) + Glue_Thickness*2
            Brick_Parameters["Brick Hight Min"] = 0
            # Brick_Parameters["Brick Hight Min"] = PAD_Thickness*2 
            Brick_Parameters["Brick Name"] = "Floating_Shield_Clad"
            Brick_Parameters["Component Name"] = Component_Name
            Brick_Parameters["Material"] = Material_Clad
            self.Brick(Brick_Parameters)



            # Create squere floating shield
            Brick_Parameters["Brick Lenght Max"] = PAD_Dist[0]/2 - PAD_Length
            Brick_Parameters["Brick Lenght Min"] = PAD_Dist[1]/2 + PAD_Length
            Brick_Parameters["Brick Width Max"] = max(Dist) + PAD_Width*2
            Brick_Parameters["Brick Width Min"] = min(Dist) - PAD_Width*2
            Brick_Parameters["Brick Hight Max"] = 2 * (PAD_Thickness + Bondwire_Height + Bondwire_Radius*2) + Glue_Thickness*2 + FloatingShieldThickness*2
            Brick_Parameters["Brick Hight Min"] = 2 * (PAD_Thickness + Bondwire_Height + Bondwire_Radius*2) + Glue_Thickness*2
            # Brick_Parameters["Brick Hight Max"] = PAD_Thickness*2 + Glue_Thickness*2 + FloatingShieldThickness*2
            # Brick_Parameters["Brick Hight Min"] = PAD_Thickness*2 + Glue_Thickness*2
            Brick_Parameters["Brick Name"] = Name_Shield
            Brick_Parameters["Component Name"] = Component_Name
            Brick_Parameters["Material"] = Material
            self.Brick(Brick_Parameters)
        else:
            pass

        if Probes == False:
            # Pick Faces for Input Waveguide Port 1
            Parameters_Port = {}
            Facer_input = []
            Parameters_Port["Option"] = "Face"

            for i in range(len(Name)):
                Facer_input.append(Name[i] + ":" + Name[i] + "_input")

            # pick input faces
            for i in range(len(Name)):
                Parameters_Port["Face Number"] = Facet_Number[0]
                Parameters_Port["Object"] = Facer_input[i]
                self.Pick(Parameters_Port)

            # Create Port Input
            Parameters_Port = {}
            Parameters_Port["Port Number"] = 1
            Parameters_Port["Coordinates"] = "Picks"
            Parameters_Port["Orientation"] = "Positive"
            Parameters_Port["Span"] = [3, SpanY, SpanZ]
            Parameters_Port["Picked Port Number"] = [1,1,1]
            Parameters_Port["Number of picks"] = 3
            Parameters_Port["Picked Port Polarity"] = ["negative", "positive", "negative"]
            Parameters_Port["Picked Component Name"] = [ "GND_L:GND_L_input", "Sig:Sig_input", "GND_R:GND_R_input"]
            Parameters_Port["Face Number"] = [Facet_Number[0], Facet_Number[0], Facet_Number[0]]
            self.WaveguidePortWithPins(Parameters_Port)



            # Pick Faces for Input Waveguide Port 2
            Parameters_Port = {}
            Facer_output = []
            Parameters_Port["Option"] = "Face"
            for i in range(len(Name)):
                Facer_output.append(Name[i] + ":" + Name[i] + "_output")

            for i in range(len(Name)):
                Parameters_Port["Face Number"] = Facet_Number[1]
                Parameters_Port["Object"] = Facer_output[i]
                self.Pick(Parameters_Port)

            #Create Port Output
            Parameters_Port = {}
            Parameters_Port["Port Number"] = 2
            Parameters_Port["Coordinates"] = "Picks"
            Parameters_Port["Orientation"] = "Positive"
            Parameters_Port["Span"] = [3, SpanY, SpanZ]
            Parameters_Port["Picked Port Number"] = [1,1,1]
            Parameters_Port["Number of picks"] = 5
            Parameters_Port["Picked Port Polarity"] = ["negative", "positive", "negative"]
            Parameters_Port["Picked Component Name"] = [ "GND_L:GND_L_output", "Sig:Sig_output", "GND_R:GND_R_output"]
            Parameters_Port["Face Number"] = [Facet_Number[1], Facet_Number[1], Facet_Number[1]]
            self.WaveguidePortWithPins(Parameters_Port)


        elif Probes == True:

            # create_GGB_Probe()
            Parameters_Probes = {}
            Parameters_Probes["Component Name"] = "Probe_Left"
            Parameters_Probes["Orientation Angle"] = -30
            Parameters_Probes["Name"] = "GGB_L"
            self.GGB_Probe(Parameters_Probes)

            # Move Probes Left
            Parameters_Probes_Move = {}
            Parameters_Probes_Move["Translate Type"] = "Translate"
            Parameters_Probes_Move["Name Object"] = "Probe_Left"
            # Parameters_Probes_Move["Position X"] = PAD_Dist[0] - PAD_Length/2 + 20
            # Parameters_Probes_Move["Position Y"] = 0
            # Parameters_Probes_Move["Position Z"] = PAD_Thickness/2
            Parameters_Probes_Move["Position"] = [PAD_Dist[0] - PAD_Length/2 + 20, 0, PAD_Thickness/2]
            Parameters_Probes_Move["Object to rotate"] = "Shape"
            self.Translation(Parameters_Probes_Move)
            


            # create_GGB_Probe()
            Parameters_Probes = {}
            Parameters_Probes["Component Name"] = "Probe_Right"
            Parameters_Probes["Orientation Angle"] = 30
            self.GGB_Probe(Parameters_Probes)

            # Move Probes Right
            Parameters_Probes_Move = {}
            Parameters_Probes_Move["Translate Type"] = "Translate"
            Parameters_Probes_Move["Name Object"] = "Probe_Right"
            # Parameters_Probes_Move["Position X"] = PAD_Dist[1] + PAD_Length/2 - 20
            # Parameters_Probes_Move["Position Y"] = 0
            # Parameters_Probes_Move["Position Z"] = PAD_Thickness/2
            Parameters_Probes_Move["Position"] = [PAD_Dist[1] + PAD_Length/2 - 20, 0, PAD_Thickness/2]
            Parameters_Probes_Move["Object to rotate"] = "Shape"
            self.Translation(Parameters_Probes_Move)


        # Set FDTD Time domain solver 
        Parameters_Solver = {}
        Parameters_Solver["Accuracy"] = Solver_Accuracy
        Parameters_Solver["Auto Impedance"] = Solver_Impedance_Status
        Parameters_Solver["Impedance"] = Solver_Impedance_Value
        Parameters_Solver["Source Port"] = Solver_Source_Port
        Parameters_Solver["Solver Mesh Type"] = Solver_Mesh_Type
        Parameters_Solver["Caclculate Modes Only"] = Solver_Calc_Modes
        self.setTimeSolver(Parameters_Solver)





    def GSGSG_Bondwire_ChipToChip_connection(self, Parameters):
        """
        This function generate and two Chips connected by bondwires setup
        The Setup is taken from IHP the Bondpads are from Aluminiam the 
        bondwires are also from Aluminium. The Glue and floating shield over the 
        bondwires was taken from KIT Paper setup. This setup is for GSGSG Pads.
        Materials:
                    Gold - Floating Shield
                    DAF Glue - For the glue between the bondwires
                    Al - For bondwires and bondpads
                    SiO2 - Layer under the bondpads on the boths chips
                    Si - Substrate Layer below the SiO2

        Parameters
        ----------
        Parameters["PAD Width GND"] : int/float
            Set the GND Bondpads width.
        Parameters["PAD Width Signal"] : int/float
            Set the signal bondpads width.
        Parameters["PAD Length"] : int/float
            Set the signal and GND bondpads lenght.
        Parameters["PAD Thickness"] : int/float
            Set the signal and GND bondpads thickness.
        Parameters["PADs Distance"] : int/float
            Set the distance from first chip bondpads to the secound chip bondpads.
        Parameters["Bonwire height"] : int/float
            Set the bondwire hight in the middle point.
        Parameters["Bonwire radius"] : int/float
            Set the bondwire ridius.
        Parameters["Floating Shield"] : boolen
            Set an floating shield with glue as dielectricum to improve the Bandwidth of the Bondwires connections.
            If set True you need to set Parameters["Floating Shield Thickness"] and Parameters["Glue Thickness"] too!
        Parameters["Glue Thickness"] : int/float
            Set the DAF glue thickness: measured from the SiO₂ layer of the first chip upward to the SiO₂ layer of the second chip.
        Parameters["Floating Shield Thickness"] : int/float
            Set the floating shield gold metal thickness.
        Parameters["Accuracy"] : int/float
            Set the FDTD Solver Accuracy. Can be :
                                            80 dB
                                            60 dB
                                            50 dB
                                            40 dB
                                            35 dB
                                            30 dB
                                            25 dB
                                            20 dB
                                            'no check'

        Returns
        -------
        None.

        """
        
        # Define Curves Parameters and Data
        Lenght = 100
        Offset = 40
        points = 100
        # Generate the Bezier and Cos points
        ObjCurves = Curves(Lenght, Offset, points)
        CosinusCurve = ObjCurves.Cosinus_Curve()
 

        # Import Materials into CST enviroment
        self.add_Al()
        self.add_Au()
        self.add_Glue()
        self.add_Si()
        self.add_SiO2()

        PAD_Width = Parameters["PAD Width GND"] 
        SIG_PAD_Width = Parameters["PAD Width Signal"] 
        PAD_Length = Parameters["PAD Length"]
        PAD_Thickness = Parameters["PAD Thickness"]
        Component_Name = ["_input", "_output"]
        Names_Chips = ["Chip_L", "Chip_R"]
        PAD_Dist_given = Parameters["PADs Distance"]
        PAD_Dist = [ -PAD_Dist_given - PAD_Length, PAD_Dist_given + PAD_Length]
        Bondwire_Height = Parameters["Bonwire height"]
        Bondwire_Radius = Parameters["Bonwire radius"]
        
        SpanY = Parameters["Port Y Span"]
        SpanZ = Parameters["Port Z Span"]
        
        Facet_Number = [4, 6]
        Solver_Accuracy = Parameters["Accuracy"] 
        Solver_Impedance_Status = True
        Solver_Impedance_Value = 50
        Solver_Source_Port = "All"
        Solver_Mesh_Type = "TLM"
        Solver_Calc_Modes = False

        # Set Chip Distance from Pad to pad and the spacing
        # this will set the Chip width too
        Dist = np.arange(-400, 500, 200)


        # Create Pads Loop
        Brick_Parameters = {}
        Material_Pads = "Al"
        Name = ["GND_L", "Sig_L", "GND_Mid", "Sig_R", "GND_R"]

        for j in range(len(Component_Name)):
            for i in range(len(Name)):
                if Name[i].split("_")[0] == "GND":
                    # Create squere Electrodes 
                    Brick_Parameters["Brick Lenght Max"] = PAD_Dist[j] + PAD_Length
                    Brick_Parameters["Brick Lenght Min"] = PAD_Dist[j] - PAD_Length
                    Brick_Parameters["Brick Width Max"] = Dist[i] + PAD_Width
                    Brick_Parameters["Brick Width Min"] = Dist[i] - PAD_Width
                    Brick_Parameters["Brick Hight Max"] = PAD_Thickness * 2
                    Brick_Parameters["Brick Hight Min"] = 0 
                    Brick_Parameters["Brick Name"] = Name[i] + Component_Name[j]
                    Brick_Parameters["Component Name"] = Name[i] 
                    Brick_Parameters["Material"] = Material_Pads
                    self.Brick(Brick_Parameters)
                else:
                    # Create squere Electrodes 
                    Brick_Parameters["Brick Lenght Max"] = PAD_Dist[j] + PAD_Length
                    Brick_Parameters["Brick Lenght Min"] = PAD_Dist[j] - PAD_Length
                    Brick_Parameters["Brick Width Max"] = Dist[i] + SIG_PAD_Width
                    Brick_Parameters["Brick Width Min"] = Dist[i] - SIG_PAD_Width
                    Brick_Parameters["Brick Hight Max"] = PAD_Thickness * 2
                    Brick_Parameters["Brick Hight Min"] = 0 
                    Brick_Parameters["Brick Name"] = Name[i] + Component_Name[j]
                    Brick_Parameters["Component Name"] = Name[i] 
                    Brick_Parameters["Material"] = Material_Pads
                    self.Brick(Brick_Parameters)




        # Create two Chips Left and Right Loop 
        Material_Substrate = "Si"
        Material_iso = "SiO2"
        Brick_Parameters = {}

        for i in range(len(Names_Chips)):
            # Create SiO2 Layer
            if PAD_Dist[i] > 0:
                Brick_Parameters["Brick Lenght Max"] = PAD_Dist[i] + PAD_Length 
                Brick_Parameters["Brick Lenght Min"] = PAD_Dist[i] - PAD_Length - 50
            else:
                Brick_Parameters["Brick Lenght Max"] = PAD_Dist[i] - PAD_Length 
                Brick_Parameters["Brick Lenght Min"] = PAD_Dist[i] + PAD_Length + 50
            Brick_Parameters["Brick Width Max"] = max(Dist) + PAD_Width*3
            Brick_Parameters["Brick Width Min"] = min(Dist) - PAD_Width*3
            Brick_Parameters["Brick Hight Max"] = 0
            Brick_Parameters["Brick Hight Min"] = -19.92
            Brick_Parameters["Brick Name"] = "SiO2_Layer" + Component_Name[i]
            Brick_Parameters["Component Name"] = Names_Chips[i]
            Brick_Parameters["Material"] = Material_iso
            self.Brick(Brick_Parameters)

            # Create Substrate
            if PAD_Dist[i] > 0:
                Brick_Parameters["Brick Lenght Max"] = PAD_Dist[i] + PAD_Length 
                Brick_Parameters["Brick Lenght Min"] = PAD_Dist[i] - PAD_Length - 50
            else:
                Brick_Parameters["Brick Lenght Max"] = PAD_Dist[i] - PAD_Length 
                Brick_Parameters["Brick Lenght Min"] = PAD_Dist[i] + PAD_Length + 50
            Brick_Parameters["Brick Width Max"] = max(Dist) + PAD_Width*3
            Brick_Parameters["Brick Width Min"] = min(Dist) - PAD_Width*3
            Brick_Parameters["Brick Hight Max"] = -19.92
            Brick_Parameters["Brick Hight Min"] = -298.86
            Brick_Parameters["Brick Name"] = "Substrate_Chip" + Component_Name[i]
            Brick_Parameters["Component Name"] = Names_Chips[i]
            Brick_Parameters["Material"] = Material_Substrate
            self.Brick(Brick_Parameters)



        # Create GND PADS Bondwires
        Names_Wires = ["Wire_L_GND_1", "Wire_L_GND_2", "Wire_Mid_GND_1", "Wire_Mid_GND_2", "Wire_R_GND_1", "Wire_R_GND_2"]
        Dist2 = [-225, -175, -25, 25, 175, 225]

        for k in range(len(Dist2)):

            Parameters_Bondwire = {}
            Parameters_Bondwire['X1'] = PAD_Dist[0]/2 
            Parameters_Bondwire['Y1'] = Dist2[k]
            Parameters_Bondwire['Z1'] = PAD_Thickness
            Parameters_Bondwire['X2'] = PAD_Dist[1]/2 
            Parameters_Bondwire['Y2'] = Dist2[k]
            Parameters_Bondwire['Z2'] = PAD_Thickness


            ParametersWire = {}
            ParametersWire["Name Wire"] = Names_Wires[k]
            ParametersWire["Coordinates"] = Parameters_Bondwire
            ParametersWire["Height"] = Bondwire_Height
            ParametersWire["Radius"] = Bondwire_Radius
            ParametersWire["Bondwire Type"] = "Spline"
            ParametersWire["Center Position"] = 0.5
            ParametersWire["Material"] = "Al"
            ParametersWire["SolidWireModel"] = True
            ParametersWire["Termination"] = "extended"
            ParametersWire["NameFolder"] = Names_Wires[k] + "_BondWire"
            self.BondWire(ParametersWire)
            
            Parameters_toSolid = {}
            Parameters_toSolid["Solid Name"] = Names_Wires[k]
            Parameters_toSolid["Curve Name"] = Names_Wires[k]
            Parameters_toSolid["Name Folder"] =  Names_Wires[k] + "_BondWire"
            Parameters_toSolid["Material"] =  "Al"
            self.ToSolid(Parameters_toSolid)
            

        
        # Create Signal PADS Bondwires
        Names_Wires = ["Wire_L_Sig", "Wire_R_Sig"]
        Dist2 = [-100, 100]


        for k in range(len(Dist2)):
            Parameters_Bondwire = {}
            Parameters_Bondwire['X1'] = PAD_Dist[0]/2 
            Parameters_Bondwire['Y1'] = Dist2[k]
            Parameters_Bondwire['Z1'] = PAD_Thickness
            Parameters_Bondwire['X2'] = PAD_Dist[1]/2 
            Parameters_Bondwire['Y2'] = Dist2[k]
            Parameters_Bondwire['Z2'] = PAD_Thickness

            
            ParametersWire = {}
            ParametersWire["Name Wire"] = Names_Wires[k]
            ParametersWire["Coordinates"] = Parameters_Bondwire
            ParametersWire["Height"] = Bondwire_Height
            ParametersWire["Radius"] = Bondwire_Radius
            ParametersWire["Bondwire Type"] = "Spline"
            ParametersWire["Center Position"] = 0.5
            ParametersWire["Material"] = "Al"
            ParametersWire["SolidWireModel"] = True
            ParametersWire["Termination"] = "extended"
            ParametersWire["NameFolder"] = Names_Wires[k] + "_BondWire"
            self.BondWire(ParametersWire)

            Parameters_toSolid = {}
            Parameters_toSolid["Solid Name"] = Names_Wires[k]
            Parameters_toSolid["Curve Name"] = Names_Wires[k]
            Parameters_toSolid["Name Folder"] =  Names_Wires[k] + "_BondWire"
            Parameters_toSolid["Material"] =  "Al"
            self.ToSolid(Parameters_toSolid)
    


        

        # Create Top Plate for Bond Wires and cladding
        if Parameters["Floating Shield"] == True:

            # Set glue thickness and floating shield thickness
            Glue_Thickness = Parameters["Glue Thickness"]
            FloatingShieldThickness = Parameters["Floating Shield Thickness"]

            # PAD_Thickness = 2.8
            Name_Shield = "Floating_Shield"
            Component_Name = "Floating_Shield"
            Material = "Au"
            Material_Clad = "DAF_Glue"


            # Create cladding
            Brick_Parameters = {}
            Brick_Parameters["Brick Lenght Max"] = PAD_Dist[0]/2 - PAD_Length
            Brick_Parameters["Brick Lenght Min"] = PAD_Dist[1]/2 + PAD_Length
            Brick_Parameters["Brick Width Max"] = max(Dist) + PAD_Width*2
            Brick_Parameters["Brick Width Min"] = min(Dist) - PAD_Width*2
            Brick_Parameters["Brick Hight Max"] = 2 * (PAD_Thickness + Bondwire_Height + Bondwire_Radius*2) + Glue_Thickness*2
            Brick_Parameters["Brick Hight Min"] = 0
            Brick_Parameters["Brick Name"] = "Floating_Shield_Clad"
            Brick_Parameters["Component Name"] = Component_Name
            Brick_Parameters["Material"] = Material_Clad
            self.Brick(Brick_Parameters)



            # Create squere floating shield
            Brick_Parameters["Brick Lenght Max"] = PAD_Dist[0]/2 - PAD_Length
            Brick_Parameters["Brick Lenght Min"] = PAD_Dist[1]/2 + PAD_Length
            Brick_Parameters["Brick Width Max"] = max(Dist) + PAD_Width*2
            Brick_Parameters["Brick Width Min"] = min(Dist) - PAD_Width*2
            Brick_Parameters["Brick Hight Max"] = 2* (PAD_Thickness + Bondwire_Height + Bondwire_Radius*2) + Glue_Thickness*2 + FloatingShieldThickness*2
            Brick_Parameters["Brick Hight Min"] = 2* (PAD_Thickness + Bondwire_Height + Bondwire_Radius*2) + Glue_Thickness*2
            Brick_Parameters["Brick Name"] = Name_Shield
            Brick_Parameters["Component Name"] = Component_Name
            Brick_Parameters["Material"] = Material
            self.Brick(Brick_Parameters)
        
        else: 
            pass


        
        # Pick Faces for Input Waveguide Port 1
        Parameters_Port = {}
        Facer_input = []
        Parameters_Port["Option"] = "Face"

        for i in range(len(Name)):
            Facer_input.append(Name[i] + ":" + Name[i] + "_input")

        # pick input faces
        for i in range(len(Name)):
            Parameters_Port["Face Number"] = Facet_Number[0]
            Parameters_Port["Object"] = Facer_input[i]
            self.Pick(Parameters_Port)

        # Create Port Input
        Parameters_Port = {}
        Parameters_Port["Port Number"] = 1
        Parameters_Port["Coordinates"] = "Picks"
        Parameters_Port["Orientation"] = "Positive"
        Parameters_Port["Span"] = [3, SpanY, SpanZ]
        Parameters_Port["Picked Port Number"] = [1,1,1,2,2,2]
        Parameters_Port["Number of picks"] = 3
        Parameters_Port["Picked Port Polarity"] = ["negative", "positive", "negative", "negative", "positive","negative"]
        Parameters_Port["Picked Component Name"] = [ "GND_L:GND_L_input", "Sig_L:Sig_L_input", "GND_Mid:GND_Mid_input", "GND_Mid:GND_Mid_input", "Sig_R:Sig_R_input", "GND_R:GND_R_input"]
        Parameters_Port["Face Number"] = [Facet_Number[0], Facet_Number[0], Facet_Number[0], Facet_Number[0], Facet_Number[0], Facet_Number[0]]
        self.WaveguidePortWithPins(Parameters_Port)



        # Pick Faces for Input Waveguide Port 2
        Parameters_Port = {}
        Facer_output = []
        Parameters_Port["Option"] = "Face"
        for i in range(len(Name)):
            Facer_output.append(Name[i] + ":" + Name[i] + "_output")

        for i in range(len(Name)):
            Parameters_Port["Face Number"] = Facet_Number[1]
            Parameters_Port["Object"] = Facer_output[i]
            self.Pick(Parameters_Port)

        #Create Port Output
        Parameters_Port = {}
        Parameters_Port["Port Number"] = 2
        Parameters_Port["Coordinates"] = "Picks"
        Parameters_Port["Orientation"] = "Positive"
        Parameters_Port["Span"] = [3, SpanY, SpanZ]
        Parameters_Port["Picked Port Number"] = [1,1,1,2,2,2]
        Parameters_Port["Number of picks"] = 5
        Parameters_Port["Picked Port Polarity"] = ["negative", "positive", "negative", "negative", "positive","negative"]
        Parameters_Port["Picked Component Name"] = [ "GND_L:GND_L_output", "Sig_L:Sig_L_output", "GND_Mid:GND_Mid_output", "GND_Mid:GND_Mid_output", "Sig_R:Sig_R_output", "GND_R:GND_R_output"]
        Parameters_Port["Face Number"] = [Facet_Number[1], Facet_Number[1], Facet_Number[1], Facet_Number[1], Facet_Number[1], Facet_Number[1]]
        self.WaveguidePortWithPins(Parameters_Port)


        # Set FDTD Time domain solver 
        
        # Set Time Solver
        Parameters_Solver = {}
        Parameters_Solver["Accuracy"] = Solver_Accuracy
        Parameters_Solver["Auto Impedance"] = Solver_Impedance_Status
        Parameters_Solver["Impedance"] = Solver_Impedance_Value
        Parameters_Solver["Source Port"] = Solver_Source_Port
        Parameters_Solver["Solver Mesh Type"] = Solver_Mesh_Type
        Parameters_Solver["Caclculate Modes Only"] = Solver_Calc_Modes
        self.setTimeSolver(Parameters_Solver)




    # def BondWire(self, NameWire, Coordinates, Height, Radius, BondwireType = "Spline", CenterPosition = 0.5, alpha = None, beta = None, Material = None, SolidWireModel = True, Termination = None, NameFolder = None):
    def BondWire(self, Parameters):
        """
        Create Bond Wire.

        Parameters
        ----------
        Parameters["Name Wire"] : str
            Set the name of the Bondwire.
        Parameters["Coordinates"] : dict of int/float 
            Dictionary with Coordinates in X,Y,Z plane to create the Bondwire:  
            For Example 
                Points = {}
                Points['X1'] = 0
                Points['Y1'] = 0
                Points['Z1'] = 0
                Points['X2'] = 5
                Points['Y2'] = 5
                Points['Z2'] = 0
                Parameters["Coordinates"] = Points
        Parameters["Height"] : int/float
            Set the hight of the middle point of the Bondwire.
        Parameters["Radius"] : str
            Set the radius of the bond wire. Bond wire is an cylinder type of object.
        Parameters["Bondwire Type"] : str, optional
            Set the type of bond wire. Defaults to "Spline". Other inputs are:
                                                                                "Spline" 
                                                                                "JEDEC4"
                                                                                "JEDEC5".
        Parameters["Center Position"] : int/float
            The center Position of the Height. This can be moved to make 
            an object that dont have the top height in the middle. Defaults to 0.5..
        Parameters["alpha"] : int/float, optional
            When using "JEDEC5" an alpha parameter need to be defined. See CST documentation!! Defaults to None..
        Parameters["beta"] : int/float, optional
            When using "JEDEC5" an betta parameter need to be defined. See CST documentation!!. Defaults to None..
        Parameters["Material"] : str, optional 
            Material for of the Bond wire. For now you need to load the material 
            in your simulation and then use this function. Otherwise the material 
            will not be found. Defaults to "PCE".
        Parameters["SolidWireModel"] : boolen, optional
            This option will turn the bondiwre into solid object. Defaults to True..
        Parameters["Termination"] : str, optional
            How the Bondwire will be temrinated. Defaults to None. Options are:
                                                                            "natural" 
                                                                            "rounded" 
                                                                            "extended"
        Parameters["NameFolder"] : str, optional
            The name of the folder. Defaults to name of the wire.

        Raises
        ------
        ValueError
            Check if some parameters are set correctly.

        Returns
        -------
        None.

        """
        
        # Set Parameters
        NameWire = Parameters["Name Wire"]
        Coordinates = Parameters["Coordinates"] 
        Height =  Parameters["Height"]
        Radius = Parameters["Radius"]
        BondwireType = Parameters["Bondwire Type"]
        CenterPosition = Parameters["Center Position"]
        Material = Parameters["Material"]
        SolidWireModel = Parameters["SolidWireModel"]
        Termination = Parameters["Termination"]
        NameFolder = Parameters["NameFolder"]

        # Check the BondwireType
        if BondwireType in ["Spline", "JEDEC4", "JEDEC5"]:
            BondwireType = BondwireType
            if BondwireType == "JEDEC5":
                alpha = Parameters["alpha"]
                beta = Parameters["beta"]
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
            Material = Material


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



    # def CurveWire(self, NameWire, Radius, Points = None, Material = None, SolidWireModel = True, Termination = None, NameFolder = None, CurveFolderName = None, CurveName = None):
    def CurveWire(self, Parameters):
        """Create curve wire based coordinates parameters points 

        Args: Parameters (dict) Dictionary with parameters needed for the function.
                Parameters["Name Wire"] : (str) Name of the wire
                Parameters["Radius"] : (int/float) Radius of the curve
                Parameters["Points"] : (dict, optional) Dictionary of X and  Y points for the curve. Defaults to None.
                Parameters["Material"] : (str, optional) Material for of the Bond wire. For now you need to load the material 
                                    in your simulation and then use this function. Otherwise the material 
                                    will not be found. Defaults to "PCE".
                Parameters["SolidWireModel"] : (bool, optional) This option will turn the Curve Wire into solid object. Defaults to True.
                Parameters["Termination"] : (str, optional) How the Bondwire will be temrinated. Defaults to None. Options are:
                                            "natural" 
                                            "rounded" 
                                            "extended"
                Parameters["NameFolder"] : (str, optional) The name of the bondwire folder. Defaults to name of the wire.
                Parameters["CurveFolderName"] : (str, optional) Name of the Curve folder name. Defaults to None.
                Parameters["CurveName"] : (str, optional) Curve name. Defaults to None.

        Raises:
            ValueError: Error massage
            ValueError: Error massage
            ValueError: Error massage

        Returns:
            str: String with VBA Code 

        """
        # Set Parameters
        NameWire = Parameters["Name Wire"]
        Radius = Parameters["Radius"]
        Points = Parameters["Points"]
        Material = Parameters["Material"] 
        SolidWireModel = Parameters["SolidWireModel"]
        Termination = Parameters["Termination"]
        NameFolder = Parameters["NameFolder"] 
        CurveFolderName = Parameters["CurveFolderName"]
        CurveName = Parameters["CurveName"]
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
            Material = Material

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
        """
        Create Cylinder object in CST enviroment.

        Parameters
        ----------
        Parameters["Cylinder Name"] : str
            Set the name of the cylinder.
        Parameters["Component Name"] : str
            Set the Name of the cylinder component in the component tree.
        Parameters["Material"] : str, optional
            Material for of the Bond wire. For now you need to load the material 
            in your simulation and then use this function.
        Parameters["Outer Radius"] : int/float
            Set the cylinder outer radius.
        Parameters["Inner Radius"] : int/float
            Set the cylinder inner radius
        Parameters["Orentation Axis"] : str
            Cylunder Orientation axis can be "X", "Y" or "Z".
            If Parameters["Orentation Axis"] = "X" then the following parameters are needed :
                Parameters["X min"] : int/float
                    Set the X-min parameter
                Parameters["X max"] : int/float
                    Set the X_max parameter
                Parameters["Z center"] : int/float
                    Set the Z-center parameter
                Parameters["Y center"] : int/float
                    Set the Y_center parameter
            If Parameters["Orentation Axis"] = "Y" then the following parameters are needed :
                Parameters["Y min"] : int/float
                    Set the Y-min parameter
                Parameters["Y max"] : int/float
                    Set the Y-max parameter
                Parameters["Z center"] : int/float
                    Set the Z-center parameter
                Parameters["X center"] : int/float
                    Set the X-center parameter
            If Parameters["Orentation Axis"] = "Z" then the following parameters are needed :
                Parameters["Z min"] : int/float
                    Set the Z-min parameter
                Parameters["Z max"] : int/float
                    Set the Z-max parameter
                Parameters["X center"] : int/float
                    Set the X-center parameter
                Parameters["Y center"] : int/float
                    Set the Y-center parameter

        Raises
        ------
        ValueError
            Values Errors.

        Returns
        -------
        None.

        """
    

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
        """
        Create Conus object in CST enviroment.

        Parameters
        ----------
        Parameters["Conus Name"] : str
            Set the Name of the conus.
        Parameters["Component Name"] : str
            Set the  Name of the conus component in the component tree.
        Parameters["Material"] : str
            Set the material for of the Conus. 
        Parameters["Top Radius"] : int/float
            Set the conus top radius.
        Parameters["Bottom Radius"] : int/float
            Set the conus bottom radius.
        Parameters["Orentation Axis"] : str
            Set the conus Orientation axis can be "X", "Y" or "Z".
            If Parameters["Orentation Axis"] = "X" then the following parameters are needed :
                Parameters["X min"] : int/float 
                    Set the X-min parameter
                Parameters["X max"] : int/float
                    Set the X-max parameter
                Parameters["Z center"] : int/float
                    Set the Z-center parameter
                Parameters["Y center"] : int/float 
                    Set the Y-center parameter
            If Parameters["Orentation Axis"] = "Y" then the following parameters are needed :
                Parameters["Y min"] : int/float
                    Set the Y-min parameter
                Parameters["Y max"] : int/float
                    Set the Y-max parameter
                Parameters["Z center"] : int/float
                    Set the Z-center parameter
                Parameters["X center"] : int/float
                    Set the X-center parameter
            If Parameters["Orentation Axis"] = "Z" then the following parameters are needed :
                Parameters["Z min"] : int/float
                    Set the Z-min parameter
                Parameters["Z max"] : int/float
                    Set the Z-max parameter
                Parameters["X center"] : int/float 
                    Set the X-center parameter
                Parameters["Y center"] : int/float
                    Set the Y-center parameter

        Returns
        -------
        None.

        """
        
        
        Name = Parameters["Conus Name"]
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
        """
        Create GGB Probes in CSt enviroment.

        Parameters
        ----------
        Parameters["Component Name"] : str
            Set the Name of the Component.
        Parameters["Orientation Angle"] : int/float
            Set the angle of tilting the GGB Probes

        Returns
        -------
        None.

        """
        
        
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
        DictComponent["Conus Name"] = "Conus_Probe_Sig"
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
        DictComponent["Conus Name"] = "Conus_Probe_Sig_Tip"
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
        # DictComponent["Position X"] = 25/2
        # DictComponent["Position Y"] = -186.1
        # DictComponent["Position Z"] = 0
        DictComponent["Position"] = [25/2, -186.1, 0]
        # DictComponent["Structure Type"] = "Shape"
        DictComponent["Object to rotate"] = "Shape"
        
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
        # DictComponent["Position X"] = -25/2
        # DictComponent["Position Y"] = 186.1
        # DictComponent["Position Z"] = 0
        DictComponent["Position"] = [-25/2, 186.1, 0]
        # DictComponent["Structure Type"] = "Shape"
        DictComponent["Object to rotate"] = "Shape"
        
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
        # Pos = {}
        # Pos["X"] = 0
        # Pos["Y"] = 0
        # Pos["Z"] = 0
        DictComponent["Center Positions"] = [0, 0, 0]
        # Scale = {}
        # Scale["X"] = 1.5
        # Scale["Y"] = 2
        # Scale["Z"] = 1
        DictComponent["Scale Factors"] = [1.5, 2, 1]
        DictComponent["Object to rotate"] = "Shape"

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
        DictComponent["Conus Name"] = "Conus_Probe_Sig"
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
        # DictComponent["Position X"] = 25/2
        # DictComponent["Position Y"] = -186.1
        # DictComponent["Position Z"] = 0
        DictComponent["Position"] = [25/2, -186.1, 0]
        # DictComponent["Structure Type"] = "Shape"
        DictComponent["Object to rotate"] = "Shape"
        
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
        # DictComponent["Position X"] = -25/2
        # DictComponent["Position Y"] = 186.1
        # DictComponent["Position Z"] = 0
        DictComponent["Position"] = [-25/2, 186.1, 0]
        # DictComponent["Structure Type"] = "Shape"
        DictComponent["Object to rotate"] = "Shape"
        
        self.Translation(DictComponent)
        
        
        # Translate probe to -30 degree
        DictComponent = {}
        DictComponent["Translate Type"] = "Rotate"
        DictComponent["Name Object"] = ComponentName
        DictComponent["Origin"] = "Free"
        # Pos["X"] = 0
        # Pos["Y"] = 0
        # Pos["Z"] = 0
        DictComponent["Center Positions"] = [0, 0, 0]
        # angle = {}
        # angle["X"] = 0
        # angle["Y"] = angleOrientation
        # angle["Z"] = 0
        DictComponent["Angle Values"] = [0, angleOrientation, 0 ]
        DictComponent["Object to rotate"] = "Shape"
        
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
            # Pos["X"] = 0
            # Pos["Y"] = 0
            # Pos["Z"] = 0
            DictComponent["Center Positions"] = [0, 0, 0]
            # angle = {}
            # angle["X"] = 0
            # angle["Y"] = angleCutPlate
            # angle["Z"] = 0
            DictComponent["Angle Values"] = [0, angleCutPlate, 0]
            DictComponent["Object to rotate"] = "Shape"
            
            
            self.Translation(DictComponent)
            

            DictComponent = {}
            DictComponent["Translate Type"] = "Translate"
            DictComponent["Name Object"] = f"{ComponentName}:Cut_Plate"
            # DictComponent["Position X"] = 0
            # DictComponent["Position Y"] = 0
            # DictComponent["Position Z"] = -144
            DictComponent["Position"] = [0, 0, -144]
            # DictComponent["Structure Type"] = "Shape"
            DictComponent["Object to rotate"] = "Shape"
            
            self.Translation(DictComponent)
            
            DictComponent = {}
            DictComponent["Name Structure to Cut"]  = f"{ComponentName}:Cut_Plate"
            DictComponent["Name Cut Structure"] = ObjectsToCut[i]
            self.Cut_structures(DictComponent)
            
            
        # Tranlate Probe to Tips Z = 0 Positiopn 
        DictComponent = {}
        DictComponent["Translate Type"] = "Rotate"
        DictComponent["Name Object"] = ComponentName
        # Pos["X"] = 0
        # Pos["Y"] = 0
        # Pos["Z"] = 0
        DictComponent["Center Positions"] = [0, 0, 0]
        # angle = {}
        # angle["X"] = 0
        # angle["Y"] = angleOrientation
        # angle["Z"] = 0
        DictComponent["Angle Values"] = [0, angleOrientation, 0]
        DictComponent["Object to rotate"] = "Shape"
        
        self.Translation(DictComponent)
        
        DictComponent = {}
        DictComponent["Translate Type"] = "Translate"
        DictComponent["Name Object"] = ComponentName
        # DictComponent["Position X"] = 0
        # DictComponent["Position Y"] = 0
        # DictComponent["Position Z"] = 124.5
        DictComponent["Position"] = [0, 0, 124.5]
        DictComponent["Object to rotate"] = "Shape"
        
        self.Translation(DictComponent)
        
       




############################################################################
# Transform to Solid
############################################################################

    # def ToSolid(self, SolidName, CurveName = "Polygon", NameFolder = None, Material = None ):
    def ToSolid(self, Parameters):
        """
        This function transfer an function or polynom to solid object.

        Parameters
        ----------
        Parameters["Solid Name"] : str
            Set the name of the solid object.
        Parameters["Curve Name"] : str
            Name of the curve. Defaults to "Polygon".
        Parameters["Name Folder"] : str
            Set the name of the folder. Defaults to Curve name.
       Parameters["Material"] : str
            Material for of the Bond wire. For now you need to load the material 
            in your simulation and then use this function. Otherwise the material 
            will not be found. Defaults to Aluminum.

        Returns
        -------
        None.

        """
        
        # Set Parameters 
        SolidName = Parameters["Solid Name"]
        CurveName = Parameters["Curve Name"]
        NameFolder = Parameters["Name Folder"]
        Material = Parameters["Material"]       

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
        """
        This function transfer curve to solid object.

        Parameters
        ----------
        Parameters["Name"] : str
            Set the Name of the solid object.
        Parameters["Component Name"] : str
            Set the name of the component in the component tree.
        Parameters["Material"] : str
            Set the material for the solid.
        Parameters["Thickness"] : int/float
            Set the thickness of the solid curve.
        Parameters["Angle"] : int/float
            Set the angle.
        Parameters["Curve Folder"] : str
            Set the curve folder. First you need to crate curve and then give here the correct name of the curve folder.
        Parameters["Curve Name"] : str
            Set the curve folder. First you need to crate curve and then give here the correct name of the curve.

        Returns
        -------
        None.

        """
        
        
        Name = Parameters["Name"] 
        NameComponent = Parameters["Component Name"]
        Material = Parameters["Material"]
        Thickness = Parameters["Thickness"]
        CurveFolder = Parameters["Curve Folder"]
        CurveName = Parameters["Curve Name"]
        
        if "Angle" not in Parameters.keys():
            Angle = 0
        else:
            Angle = Parameters["Angle"]
        
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
                
        
        

    # def RibWaveguide_ToSolid(self, SolidName, WaveguideName = "Rib_Waveguide", WG_Hight = None, Angle = None, NameFolder = None, Material = None, WGFolderName = None, WGName = None ):
    def RibWaveguide_ToSolid(self, Parameters):
        """
        This is ToSolid function that will allow the use to create the RibWaveguide. 
        TODO: Murrge the TOSolid and RibWaveguideToSolid functions later on. 

        Parameters
        ----------
        Parameters["Waveguide Name"] : TYPE
            Set the Waveguide Name. Defaults to "Rib_Waveguide".
        Parameters["Wavaguide Hight"] : int/float
            Set the hight of the waveguide. Defaults to None.
        Parameters["Angle"] : int/float
            Set the side angle of the waveguide. Defaults to None.
        Parameters["Name Folder"] : str
            Set the folder name. Defaults to None.
        Parameters["Material"] : str
            Material for of the Bond wire. For now you need to load the material 
            in your simulation and then use this function. Otherwise the material 
            will not be found. Defaults to None.
        Parameters["Waveguide Folder Name"] : str
            Set the name of the folder where the poligon 3D or 2D is created. Defaults to None.
        Parameters["Waveguide Name"] : TYPE
            DESCRIPTION.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        # Parameters set
        WaveguideName = Parameters["Waveguide Name"]
        WG_Hight = Parameters["Waveguide Hight"]
        Angle = Parameters["Angle"]
        NameFolder = Parameters["Name Folder"]
        Material = Parameters["Material"]
        WGFolderName = Parameters["Waveguide Folder Name"] 
        WGName = Parameters["Waveguide Name"]
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
        """
        Set the Waveguide Port.

        Parameters
        ----------
        Parameters["Orientation"] : str
            Set the Port Orientation it can be "Positive" or "Negative". For 
            this function an 2 ports will be defined so please give an array with two Oriantations
            like  Parameters["Orientation"] = ["Positive", "Positive"].
        Parameters["Coordinates"] : str
            Set the Coordinates type, "Picks" is the best one!. It can be on of:
                                                                                "Free"
                                                                                "Full"
                                                                                "Picks"
        Parameters["Span"] : int/float array
            Set an array of port span [[Ymin, Ymax],[Zmin, Zmax]].
        Parameters["Potential"] : int array
            Set the potential of the ports. For example [1,2].
        Parameters["Port Number"] : int array
            Set the array with port numbers.

        Returns
        -------
        None.

        """
        
        Coordinates_Types = ["Free", "Full", "Picks"]

        # Parameters to determin port position 
        Orientation = Parameters["Orientation"]
        Coordinates = Parameters["Coordinates"]
        # Choose coordinates
        Span11 = Parameters["Span"][0][0]
        Span12 = Parameters["Span"][0][1]
        Span21 = Parameters["Span"][1][0]
        Span22 = Parameters["Span"][1][1]
        # Polarity = Parameters["Polarity"]
        # SolidName = Parameters["Solid Name"]
        # PickID = Parameters["Face ID"]
        #Return Ports in Dictionary for one Object 2 Ports can be define 
        Port = {}
        Port["1"] = None
        Port["2"] = None

        # check coordinates type
        if Coordinates not in Coordinates_Types:
            raise ValueError("The Coordinates type can be one of 'Free' , 'Full' or 'Picks'")  
        else:
            pass
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
        """
        Set the Waveguide Port for  

        Parameters
        ----------
        Parameters["Orientation"] : str
            Set the port Orientation it can be "Positive" or "Negative". For 
            this function an 2 ports will be defined so please give an array with two Oriantations
            like  Parameters["Orientation"] = ["Positive", "Positive"].
        Parameters["Coordinates"] : str
            Set the Coordinates type, "Picks" is the best one!. It can be on of:
                                                                                "Free"
                                                                                "Full"
                                                                                "Picks".
        Parameters["Span"] :  int/float array
            Set an array of port span [[Ymin, Ymax],[Zmin, Zmax]].
        Parameters["Potential"] : int array
            Set the potential of the ports. For example [1,2].
        Parameters["Port Number"] : int array
            Set the array with port numbers.
        Parameters["Polarity"] : str array 
            Set the Port Polarity it can be be "Positive" or "Negative". For 
            this function an 2 ports will be defined so please give an 
            array with two Polaritys. For example Parameters["Polarity"] = ["Positive", "Positive"]
        Parameters["Solid Name"] : str
            Set the name of the Object on witch the Waveguide port will be created. For example "WG:WG1".
        Parameters["Face ID"] : int array
            Set the ID of the two picked faces. For example Parameters["Face ID"] = [2,4]
        PickParams["Face Number"] : int array
            Set the number of the picked face of the structure

        Returns
        -------
        None.

        """
       
        
        Coordinates_Types = ["Free", "Full", "Picks"]
        # Parameters to determin port position 
        Orientation = Parameters["Orientation"]
        Coordinates = Parameters["Coordinates"]
        # Choose coordinates
        PortNumber = Parameters["Port Number"]
        Span = Parameters["Span"]
        set_number = Parameters["Picked Port Number"] 
        Number_of_picks = Parameters["Number of picks"]
        
        # Check the corrdinates 
        # check coordinates type
        if Coordinates not in Coordinates_Types:
            raise ValueError("The Coordinates type can be one of 'Free' , 'Full' or 'Picks'")  
        else:
            pass
        
        
        lines = []

        lines = []
        if Number_of_picks > 2:
            potential = Parameters["Picked Port Polarity"]
            marked_face = Parameters["Picked Component Name"]
            facet_Id = Parameters["Face Number"]
            for i in range(len(marked_face)):
                lines.append(f'.AddPotentialPicked "{set_number[i]}", "{potential[i]}", "{marked_face[i]}", "{facet_Id[i]}"')

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
        """
        Function to move the waveguide ports.

        Parameters
        ----------
        Parameters["Port Number"] : int
            Set the port number.
        Parameters["Distance"] : int/float
            Set the distance to move the waveguide.
        Parameters["Span"] : dict
            Set an Dictionary with array to set the range in "X", "Y" and "Z"
                Parameters["Span"][0][0] - Yrange min
                Parameters["Span"][0][1] - Yrange max
                Parameters["Span"][1][0] - Zrange min
                Parameters["Span"][1][1] - Zrange max.
        

        Returns
        -------
        None.

        """
        
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
        """
        Function to set discrete port

        Parameters
        ----------
        Parameters["Discrete Port Number"] : int
            Set the Port number.
        Parameters["Discrete Port Type"] : str
            Set the type of discrete port, can be:
                                               "Voltage"
                                               "S-Parameters"
                                               "Current".
        Parameters["Port Impedance"] : int float
            Set the Impedance of the port.
        Parameters["Port Voltage"] : int/float
            Set Port voltage.
        Parameters["Port Current"] : int/float
            Set Port current.
        Parameters["Port Radius"] : int/float
            Set Port radius.
        Parameters["Discrete Port Coordinates"] : dict
            Set an dictionary with port coordiantes.
                Parameters["Discrete Port Coordinates"]["X1"]
                Parameters["Discrete Port Coordinates"]["Y1"]
                Parameters["Discrete Port Coordinates"]["Z1"]
                Parameters["Discrete Port Coordinates"]["X2"]
                Parameters["Discrete Port Coordinates"]["Y2"]
                Parameters["Discrete Port Coordinates"]["Z2"]

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        """

        PortNumber = Parameters["Discrete Port Number"]
        if Parameters["Discrete Port Type"] in ["Voltage", "S-Parameters", "Current"]:
            PortType = Parameters["Discrete Port Type"]
        else:
            raise ValueError("Port Type can be on of 'Voltage''S-Parameters' or 'Current'")
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
        """
        Pick function

        Parameters
        ----------
        PicParams["Option"] : str
            Set the Pcik option. For now only "Centerpoint" and "Face" can be choosen.
        PicParams["Object"] : str
            Set the name of the object on witch the pick will be executed..
        PicParams["Face Number"] : int
            Set the ID of the picked Face for example.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        obj_list = ["Centerpoint", "Face"]
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
        """
        Clear all pciks

        Returns
        -------
        None.

        """
        
        vba_code = f"""
                Pick.ClearAllPicks
                """ 
        self.prj.model3d.add_to_history("delete pick", vba_code)    





                
                    


############################################################################
# Solvers
############################################################################

    def setTimeSolver(self, Parameters):
        """
        Set Time solver Parameters

        Parameters
        ----------
        Parameters["Accuracy"] : int
            Set the accury of the solver. It can be 10
                                                    15
                                                    20
                                                    25
                                                    30
                                                    35
                                                    40
        Parameters["Caclculate Modes Only"] : boolen
            Set to True if you want to calculate the the Ports modes only.
            False to calculate the hole structure.
        Parameters["Auto Impedance"] : boolen
            Set Port Impedance. True if you want to set manually
            False otherwise.
        Parameters["Impedance"] : int/float
            Set the Port Impedance.
        Parameters["Source Port"] : int
            Set the Source Port. For example 1 or 2.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

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
            self.prj.model3d.add_to_history(f"set {MeshType} time solver", vba_code)  


        elif MeshType == "TLM":
            vba_code = f"""
                        Mesh.SetCreator "High Frequency"
                        With Solver
                        .Reset
                        .Method "Hexahedral TLM"
                        .SteadyStateLimit "-{Accuracy}"
                        .StimulationPort "{Source}"
                        .StimulationMode "All"
                        .MeshAdaption "False"
                        .SParaSymmetry "False"
                        .AutoNormImpedance "True"
                        .NormingImpedance "{Impedance}"
                        .StoreTDResultsInCache  "False"
                        .RunDiscretizerOnly "False"
                        .SuperimposePLWExcitation "False"
                        End With
                        """
            self.prj.model3d.add_to_history(f"set {MeshType} time solver", vba_code) 
        else:
            raise ValueError("Mesh Type is not recognized, only 'TLM' or 'FIT' are supported! ")






    def setFreqSolver(self):
        """
        Set Frequency solver Parameters
        TODO: This function is hard codet and need to be worked on. 
        You cannot set parameters!

        Returns
        -------
        None.

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
        """
        Start Time Simulation

        Returns
        -------
        None.

        """
        
        vba_code = f"""
                    With Solver
                    .Start
                    End With
                    """
        self.prj.model3d.add_to_history("start time solver", vba_code)



    # Set Optical Simulation Domain 
    def setOpticalSimulationProperties(self, Parameters):
        """
        Function to set an Optical Simulation Properties. 
        TODO : The function need to be expands so more parameters such as mesh properties can be added!

        Parameters
        ----------
        Parameters['Dimensions'] : str
            Set unit Leghts.
        Parameters['Frequency'] : str
            Set frequency Units.
        Parameters['Time'] : str
            Set the temperature Units.
        Parameters['Temperature'] : str
            Set the temperature Units.
        Parameters["Type Background"]  : str
            Set the background Type.
        Parameters["Xmin Background"] : int/float
            Set the backgroun X-min span.
        Parameters["Xmax Background"] : int/float
            Set the backgroun X-max span.
        Parameters["Ymin Background"] : int/float
            Set the backgroun Y-min span.
        Parameters["Ymax Background"] : int/float
            Set the backgroun Y-max span.
        Parameters["Zmin Background"] : int/float
            Set the backgroun Z-min span.
        Parameters["Zmax Background"] : int/float
            Set the backgroun Z-max span.
        Parameters["Min Wavelength"] : int/float
            Set the min wavelenght
        Parameters["Max Wavelength"] : int/float
            Set the max wavelenght
        Parameters["Xmin Boundary"] : int/floar
            Set the X-min boundary
        Parameters["Xmax Boundary"] : int/floar
            Set the X-max boundary
        Parameters["Ymin Boundary"] : int/floar
            Set the Y-min boundary
        Parameters["Ymax Boundary"] : int/floar
            Set the Y-max boundary
        Parameters["Zmin Boundary"] : int/floar
            Set the Z-min boundary
        Parameters["Zmax Boundary"] : int/floar
            Set the Z-max boundary
        Parameters["Xsymmetry Boundary"] : int/float
            Set symetrie to the X boundary
        Parameters["Ysymmetry Boundary"] : int/float
            Set symetrie to the Y boundary
        Parameters["Zsymmetry Boundary"] : int/float 
            Set symetrie to the Z boundary
            
        Returns
        -------
        None.

        """
        
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

        # # Mesh Settings
        # MeshType = Parameters["Mesh Type"] 
        # CellNear = Parameters["Mesh Cells Near Object"] 
        # CellFar = Parameters["Mesh Cells far Object"] 


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
                    End With

                    With Mesh
                        .MeshType "PBA"
                    End With

                ChangeSolverType("HF Time Domain")
                """
        self.prj.model3d.add_to_history(f"set optical simulation", vba_code)





      
    # Set Optical Simulation Domain 
    def setElectricalSimulationProperties(self, Parameters):
        """
        Set electrical simulation properties
        TODO : The function need to be expands so more parameters such as mesh properties can be added!

        Parameters
        ----------
        Parameters['Dimensions'] : str
            Set unit Leghts.
        Parameters['Frequency'] : str
            Set frequency Units.
        Parameters['Time'] : str
            Set the temperature Units.
        Parameters['Temperature'] : str
            Set the temperature Units.
        Parameters["Type Background"]  : str
            Set the background Type.
        Parameters["Xmin Background"] : int/float
            Set the backgroun X-min span.
        Parameters["Xmax Background"] : int/float
            Set the backgroun X-max span.
        Parameters["Ymin Background"] : int/float
            Set the backgroun Y-min span.
        Parameters["Ymax Background"] : int/float
            Set the backgroun Y-max span.
        Parameters["Zmin Background"] : int/float
            Set the backgroun Z-min span.
        Parameters["Zmax Background"] : int/float
            Set the backgroun Z-max span.
        Parameters["Min Frequency"] : int/float
            Set the min frequency
        Parameters["Max Frequency"] : int/float
            Set the max frequency
        Parameters["Xmin Boundary"] : int/floar
            Set the X-min boundary
        Parameters["Xmax Boundary"] : int/floar
            Set the X-max boundary
        Parameters["Ymin Boundary"] : int/floar
            Set the Y-min boundary
        Parameters["Ymax Boundary"] : int/floar
            Set the Y-max boundary
        Parameters["Zmin Boundary"] : int/floar
            Set the Z-min boundary
        Parameters["Zmax Boundary"] : int/floar
            Set the Z-max boundary
        Parameters["Xsymmetry Boundary"] : int/float
            Set symetrie to the X boundary
        Parameters["Ysymmetry Boundary"] : int/float
            Set symetrie to the Y boundary
        Parameters["Zsymmetry Boundary"] : int/float 
            Set symetrie to the Z boundary

        Returns
        -------
        None.

        """
       
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

        # # Mesh Settings
        # MeshType = Parameters["Mesh Type"] 
        # CellNear = Parameters["Mesh Cells Near Object"] 
        # CellFar = Parameters["Mesh Cells far Object"] 




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
                        .FrequencyRange "{FreqMin}", "{FreqMax}"
                    End With
                        Mesh.MinimumCurvatureRefinement "150"
                        With MeshSettings
                            .SetMeshType "HexTLM"
                            .Set "StepsPerWaveNear", "30"
                        .Set "StepsPerWaveFar", "20"
                            .Set "StepsPerBoxNear", "20"
                        .Set "StepsPerBoxFar", "20"
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



    def ChangeSolverType(self, Parameters):
        """
        Change solver type
        TODO: Include the rest of the solver types. For now only time is available!

        Parameters
        ----------
        Parameters['Type'] : str
            Solver type. Can be Parameters['Type'] = "Time".

        Returns
        -------
        None.

        """
       
        time = ["TIME", "time", "Time", "t", "T"]
        Type = Parameters["Type"]
        if Type in time:
            vba_code = f"""
                        ChangeSolverType "HF Time Domain"
                        """
            self.prj.model3d.add_to_history("change solver", vba_code)
        
        else: 
            pass
    



    def setDomainSolverType(self, Parameters):
        """
        Set Domein solver type


        Parameters
        ----------
        Parameters['Domain'] : str
            Set the solver domain. Can be:
                                        "Time"
                                        "Freq"
                                        "EigenMode"
                                        "Integral"
                                        "Asymptotic"
                                        "Multilayer".

        Returns
        -------
        None.

        """
        

        #Set Parameters
        Domain = Parameters["Domain"]

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
        """
        Set Monitor

        Parameters
        ----------
        Parameters["Monitor Wavelength"] : int/float
            Set the wavelength. For this function and only for this function you need to give the exact number. So 1.55 um will be Parameters["Wavelength"]  = 1.55e-6
        Parameters["Monitor Frequency"] : int/float
            Set Frequency of the Monitor.
        Parameters["Domain"] : str
            Set the domain can be "Frequency" or "Wavelength".
        Parameters["Monitor Type"] : str
            Set the monitor Type. Can be one of :
            "Efield", "Hfield", "Surfacecurrent", "Powerflow", "Current",
            "Powerloss", "Eenergy", "Elossdens", "Lossdens", "Henergy", 
            "Farfield", "Fieldsource", "Spacecharge", "ParticleCurrentDensity", "Electrondensity" .

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

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
        """
        Sets the type of the mesh. The user can define the mesh cells per wavelength near and far
        from the simulations object.

        Parameters
        ----------
        Parameters["Mesh Type"] : TYPE
            Set the type of the Mesh. You can choose between:  
                                    PBA - Hexahedral mesh with Perfect Boundary Approximation
                                    HexahedralTLM
                                    CFD
                                    CFDNew
                                    Staircase - Hexahedral mesh with staircase cells
                                    Tetrahedral - Tetrahedral mesh
                                    Surface - Surface mesh
                                    SurfaceMLS - urface multi layer mesh
                                    Planar - Planar 2D mesh.
        Parameters["Mesh Cells Near Object"] : int
            Set the cells per Wavelength near the simulation object.
        Parameters["Mesh Cells far Object"] : int
            Set the Cells per Wavelength far from the simulation object

        Returns
        -------
        None.

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
        """
        Rotation Translation function

        Parameters
        ----------
        Parameters["Name Object"] : str
            Set the name of the Object that you want to rotate.
        Parameters["Angle X"] : int/float
            Set the X Angle of rotation.
        Parameters["Angle Y"] : int/float
            Set the Y Angle of rotation.
        Parameters["Angle Z"] : int/float
            Set the Z Angle of rotation.

        Returns
        -------
        None.

        """
        

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
        self.prj.model3d.add_to_history(f"rotational translation", vba_code)




    def Translation(self, Parameters):
        """
        Translation function

        Parameters
        ----------
        Parameters["Translate Type"] : str
            Set the type of translation. Can be:
                                            "Translate"
                                            "Scale"
                                            "Rotate".
        If Parameters["Translate Type"] = "Translate"
            Parameters["Name Object"] : str 
                Set the name of the Object that you want to rotate.
            Parameters["Object to rotate"] : str 
                Set the object you wanna rotate. If Solid
                                                Parameters["Object to rotate"] = "Shape", 
                                                If Curve 
                                                Parameters["Object to rotate"] = "Curve"
            Parameters["Position"] : int/float array 
                Set array with position parameters. For ["X", "Y" , "Z"] directions 
                For Example -> Parameters["Position"] = [1,2,3]      
        If Parameters["Translate Type"] = "Scale"
             Parameters["Name Object"] : str
                 Set the name of the Object that you want to rotate
             Parameters["Center Positions"] : int/float array
                 Set an array with Center position for scale operation. For ["X", "Y" , "Z"] directions 
                 For Example -> Parameters["Center Positions"] = [1,2,3]     
             Parameters["Scale Factors"] : int/float array 
                 Set array with Scale factor. For ["X", "Y" , "Z"] directions 
                 For Example -> Parameters["Scale Factors"] = [1,2,3]     
        If Parameters["Translate Type"] = "Rotate"
            Parameters["Name Object"] : str
                Set the name of the Object that you want to rotate
            Parameters["Object to rotate"] : str
                Set the object you wanna rotate. If Solid
                                                Parameters["Object to rotate"] = "Shape"
                                                If Curve 
                                                Parameters["Object to rotate"] = "Curve"
            Parameters["Center Positions"] : int/float array 
                Set the Center position for scale operation. For ["X", "Y" , "Z"] directions 
                For Example -> Parameters["Center Positions"] = [1,2,3]
            Parameters["Angle Values"] : int/float array 
                Set the Scale factor. For ["X", "Y" , "Z"] directions 
                For Example -> Parameters["Angle Values"] = [1,2,3]

        Returns
        -------
        None.

        """
        
        Type = Parameters["Translate Type"]
        Type_list = ["Translate", "Scale", "Rotate"]
        
        
        

        if Type in Type_list:
            if Type == "Translate":
                Name = Parameters["Name Object"]
                PosX = Parameters["Position"][0]
                PosY = Parameters["Position"][1]
                PosZ = Parameters["Position"][2]
                OBJ = Parameters["Object to rotate"]
                if OBJ != "Curve":
                    auto_destination = '     .AutoDestination "True"\n'
                else:
                    auto_destination = '\n'
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
                            {auto_destination}    
                            .Transform "{OBJ}", "Translate
                            End With
                            """
                self.prj.model3d.add_to_history(f"translation {Name}", vba_code)
                
            elif Type == "Scale":  
                Name = Parameters["Name Object"]
                CenterPos = Parameters["Center Positions"]
                ScaleFactor = Parameters["Scale Factors"]
                OBJ = Parameters["Object to rotate"]
                if OBJ != "Curve":
                    auto_destination = '     .AutoDestination "True"\n'
                else:
                    auto_destination = '\n'
                if "Origin" in Parameters.keys():
                    Origin = Parameters["Origin"]
                else:
                    Origin = "Free"
                vba_code2 = f"""                                        
                            With Transform 
                                 .Reset 
                                 .Name "{Name}" 
                                 .Origin "Free" 
                                 .Center "{CenterPos[0]}", "{CenterPos[1]}", "{CenterPos[2]}" 
                                 .ScaleFactor "{ScaleFactor[0]}", "{ScaleFactor[1]}", "{ScaleFactor[2]}" 
                                 .MultipleObjects "False" 
                                 .GroupObjects "False" 
                                 .Repetitions "1" 
                                 .MultipleSelection "False" 
                                 {auto_destination}    
                                 .Transform "{OBJ}", "Scale" 
                            End With
                         """
                self.prj.model3d.add_to_history(f"translation scale {Name}", vba_code2)
            elif Type == "Rotate":
                Name = Parameters["Name Object"]
                CenterPos = Parameters["Center Positions"]
                Angle = Parameters["Angle Values"]
                OBJ = Parameters["Object to rotate"]
                if OBJ != "Curve":
                    auto_destination = '     .AutoDestination "True"\n'
                else:
                    auto_destination = '\n'
                if "Origin" in Parameters.keys():
                    Origin = Parameters["Origin"]
                else:
                    Origin = "Free"

                vba_code3 = f"""
                With Transform
                    .Reset
                    .Name "{Name}"
                    .Origin "{Origin}"
                    .Center "{CenterPos[0]}", "{CenterPos[1]}", "{CenterPos[2]}"
                    .Angle "{Angle[0]}", "{Angle[1]}", "{Angle[2]}"
                    .MultipleObjects "False"
                    .GroupObjects "False"
                    .Repetitions "1"
                    .MultipleSelection "False"
                    {auto_destination}     
                    .Transform "{OBJ}", "Rotate"
                End With
                """
                self.prj.model3d.add_to_history(f"translation scale {Name}", vba_code3)
                        
        
    
    def Translation_Scale(self, Parameters):
        """
        Scale tranlation function

        Parameters
        ----------
        Parameters["Name Object"] : str
            Set the name of the Object that you want to rotate
        Parameters["Center Positions"] : dict
            Set an dictionary with Center position for scale operation
                Parameters["Center Positions"]["X"] : int/float 
                    Set Center position X
                Parameters["Center Positions"]["Y"] : int/float 
                    Set Center position Y
                Parameters["Center Positions"]["Z"] : int/float
                Set the Center position Z
        Parameters["Scale Factors"] : dict 
            Set the scale factor S
            Parameters["Scale Factors"]["X"] : int/float
                Set the scale factor X
            Parameters["Scale Factors"]["Y"] : int/float
                Set the scale factor Y
            Parameters["Scale Factors"]["Z"] : int/float
                Set the scale factor Z

        Returns
        -------
        None.

        """
        

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
        self.prj.model3d.add_to_history(f"translation scale {Name}", vba_code)
        
        
        
        
    def Translation_Rotation(self, Parameters):
        """
        Rotation tranlation function

        Parameters
        ----------
        Parameters["Name Object"] : str
            Set the name of the Object that you want to rotate
        Parameters["Center Positions"] : dict
            Set the Center position for scale operation
                Parameters["Center Positions"]["X"] : int/float 
                    Set the center position X
                Parameters["Center Positions"]["Y"] : int/float 
                    Set the center position Y
                Parameters["Center Positions"]["Z"] : int/float 
                    Set the center position Z
        Parameters["Angle Values"] : dict
            Set the scale factor
                Parameters["Angle Values"]["X"] : int/float 
                    Set the Angle Value X
                Parameters["Angle Values"]["Y"] : int/float 
                    Set the Angle Value Y
                Parameters["Angle Values"]["Z"] : int/float 
                    Set the Angle Value Z

        Returns
        -------
        None.

        """
        
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
        self.prj.model3d.add_to_history(f"translation scale {Name}", vba_code)
    



    def Cut_structures(self, Parameters):
        """
        Cut function

        Parameters
        ----------
        Parameters["Name Cut Structure"] : str
            Set the name of the structure to cut.
        Parameters["Name Structure to Cut"] : str
            Set the name Structure to cut the Parameters["Name Cut Structure"] to.

        Returns
        -------
        None.

        """
        
        CutElement = Parameters["Name Cut Structure"]
        CutterElement = Parameters["Name Structure to Cut"]
        vba_code = f"""
                    Solid.Subtract "{CutElement}", "{CutterElement}"
        """
        self.prj.model3d.add_to_history(f"Cur Structure {CutElement}", vba_code)
        
        




#################################################################################
# Curves Class
#################################################################################

class Curves:
    """
    Curve class is given an Euler Bezier and Cosinis curves (S-Curves). It was created to use the mathematical expressions 
    to assist the above CST wire to solid functions.
    """
    
    def __init__(self, Span_X, Span_Y, NumberOfPoints):
        """
        

        Parameters
        ----------
        Span_X : int/float
            Set the span of the X (Lenght) of the S-Curve
        Span_Y : int/float
            Set the span of the Y (Lenght) of the S-Curve
        NumberOfPoints : int
            Set the number of points. This will be the number of points of your curve.

        Returns
        -------
        None.

        """
        
        self.Span_X = Span_X
        self.Span_Y = Span_Y
        self.NumberOfPoints = NumberOfPoints
        self.t = np.arange(0, 1, 1/self.NumberOfPoints)
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
        """
        The cosinus function is SpanY*(cos((pi/(2*SpanX))*t)^2)  ---->> t in the range of 0 to SpanY

        Returns
        -------
        curve : 2 dimentional np array
            curve[:,0] - X Param of the Bezier Curve
            curve[:,1] - Y Param of the Bezier Curve.

        """
        
        P = self.Span_Y
        L = self.Span_X
        stepSize = L/len(self.t)
        t = np.arange(0,L, stepSize) 
        Func = P*(np.cos((np.pi/(2*L))*t)**2)  
        curve = np.vstack((t, Func[::-1])).T
        
        return curve
    
    
    
    
    def Euler_Curve(self):
        """
        Construct an Euler S-bend from (0, 0) to (Span_X, Span_Y)
        Using Fresnel integrals and linear curvature.

        Returns
        -------
        curve : 2 dimentional np array
            curve[:,0] - X Param of the Bezier Curve
            curve[:,1] - Y Param of the Bezier Curve.


        """
        """ 
        Returns a smooth, continuous-curvature Euler S-bend from (0, 0) to (Span_X, Span_Y)
        using Fresnel integrals and linear curvature.
        """
        Span_X = self.Span_X
        Span_Y = self.Span_Y
        num_pts = self.NumberOfPoints
        
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



#################################################################################
# Post Process Class
#################################################################################


class PostProcess:
    """
    This class is used to export and also plot data that is generated with CST
    """
    
    def __init__(self, CST_obj: CST_Commands):
            """
            Initialize PostProcess with a CST_Commands instance.

            Parameters
            ----------
            CST_obj : CST_Commands
                Instance of CST_Commands containing self.prj
            """
            # Store only the CST project, not the full CST object
            self.prj = CST_obj.prj

    
    def Export_SParameters(self, Parameters):
        """
        

        Parameters
        ----------
        Parameters["Name"] : str
            Set the name of the S_parameter file. For example Parameters["Name"] = "Export_S_Params".
        Parameters["Path"] : str
            Set the path to where to save the file. For example Parameters["Path"] = "C:/../../CST_Files".
        Parameters["Impedance"] : int/float
            Set the impedance for the S-Parameters. Per defoult is set to 50 Ohms.
        Parameters["Format"] : str
            Set the format that the S-Parameters will be exportet. You can choose between: "MA", "DB" and "RI".
       

        Raises
        ------
        ValueError
            Error choosing the wrong format.

        Returns
        -------
        None.

        """

        name = Parameters["Name"]
        path = Parameters["Path"]
        impedance = Parameters["Impedance"]
        format = Parameters["Format"]


        # Add patha nd name
        path_and_name = path + "/" + name
        
        # Check the format 
        form_list = ["MA", "RI", "DB"]
        if format not in form_list:
            raise ValueError("Format is invalid. You can choose between 'MA', 'RI' or 'DB'")
        else:
            vba_code = f"""
                            With TOUCHSTONE
                                .Reset
                                .FileName "{path_and_name}"
                                .Impedance {impedance}
                                .ExportType "S"
                                .Format "{format}"
                                .FrequencyRange "Full"
                                .Renormalize True
                                .UseARResults False
                                .SetNSamples 100
                                .Write
                            End With
            """
            self.prj.model3d.add_to_history(f"Export S-Parameters", vba_code)



    def read_snp(self, Parameters):
        """
        Reads an S-parameter Touchstone file (.s2p, .s3p, .s4p, etc.) and returns
        frequency array and S-matrix array.

        Parameters
        ----------
        Parameters["File Name"] : str
            Set the name of the S_parameter file. For example an .s2p file Parameters["File Name"] = "Export_S_Params.s2p".
        Parameters["Path"] : str
            Set the path to where to save the file. For example Parameters["Path"] = "C:/../../CST_Files".

        Returns
        -------
        freq : np.ndarray
            Frequency points (Hz or as in file).
        S_params : np.ndarray
            Array of shape (num_points, N, N) containing complex S-matrices.
        """
        
        name = Parameters["File Name"]
        path = Parameters["Path"]
        freq = []
        S_params = []
        
        # Add patha nd name
        filename = path + "/" + name

        with open(filename, 'r') as f:
            lines = f.readlines()

        # Remove comments and empty lines
        data_lines = [line for line in lines if not line.startswith('!') and line.strip() != '']

        # Detect number of ports from the first line starting with '#' (the format line)
        header_line = next(line for line in data_lines if line.startswith('#'))
        tokens = header_line.lower().split()
        # Format line example: "# GHZ S MA R 50"
        # tokens: ['#', 'ghz', 's', 'ma', 'r', '50']
        # Number of ports = sqrt((len(data_values)/2))
        # We'll deduce ports later from data

        # Remove header line from data
        data_lines = [line for line in data_lines if not line.startswith('#')]

        # Flatten all numbers
        all_numbers = []
        for line in data_lines:
            all_numbers.extend(list(map(float, line.strip().split())))

        # Detect number of ports automatically
        # Each frequency point has 1 frequency + 2 * N^2 numbers
        total_numbers = len(all_numbers)
        freq_guess = all_numbers[0]  # first number is frequency
        # Solve N from total numbers: total_numbers / num_points = 1 + 2*N^2
        # We first estimate number of points: consecutive freq differences are >0
        # Simplest approach: assume last frequency = last block, or just compute sqrt
        # Let's make it simple: user provides N or infer from length
        # We'll infer N assuming uniform number of points:
        # Count number of frequencies by detecting differences > 0 (simple)
        freqs = []
        i = 0
        while i < len(all_numbers):
            freqs.append(all_numbers[i])
            i += 1
            # Cannot infer N yet from single pass; we assume N=2 or N=4 user can change
            break
        # Instead, require user to provide ports or infer from filename
        # We'll infer from filename extension: s2p -> N=2, s4p -> N=4, etc.
        import re
        match = re.search(r's(\d+)p', filename.lower())
        if match:
            N = int(match.group(1))
        else:
            raise ValueError("Cannot detect number of ports from filename. Use s2p, s3p, s4p format.")

        # Compute number of points
        block_size = 1 + 2 * N**2
        num_points = len(all_numbers) // block_size

        for i in range(num_points):
            start = i * block_size
            block = all_numbers[start:start + block_size]

            freq.append(block[0])

            # Build NxN S-matrix
            S_matrix = np.zeros((N, N), dtype=complex)
            for r in range(N):
                for c in range(N):
                    mag = block[1 + 2*(r*N + c)]
                    angle = block[2 + 2*(r*N + c)]
                    S_matrix[r, c] = mag * np.exp(1j * np.deg2rad(angle))
            S_params.append(S_matrix)

        return np.array(freq), np.array(S_params)
        
    


            


        