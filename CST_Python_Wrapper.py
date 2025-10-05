import numpy as np 
import matplotlib.pyplot as plt
import os 
import sys
import cst
import cst.interface
import json


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
        list_of_materials = ["Si", "LiNbO3", "SiO2", "Au", "Gold"]

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
        else:
            raise ValueError(
                            """
                            For now only the following materials can be added with this python script:
                            "Si"
                            "LiNbO3"
                            "SiO2"
                            "Au or Gold"
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
                        .Epsilon "2.0851"
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
        self.prj.model3d.add_to_history("create brick", vba_code)


        

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
            

    def Poligon_3D(WGName, Points):
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

        WGName  = WGName
        lines = []

        for i in range(1, len(Points['X'])):
            lines.append(f'.Point "{Points["X"][i]}", "{Points["Y"][i]}"')

        # join them with newlines
        line_block = "\n".join(lines)

        vba_code = f"""
                    With Polygon3D
                    .Reset
                    .Name "{WGName}"
                    .Curve "{WGName}"
                    {line_block}
                    .Create
                    End With
                    """
        self.prj.model3d.add_to_history("3d polygon", vba_code) 

    
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
                        .Folder "' + NameFolder + '"
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
            self.prj.model3d.add_to_history("create bondwire", vba_code)
            

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
            self.prj.model3d.add_to_history("create bondwire", vba_code)




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
            vba_code3 = f"""
                    .Zrange "{Z_min}", "{Z_max}" 
                    .Xcenter "{Xcenter}" 
                    .Ycenter "{Ycenter}" 
                    .Segments "0" 
                    .Create 
                End With
                """
            full_vba_code = vba_code + vba_code3
        
        else:
            raise ValueError("Orentation Axis can be only X, Y and Z!")
        
        self.prj.model3d.add_to_history("create cylinder", full_vba_code)




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
                    .SolidName "component1:{SolidName}"
                    .Name "{CurveName}"
                    .Folder "{NameFolder}"
                    .Material "{Material}"
                    .KeepWire "False"
                    .ConvertToSolidShape
                    End With
                    """
        self.prj.model3d.add_to_history("to solid", vba_code)







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
        potential = Parameters["Picked Port Polarity"]
        marked_face = Parameters["Picked Component Name"]
        Number_of_picks = Parameters["Number of picks"]
        lines = []

        lines = []
        if Number_of_picks > 2:
            for i in range(len(potential)):
                lines.append(f'.AddPotentialPicked "{set_number}", "{potential[i]}", "{marked_face[i]}", "1"')

            # join them with newlines
            line_block = "\n".join(lines)

            vba_code = f"""
                        With Port
                        .Reset
                        .PortNumber "{PortNumber}"
                        .NumberOfModes "5"
                        .Orientation "{Orientation}"
                        .Coordinates "{Coordinates}"
                        .PortOnBound "False"
                        .ClipPickedPortToBound "False"
                        .XrangeAdd "{Span[0]}", "{Span[0]}"
                        .YrangeAdd "{Span[1]}", "{Span[1]}"
                        .ZrangeAdd "{Span[2]}", "{Span[2]}"
                        .AdjustPolarization "True"
                        .PolarizationAngle "0"
                        .SingleEnded "False"
                        {line_block}
                        .Create
                        End With
                        """
            self.prj.model3d.add_to_history("create waveguide port", vba_code)

        else:
            vba_code = f"""
                        With Port
                        .Reset
                        .PortNumber "{PortNumber}"
                        .NumberOfModes "5"
                        .ReferencePlaneDistance "0"
                        .Coordinates "{Coordinates}"
                        .Orientation "{Orientation} "
                        .PortOnBound "False"
                        .ClipPickedPortToBound "False"
                        .XrangeAdd "{Span[0]}", "{Span[0]}"
                        .YrangeAdd "{Span[1]}", "{Span[1]}"
                        .ZrangeAdd "{Span[2]}", "{Span[2]}"
                        .AdjustPolarization "True"
                        .PolarizationAngle "0"
                        .SingleEnded "False"
                        .AddPotentialPicked "{set_number}", "{potential}", "{marked_face}", "1"
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
        Accuracy = Parameters["Accuracy"]
        ModesOnly = Parameters["Caclculate Modes Only"]
        AutoImpedance = Parameters["Auto Impedance"]
        Impedance = Parameters["Impedance"]
        Source = Parameters["Source Port"]

        vba_code = f"""
            Mesh.SetCreator "High Frequency"
            With FDSolver
            .Reset
            .SetMethod "Tetrahedral", "General purpose"
            .OrderTet "Second"
            .OrderSrf "First"
            .Stimulation "All", "All" 
            .ResetExcitationList 
            .AutoNormImpedance "False"
            .NormingImpedance "50"
            .ModesOnly "False"
            .ConsiderPortLossesTet "True"
            .SetShieldAllPorts "False"
            .AccuracyHex "1e-6"
            .AccuracyTet "1e-5"
            .AccuracySrf "1e-3"
            .LimitIterations "False"
            .MaxIterations "0"
            .SetCalcBlockExcitationsInParallel "True", "True", ""
            .StoreAllResults "False"
            .StoreResultsInCache "False"
            .UseHelmholtzEquation "True"
            .LowFrequencyStabilization "True"
            .Type "Direct"
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
            .SetNumberOfResultDataSamples "5001"
            .SetResultDataSamplingMode "Automatic"
            .SweepWeightEvanescent "1.0"
            .AccuracyROM "1e-4"
            .AddSampleInterval "", "", "5", "Automatic", "True"
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
            Name = MonitorType + ' (wl=' + str(Wavelength) + ')'
            vba_code = f"""
                With Monitor
                .Reset
                .Name "{Name}"
                .Domain "Wavelength"
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
            self.prj.model3d.add_to_history("create Efield monitor", vba_code)
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

        Name = Parameters["Name Object"]
        PosX = Parameters["Position X"]
        PosY = Parameters["Position Y"]
        PosZ = Parameters["Position Z"]
        Shape = Parameters["Structure Type"]

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
                    .Transform "{Shape}", "Translate"
                    End With
                    """
        self.prj.model3d.add_to_history("translation", vba_code)


