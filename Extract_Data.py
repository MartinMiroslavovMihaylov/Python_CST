import numpy as np 
import matplotlib.pyplot as plt
import os 
import sys
# Add the directory containing the project to sys.path
current_path = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_path)
from VBA_Test_Syntax import *
from Components import *


def Port_Data(Parameters):

    PortNumber = Parameters["Port Number"]
    ModeNumber = Parameters["Mode Number"]

    data = 'Sub Main () ' \
        '\n.GetFcutoff ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetFrequency ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetModeType ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetBeta ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetAlpha ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetWaveImpedance ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetLineImpedance ' + '"' + str(PortNumber) + '"' + ',' + '"' +str(ModeNumber) + '"' + \
        '\n.GetNumberOfModes ' + '"' + str(PortNumber) + \
    '\nEnd Sub'
    data = ''.join(data)

    return Port
    

