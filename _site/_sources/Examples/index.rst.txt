Examples
========

Python code example how to start and call the CST from Python. 
To call the following in your Python IDLE.

.. code-block:: python

	import numpy as np 
	import matplotlib.pyplot as plt
	import pandas as pandas
	import sys
	import os
	# Add the directory containing the project to sys.path
	current_path = os.path.dirname(os.path.abspath('__file__'))
	sys.path.append(current_path)
	from CST_Constructor import CST_Commands, Curves


	# #Open new CST Project 
	# Call CST script 
	obj = CST_Commands()
	
	# #List of all the available CSt Projects that can be called from this scrip
	# obj.New_Project("CS")
	# obj.New_Project("DS")
	# obj.New_Project("EMS")
	# obj.New_Project("FD3D")
	# obj.New_Project("MPS")
	obj.New_Project("MWS")
	# obj.New_Project("PCBS")
	# obj.New_Project("PS")

	
	# #Open existing CST Project
	# obj.Open_Project("C:/...../Test_Save_Code.cst")
	
	# #Open existing CST Project
	# obj.Open_Project("C:/...../Test_Save_Code.cst")
	
	
	# #Save CST Project
	# obj.Save_Project("C:/...../CST", "Test_Save_Code", False)



Full Code Examples:
-------------------

All detailed examples can be found below:

.. toctree::
   :maxdepth: 1
   :caption: Available Examples

   Example_GSG_Pads_with_Bondwires
   Example_GSGSG_Pads_with_Bondwires
   Example_MZM
   Example_Phase_Modulator
   Example_Waveguide_FDE
