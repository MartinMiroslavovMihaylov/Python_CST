
Install
=======
These Python libraries can be used with Python version 3.11.0. Please check the CST Help menu to verify which Python versions are compatible with your CST installation.

For the Moment you will need the following Python Libraries. ::

   pip install numpy
   pip install pandas
   pip install matplotlib

   
CST Installation and file location 
==============================================================

The CST Studio Suite installation comes with a working Python interpreter, which requires no further setup to start using it with the CST Python Libraries. Various packages like numpy and scipy are pre-installed. You can find it under:

	Windows: “<CST_STUDIO_SUITE_FOLDER>\AMD64\python\python”
	Linux: “<CST_STUDIO_SUITE_FOLDER>/LinuxAMD64/python/python”

where <CST_STUDIO_SUITE_FOLDER> should be replaced with the path to the CST Studio Suite installation on your system.

Custom Python interpreter
To make sure your interpreter is able to load the CST Python Libraries, a minimal package can be installed that links to the actual CST Studio Suite installation:

pip install --no-index --find-links "<CST_STUDIO_SUITE_FOLDER>/Library/Python/repo/simple" cst-studio-suite-link
where <CST_STUDIO_SUITE_FOLDER> should be replaced with the path to the CST Studio Suite installation on your system. Do not install this package if you don’t have CST Studio Suite installed on the same machine as the python interpreter.

Testing your environment
You have successfully set up your Python environment if you are able to execute the following code without an error:

import cst
print(cst.__file__)  # should print '<PATH_TO_CST_AMD64>\python_cst_libraries\cst\__init__.py'

