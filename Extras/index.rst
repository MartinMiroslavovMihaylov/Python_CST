Curves ClassExamples
========

The class constructor gives the user the possibility to create S-bends using different mathematical functions.
The mathematical functions are implemented in the Curves class and can be called using the following code:


.. code-block:: python
    
	from CST_Constructor import Curves
	
	
	# Define Curves Parameters and Data
	Lenght = 100
	Offset = 40
	points = 100

		
	# Generate the Bezier and Cos points
	ObjCurves = Curves(Lenght, Offset, points)
	BezierCuve = ObjCurves.Bezier_Curve()
	CosinusCurve = ObjCurves.Cosinus_Curve()
	EulerCurve = ObjCurves.Euler_Curve()


	plt.figure()
	plt.plot(BezierCuve[:,0], BezierCuve[:,1], color = "red", label = "Bezier Curve")
	plt.plot(CosinusCurve[:,0], CosinusCurve[:,1], color = "blue", label = "Cosinus Curve")
	plt.plot(EulerCurve[:,0], EulerCurve[:,1], color = "green", label = "Euler Curve")
	plt.xlabel("Length S-Bend/$\mu m$")
	plt.ylabel("Offset S-Bends/ $\mu m$")
	plt.legend(loc = "best")
	plt.grid()
	plt.show()


