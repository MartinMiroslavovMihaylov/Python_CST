# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 15:49:14 2024

@author: marti
"""


import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.integrate as integrate
from scipy.optimize import fsolve
plt.rcParams.update({"font.size":22})


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
        
        def IntegralY(A,R):
            return np.sin( (sigma**2/2) + (A*Sigma/R) )
            
        def IntegralX(A,R):
            return np.sin( (sigma**2/2) + (A*Sigma/R) )
            
            
        L = 50

        A = np.sqrt(L/( (1/Rmin) - (1/Rmax)))
        Ylen = self.Span_Y
        R_eff = 4
        sigma = 90
        x = []
        y = []
        for i in range(len(self.t)):
            x.append(quad(A*IntegralX(sigma), 0, L/A))
            y.append(quad(A*IntegralY(SIGMA), 0, L/A))
            
            
        x = np.array(x,dtype=float)
        y = np.array(y,dtype=float)  
        curve =  np.vstack((x, y)).T
    
    
    

# obj = Curves(50,20,100)
# curve = obj.Bezier_Curve()
# curve2 = obj.Cosinus_Curve()


# plt.figure()
# plt.plot(curve[:,0], curve[:,1], label = "Bezier")
# plt.plot(curve2[:,0], curve2[:,1], label = "Cosinus")
# plt.xlabel("L/ $\mu m$")
# plt.ylabel("P/ $\mu m$")
# plt.legend(loc = "best")
# plt.grid()
# plt.show()


# import matplotlib.pyplot as plt

# def euler_method(func, y0, t0, tn, h):
#     """
#     Approximate the solution of the ODE using Euler's method.

#     Parameters:
#         func: The differential equation function, f(t, y).
#         y0: Initial value of the dependent variable.
#         t0: Initial value of the independent variable.
#         tn: Final value of the independent variable.
#         h: Step size.

#     Returns:
#         Two lists containing the values of t and y at each step.
#     """
#     t_values = [t0]
#     y_values = [y0]

#     while t_values[-1] < tn:
#         t = t_values[-1]
#         y = y_values[-1]
#         y_next = y + h * func(t, y)
#         t_values.append(t + h)
#         y_values.append(y_next)

#     return t_values, y_values

# # Example: Solve the ODE dy/dt = -2 * y with y(0) = 1
# def differential_equation(t, y):
#     return -2 * y

# t0, tn = 0, 4
# y0 = 20
# h = 0.1

# t_values, y_values = euler_method(differential_equation, y0, t0, tn, h)

# # Plot the solution
# plt.plot(t_values, y_values, label='Euler Method')
# plt.title("Euler's Method for dy/dt = -2y")
# plt.xlabel('t')
# plt.ylabel('y(t)')
# plt.legend()
# plt.show()








# import scipy.integrate as integrate
# from scipy.optimize import fsolve




# R_eff = 6  # effective radius of the bend
# A = 1.3  # clothoid parameter


# L_max = 0  # starting point of L_max
# precision = 0.05  # increasement of L_max at each iteration
# tolerance = 0.01  # difference tolerance of the derivatives

# # determine L_max
# while True:
#     L_max = L_max + precision  # update L_max
#     Ls = np.linspace(0, L_max, 50)  # L at (x1,y1)
#     x1 = np.zeros(len(Ls))  # x coordinate of the clothoid curve
#     y1 = np.zeros(len(Ls))  # y coordinate of the clothoid curve

#     # compute x1 and y1 using the above integral equations
#     for i, L in enumerate(Ls):
#         y1[i], err = integrate.quad(lambda theta: A * np.sin(theta**2 / 2), 0, L / A)
#         x1[i], err = integrate.quad(lambda theta: A * np.cos(theta**2 / 2), 0, L / A)

#     # compute the derivative at L_max
#     k = -(x1[-1] - x1[-2]) / (y1[-1] - y1[-2])
#     xp = x1[-1]
#     yp = y1[-1]
#     # check if the derivative is continuous at L_max
#     R = np.sqrt(
#         ((R_eff + k * xp - yp) / (k + 1) - xp) ** 2
#         + (-(R_eff + k * xp - yp) / (k + 1) + R_eff - yp) ** 2
#     )
#     if np.abs(R - A**2 / L_max) < tolerance:
#         break

# # after L_max is determined, R_min is also determined
# R_min = A**2 / L_max

# # getting the coordinates of the second clothoid curve by mirroring the first curve with respect to y=-x+R_eff
# x3 = np.flipud(R_eff - y1)
# y3 = np.flipud(R_eff - x1)


# # solve for the parameters of the circular curve
# def circle(var):
#     a = var[0]
#     b = var[1]
#     Func = np.empty((2))
#     Func[0] = (xp - a) ** 2 + (yp - b) ** 2 - R_min**2
#     Func[1] = (R_eff - yp - a) ** 2 + (R_eff - xp - b) ** 2 - R_min**2
#     return Func


# a, b = fsolve(circle, (0, R_eff))

# # calculate the coordinates of the circular curve
# x2 = np.linspace(xp + 0.01, R_eff - yp - 0.01, 50)
# y2 = -np.sqrt(R_min**2 - (x2 - a) ** 2) + b


# # obtain the coordinates of the whole Euler bend by concatenating three pieces together
# x_euler = np.concatenate((x1, x2, x3))
# y_euler = np.concatenate((y1, y2, y3))

# # the conventional circular bend is simply given by a circle
# x_circle = np.linspace(0, R_eff, 100)
# y_circle = -np.sqrt(R_eff**2 - (x_circle) ** 2) + R_eff

# # plotting the shapes of the Euler bend and the circular bend
# plt.plot(x_euler, y_euler, label="Euler bend")
# plt.plot(x_circle, y_circle, "--", label="Circular bend")
# plt.axis("equal")
# plt.ylim(-1, R_eff + 1)
# plt.legend()
# plt.show()