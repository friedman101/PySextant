import sympy as sp
from math import *

r_sun1x = sp.symbols('r_sun1x')
r_sun1y = sp.symbols('r_sun1y')
r_sun1z = sp.symbols('r_sun1z')
r_sun2x = sp.symbols('r_sun2x')
r_sun2y = sp.symbols('r_sun2y')
r_sun2z = sp.symbols('r_sun2z')
r_mex = sp.symbols('r_mex')
r_mey = sp.symbols('r_mey')
r_mez = sp.symbols('r_mez')
r_sun1 = sp.Matrix([[r_sun1x],[r_sun1y],[r_sun1z]])
r_sun2 = sp.Matrix([[r_sun2x],[r_sun2y],[r_sun2z]])
r_me = sp.Matrix([[r_mex], [r_mey], [r_mez]])
eq1 = sp.Eq( r_me.dot(r_me), 1 )
# pi/2 - acos(dot(r_sun, r_me)) = theta_meas
# theta1 = cos(pi/2 - theta_meas)
theta1 = sp.symbols('theta1')
theta2 = sp.symbols('theta2')
eq2 = sp.Eq( r_sun1.dot(r_me), theta1 )
eq3 = sp.Eq( r_sun2.dot(r_me), theta2 )

sp.solve([eq1, eq2, eq3], [r_mex, r_mey, r_mez])