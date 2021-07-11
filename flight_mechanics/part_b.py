"""Part B solvable equations as per the derivations sheet"""

import numpy as np
from sympy import *


'''
1.1 - Body Axes Coordinate System

X - force forwards
Y - starboard force (right wing facing sideways)
Z - force downwards
'''


'''
1.2 - Inertial Rates of Change
'''


def Force_on_aircraft_in_body_reference_frame(m, V_B, V_dot_B, omega_B):
    """(14) Force (1x3) on aircraft in body reference frame
    
    V_B, V_dot_B, omega_B, returns - (1x3) matrices"""
    return m * (V_dot_B + omega_B.cross(V_B))


def Moment_on_aircraft_in_body_reference_frame(h_B, h_dot_B, omega_B):
    """(15) Moment (1x3) on aircraft in body reference frame
    
    h_B, h_dot_B, omega_B, returns - (1x3) matrices"""
    return h_dot_B + omega_B.cross(h_B)


'''
2.1 - Mass Moments of Inertia
'''


def v_from_omega_r(w, r):
    """Velocity (1x3) vector from  angular velocity (1x3) omega and (1x3) r"""
    return w.cross(r)

def h_angular_momentum(r, v):
    """(16) Angular momentum h (1x3)"""
    m = Symbol("m")
    return integrate(r.cross(v), m)

#m = symbols("m")
#h_angular_momentum(Matrix([12,2,3]), Matrix([100,2,3]))
'''
h, I, w = symbols("h, I, w")
expr = solve(h - I * w, h)[0]
result = expr.subs(I, Matrix([[1, 2], [3, 4]]))
result = result.subs(w, Matrix([[1, 2], [3, 4]]))
print(result)
'''

def h_from_I_omega(I, omega):
    """(30) Find h from I*oemga"""
    I, omega = symbols("I, omega")
    return I * omega


def omega_from_h_I(h, I):
    """(30) But inversed"""
    h, I = symbols("h, I")
    return I**-1 * h


'''
2.2 - Force Equations
'''


def X_force(omega_B, V_B, m):
    """(34) returns X body force, make sure V_B differentiates to correct value, use V_B equation not value!"""
    t = Symbol("t")
    return m * (diff(V_B[0], t) + omega_B[1] * V_B[2] - omega_B[2] * V_B[1])


# Example X_force
#t = Symbol("t")
#print(X_force(Matrix([10, 1, 3]), Matrix([10 * t ** 2, 20 * t, 3]), 2))


def Y_force(omega_B, V_B, m):
    """(35) returns Y body force, make sure V_B differentiates to correct value, use V_B equation not value!"""
    t = Symbol("t")
    return m * (diff(V_B[1], t) + omega_B[2] * V_B[0] - omega_B[0] * V_B[2])


def Z_force(omega_B, V_B, m):
    """(36) returns Z body force, make sure V_B differentiates to correct value, use V_B equation not value!"""
    t = Symbol("t")
    return m * (diff(V_B[2], t) + omega_B[0] * V_B[1] - omega_B[1] * V_B[0])


'''
2.3 - Moment Equations

I the mass moment of intertia does not vary with time, only H and omega
'''


'''
2.4 - Symmetric Aircraft Approximation

For an aircraft with a plane of symetry assume I_x_y = I_y_z
'''


'''
2.5 - Linearisation of Forces

Assume perturbations and angular veolcities are small, products of 2 small numbers are assumed negligable during linearisation
'''


def delta_X_force(find, **kwargs):
    """(55) Change in X direction force due to perturbations"""
    delta_X, m, u_dot = symbols("delta_X, m, u_dot")
    expr = solve(m * u_dot - delta_X, find)[0]
    return expr.subs(kwargs.items())


def delta_Y_force(find, **kwargs):
    """(58) Change in Y direction force due to perturbations"""
    delta_Y, m, v_dot, r, U_infinity = symbols("delta_Y, m, v_dot, r, U_infinity")
    expr = solve(m * (v_dot + r * U_infinity) - delta_Y, find)[0]
    return expr.subs(kwargs.items())


def delta_Z_force(find, **kwargs):
    """(61) Change in Z direction force due to perturbations"""
    delta_Z, m, w_dot, q, U_infinity = symbols("delta_z, m, w_dot, q, U_infinity")
    expr = solve(m * (w_dot - q * U_infinity) - delta_Z, find)[0]
    return expr.subs(kwargs.items())


'''
2.6 - Linearisation of Moments
'''


'''
3.1 - Gravity Force

For an aircraft at a climb angle of gamma_0 of equilibrium, gravity force can be calculated for X and Z directions
'''


def X_g(find, **kwargs):
    """(70) forwards gravity component, angles in radians"""
    X_g, m, g, gamma_0, theta = symbols("X_g, m, g, gamma_0, theta")
    expr = solve(- m * g * sin(gamma_0 + theta) - X_g, find)[0]
    return expr.subs(kwargs.items())


# Example X_g
#print(X_g(find="theta", X_g = 10, m = 10, g=9.81, gamma_0=0.1))


def Z_g(find, **kwargs):
    """(71) downwards gravity component, angles in radians"""
    Z_g, m, g, gamma_0, theta = symbols("Z_g, m, g, gamma_0, theta")
    expr = solve(m * g * cos(gamma_0 + theta) - Z_g, find)[0]
    return expr.subs(kwargs.items())


def delta_X_g(find, **kwargs):
    """(75) change in forwards gravity component, angles in radians"""
    delta_X_g, m, g, gamma_0, theta = symbols("delta_X_g, m, g, gamma_0, theta")
    expr = solve(- m * g * cos(gamma_0 * theta) - delta_X_g)
    return expr.subs(kwargs.items())


def delta_Z_g(find, **kwargs):
    """(80) change in downards gravity component, angles in radians"""
    delta_Z_g, m, g, gamma_0, theta = symbols("delta_Z_g, m, g, gamma_0, theta")
    expr = solve(- m * g * sin(gamma_0 * theta) - delta_Z_g)
    return expr.subs(kwargs.items())


'''
3.2 - Aerodynamic Derivatives

It is assumed that aerodynamic forces and moments are functions of angular and linear disturbance velocities and their derivatives

Circle above quantity indicates its a dimensional aerodynamic derivative, a subscript of 0 represents the equilibrium flight condition
'''