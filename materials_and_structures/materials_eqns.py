from sympy import *


a_crit, a_i, A, Q, delta_sigma, m, N_f = symbols(
    "a_crit, a_i, A, Q, delta_sigma, m, N_f")

'''
A               Paris law constant
m               Slope of log log linear stage two of the Paris law graph
delta_sigma     Stress range (ignoring compressive stresses)
Q               Shape factor
a_crit          Critical crack length for fast failure via plastic collapse
a_i             Inital crack length
N_f             Number of cycles to failure
'''

crack_lifetime_eqn = a_crit ** (1 - m / 2) / (1 - m / 2) - a_i ** (
    1 - m / 2) / (1 - m / 2) - A * (Q * delta_sigma * pi ** 0.5) ** m * N_f
