"""Part A solvable equations as per the derivations sheet"""

import numpy as np
from sympy.solvers import solve
from sympy import symbols
from sympy import *


'''
1.1 - Coefficients
'''


def C_L_star(find, **kwargs):
    """(1) Total aircraft lift coefficient"""
    C_L_star, L_star, rho, V, S = symbols("C_L_star, L_star, rho, V, S")
    expr = solve(L_star / (0.5 * rho * V**2 * S) - C_L_star, find)[0]
    return expr.subs(kwargs.items())


# Find V
V = C_L_star(find="V", C_L_star=2, L_star=100, rho=1.2, S=15)
pprint(C_L_star(find="V")) # Returns negative value sometimes as + or - for a square root
print(V)
print(V.evalf())


def C_L(find, **kwargs):
    """(1) Main wing lift coefficient"""
    C_L, L, rho, V, S = symbols("C_L, L, rho, V, S")
    expr = solve(L / (0.5 * rho * V**2 * S) - C_L, find)[0]
    return expr.subs(kwargs.items())


def C_W(find, **kwargs):
    """(1) Weight coefficient"""
    C_W, rho, V, S, W = symbols("C_W, rho, V, S, W")
    expr = solve(W / (0.5 * rho * V**2 * S) - C_W, find)[0]
    return expr.subs(kwargs.items())


def C_L_T(find, **kwargs):
    """(1) Tailplane lift coefficient"""
    C_L_T, L_T, rho, V, S_T = symbols("C_L_T, L_T, rho, V, S_T")
    expr = solve(L_T / (0.5 * rho * V**2 * S_T) - C_L_T, find)[0]
    return expr.subs(kwargs.items())


def C_M_0(find, **kwargs):
    """(1) Moment coefficient"""
    C_M_0, M, rho, V, S, c = symbols("C_M_0, M, rho, V, S, c")
    expr = solve(M / (0.5 * rho * V**2 * S * c) - C_M_0, find)[0]
    return expr.subs(kwargs.items())


'''
1.3 - Moments
'''


def K(find, **kwargs):
    """(11) Tail volume fraction

    l - is distance between aerodynamic centres of wing and tailplane"""
    K, S_T, l, S, c = symbols("K, S_T, l, S, c")
    expr = solve((S_T * l) / (S * c) - K, find)[0]
    return expr.subs(kwargs.items())


'''
1.4 - Forces
'''


def C_L_T_from_lift_coefficients(find, **kwargs):
    """(15) Calculate tail lift coefficient from coefficients and angles

    a_1 - C_L_T partial derivative wrt alpha_T_eff (tail plane aoa)
    a_2 - C_L_T partial derivative wrt eta (elevator angle)
    a_3 - C_L_T partial derivative wrt beta (elevator trim tab angle)
    """
    C_L_T, a_1, a_2, a_3, alpha_T_eff, eta, beta = symbols("C_L_T, a_1, a_2, a_3, alpha_T_eff, eta, beta")
    expr = solve(a_1 * alpha_T_eff + a_2 * eta + a_3 * beta - C_L_T, find)[0]
    return expr.subs(kwargs.items())


def C_M_H(find, **kwargs):
    """(16) Tailplane hinge moment coefficient

    b_1 - C_M_H partial derivative wrt alpha_T_eff (effective tailplane aoa)
    b_2 - C_M_H partial derivative wrt eta (elevator angle)
    b_3 - C_M_H partial derivative wrt beta (elevator trim tab angle)"""
    C_M_H, b_1, b_2, b_3, alpha_T_eff, eta, beta = symbols("C_M_H, b_1, b_2, b_3, alpha_T_eff, eta, beta")
    expr = solve(b_1 * alpha_T_eff + b_2 * eta + b_3 * beta - C_M_H, find)[0]
    return expr.subs(kwargs.items())


def alpha_T_eff(find, **kwargs):
    """(17) Effective tailplane angle of attack

    eta_alpha - derivative of eta with respect to alpha
    alpha_0 - wing zero lift angle
    alpha_s - tail setting angle
    e_T - tail span efficiency factor
    """
    alpha_T_eff, eta, eta_alpha, alpha, alpha_0, alpha_S, C_L_T, A_T, e_T = symbols("alpha_T_eff, eta, eta_alpha, alpha, alpha_0, alpha_S, C_L_T, A_T, e_T")
    expr = solve((1 - eta_alpha) * alpha + eta_alpha * alpha_0 + alpha_S - C_L_T / (np.pi * A_T * e_T) - alpha_T_eff, find)[0]
    return expr.subs(kwargs.items())


'''
1.5 - Stick Fixed Condition

Trim tab angle beta is zero, pilot applies a force on stick to overcome elevator hinge moment
Use (16) and (17) with beta = 0
'''


def eta_stick_fixed(find, **kwargs):
    """(20) Eta the elevator position in stick fixed mode"""
    eta_stick_fixed, C_L_T, a_1, a_2, alpha_T_eff = symbols("eta_stick_fixed, C_L_T, a_1, a_2, alpha_T_eff")
    expr = solve((C_L_T - a_1 * alpha_T_eff) / a_2 - eta_stick_fixed, find)[0]
    return expr.subs(kwargs.items())


'''
1.6 - Stick Free Condition

Elevator hinge moment overcome by trim tab, the elevator is then free to float
'''


def eta_stick_free(find, **kwargs):
    """(22) Eta the elevator position in stick free mode"""
    eta_stick_free, b_1, b_2, b_3, alpha_T_eff, beta = symbols("eta_stick_free, b_1, b_2, b_3, alpha_T_eff, beta")
    expr = solve(((-b_1 * alpha_T_eff) - (b_3 * beta)) / b_2 - eta_stick_free, find)[0]
    return expr.subs(kwargs.items())


def a_1_bar_eqn(find, **kwargs):
    """(24) a_1_bar defined constant"""
    a_1_bar, a_1, a_2, b_1, b_2 = symbols("a_1_bar, a_1, a_2, b_1, b_2")
    expr = solve(a_1 - a_2 * b_1 / b_2 - a_1_bar, find)[0]
    return expr.subs(kwargs.items())


def a_3_bar(find, **kwargs):
    """(24) a_3_bar defined constant"""
    a_3_bar, a_2, a_3, b_2, b_3 = symbols("a_3_bar, a_2, a_3, b_2, b_3")
    expr = solve(a_3 - a_2 * b_3 / b_2 - a_3_bar, find)[0]
    return expr.subs(kwargs.items())


def C_L_T_stick_free(find, **kwargs):
    """(25) C_L_T_stick_free"""
    C_L_T_stick_free = a_1_bar, a_3_bar, alpha_T_eff, beta = symbols("C_L_T_stick_free, a_1_bar, a_3_bar, alpha_T_eff, beta")
    expr = solve(a_1_bar * alpha_T_eff + a_3_bar * beta - C_L_T_stick_free, find)[0]
    return expr.subs(kwargs.items())


'''
2.1 - Trim drag

Trim drag is the component of lift created by the flight control surfaces
'''


def C_D_L(find, **kwargs):
    """(30) Main wing lift induced drag

    A - aspect ratio main wing b^2/S (wingspan squared over wing area)
    e - span efficiency factor main wing"""
    C_D_L, C_L, A, e = symbols("C_D_L, C_L, A, e")
    expr = solve(C_L**2 / (np.pi * A * e) - C_D_L, find)[0]
    return expr.subs(kwargs.items())


def C_D_L_star(find, **kwargs):
    """(30) Total lift induced drag, A and e appear to have same usage as C_D_L for only main wing

    A - aspect ratio main wing b^2/S (wingspan squared over wing area)
    e - span efficiency factor main wing"""
    C_D_L_star, C_L_star, A, e = symbols("C_D_L_star, C_L_star, A, e")
    expr = solve(C_L_star**2 / (np.pi * A * e) - C_D_L_star, find)[0]
    return expr.subs(kwargs.items())


def C_D_L_T(find, **kwargs):
    """(30) Tailplane lift induced drag

    A_T - aspect ratio of the tail
    e_T - spanwise efficiency factor of the tail"""
    C_D_L_T, C_L_T, A_T, e_T = symbols("C_D_L_T, C_L_T, A_T, e_T")
    expr = solve(C_L_T**2 / (np.pi * A_T * e_T) - C_D_L_T, find)[0]
    return expr.subs(kwargs.items())


def C_D_m(find, **kwargs):
    """(32) Trim drag coefficient from the 3 lift coefficient terms, areas, aspect ratios, and spanwise efficiency factors"""
    C_D_m, C_L, C_L_star, C_L_T, S, S_T, A, A_T, e, e_T = symbols("C_D_m, C_L, C_L_star, C_L_T, S, S_T, A, A_T, e, e_T")
    expr = solve((C_L**2 - C_L_star**2) / np.pi * A * e + C_L_T**2 * S_T / (np.pi * A_T * e_T * S) - C_D_m, find)[0]
    return expr.subs(kwargs.items())


def C_L_from_total_and_tailplane(find, **kwargs):
    """(35) Get main wing lfit coefficient from total lift and tailplane lift coefficients"""
    C_L, C_L_star, C_L_T, S, S_T = symbols("C_L, C_L_star, C_L_T, S, S_T")
    expr = solve(C_L_star - C_L_T * S_T / S - C_L, find)[0]
    return expr.subs(kwargs.items())


def sigma(find, **kwargs):
    """(43) Sigma constant for trim drag eqn (45) and others"""
    sigma, S, S_T, A, A_T, e, e_T = symbols("sigma, S, S_T, A, A_T, e, e_T")
    expr = solve(1 + (S * A * e) / (S_T * A_T * e_T) - sigma, find)[0]
    return expr.subs(kwargs.items())


def C_D_m_no_C_L(find, **kwargs):
    """(45) Trim drag coefficient when C-L of the main wing is unknown"""
    C_D_m, C_D_L_star, C_L_T, C_L_star, sigma, S, S_T = symbols("C_D_m, C_D_L_star, C_L_T, C_L_star, sigma, S, S_T")
    expr = solve(C_D_L_star * C_L_T * S_T * (sigma * C_L_T * S_T / (C_L_star * S) - 2) / (C_L_star * S) - C_D_m, find)[0]
    return expr.subs(kwargs.items())


def C_D_eqn_46(find, **kwargs):
    """(46) find C_D"""
    C_D, C_D_0, C_D_L_star, C_D_m = symbols("C_D, C_D_0, C_D_L_star, C_D_m")
    expr = solve(C_D_0 + C_D_L_star + C_D_m - C_D, find)[0]
    return expr.subs(kwargs.items())


def C_D_eqn_47(find, **kwargs):
    """(47) find C_D"""
    C_D, C_D_0, C_D_L, C_D_L_T = symbols("C_D, C_D_0, C_D_L, C_D_L_T")
    expr = solve(C_D_0 + C_D_L + C_D_L_T - C_D, find)[0]
    return expr.subs(kwargs.items())


'''
2.2 - Minimum Trim Drag

Minimum trim drag occurs when C_D_m is differentiated wrt C_L_T and set to zero
'''


def C_L_M_minimum(find, **kwargs):
    """(54) Calculate the minimal value of trim drag at a certain C_L_T, this does not mean the aircraft is trimmed at this value"""
    C_L_M_minimum, C_D_L_star, sigma = symbols("C_L_M_minimum, C_D_L_star, sigma")
    expr = solve((- C_D_L_star / sigma) - C_L_M_minimum, find)[0]
    return expr.subs(kwargs.items())


'''
3.1 -  Loiter/Cruise Drag at Minimum Trim Drag Condition
'''


def C_D_at_minimum_drag(find, **kwargs):
    """(59) Total drag at the minimum trim drag"""
    C_D_at_minimum_drag, C_D_0, C_D_L_star, sigma = symbols("C_D_at_minimum_drag, C_D_0, C_D_L_star, sigma")
    expr = solve((C_D_0 + C_D_L_star - C_D_L_star / sigma) - C_D_at_minimum_drag, find)[0]
    return expr.subs(kwargs.items())


'''
3.2 -  Maximising  Loiter/Cruise  Performance  at  Minimum Trim Drag Condition

(C_L**d / C_D) max is the optimal cruise performance, d is power dependning on required efficiency and propulsion
'''


def C_L_star_max_cruise_performance(find, **kwargs):
    """(72) C_L_star at the maximum cruise performance condition"""
    C_L_star_max_cruise_performance, A, C_D_0, d, e, sigma = symbols("C_L_star_max_cruise_performance, A, C_D_0, d, e, sigma")
    expr = solve(((C_D_0 * d)/((1 - 1 / sigma) * (2 - d)/ (np.pi * A * e)))**0.5 - C_L_star_max_cruise_performance, find)[0]
    return expr.subs(kwargs.items())


def C_D_max_cruise_performance(find, **kwargs):
    """(72) C_D at the maximum cruise performance condition"""
    C_D_max_cruise_performance, C_D_0, d = symbols("C_D_max_cruise_performance, C_D_0, d")
    expr = solve((C_D_0 * (1 + d / (2 - d))) - C_D_max_cruise_performance, find)[0]
    return expr.subs(kwargs.items())


'''
3.3 - Optimal lift-to-drag ratio for the maximum flight efficiency
'''


def C_L_star_over_C_D_opt(find, **kwargs):
    """(85) Optimal total lift-to-drag ratio for maximum fuel efficiency"""
    C_L_star_over_C_D_opt, A, A_T, C_D_0, e, e_T, d, S, S_T = symbols("C_L_star_over_C_D_opt, A, A_T, C_D_0, e, e_T, d, S, S_T")
    expr = solve((((S_T * np.pi * A_T * e_T / S + np.pi * A * e) * d * (2 - d) / (C_D_0 * 4))**0.5) - C_L_star_over_C_D_opt, find)[0]
    return expr.subs(kwargs.items())


'''
3.4 - Optimum CG Position

We need to ensure aircraft is trimmed at the minimum trim drag condition
'''


def h_opt(find, **kwargs):
    """(90) CG pos for trimmed flight at the minimum trim drag position"""
    h_opt, h_0, C_M_0, C_L_star, l, c, sigma = symbols("h_opt, h_0, C_M_0, C_L_star, l, c, sigma")
    expr = solve(h_0 - C_M_0 / C_L_star + l / (c * sigma) - h_opt, find)[0]
    return expr.subs(kwargs.items())


'''
4.1 - Longitudinal Static Stability
'''


def h_N(find, **kwargs):
    """(95) Nuetral point position from the leading edge"""
    h_N, h_0, C_L_T_alpha, C_L_star_alpha, K = symbols(
        "h_N, h_0, C_L_T_alpha, C_L_star_alpha, K")
    expr = solve((h_0 + K * C_L_T_alpha / C_L_star_alpha) - h_N, find)[0]
    return expr.subs(kwargs.items())


def h_S(find, **kwargs):
    """(97) Aircraft static margin

    C_L_T_alpha - derivative of C_L_T wrt alpha
    C_L_star_alpha - derivative of C_L_star wrt alpha"""
    h_S, h, h_0, C_L_T_alpha, C_L_star_alpha, K = symbols(
        "h_S, h, h_0, C_L_T_alpha, C_L_star_alpha, K")
    expr = solve((h_0 - h + K * C_L_T_alpha / C_L_star_alpha) - h_S, find)[0]
    return expr.subs(kwargs.items())


'''
4.2 - Tailplane Stability Derivative - Stick Fixed
'''


def k(find, **kwargs):
    """(102) little k is tailplane lift curve slope with tail effective angle of attack for the stick fixed condition"""
    k, a_1, A_T, e_T = symbols("k, a_1, A_T, e_T")
    expr = solve(a_1 * np.pi * A_T * e_T /
                 (np.pi * A_T * e_T + a_1) - k, find)[0]
    return expr.subs(kwargs.items())


def H_S_stack_fixed(find, **kwargs):
    """(105) Stability static margin stick fixed"""
    H_S_stack_fixed, C_L_star_alpha, h, h_0, k, K, eta_alpha = symbols(
        "H_S_stack_fixed, C_L_star_alpha, h, h_0, k, K, eta_alpha")
    expr = solve((h_0 - h + K * k * (1 - eta_alpha) / C_L_star_alpha), find)[0]
    return expr.subs(kwargs.items())


def C_L_star_alpha(find, **kwargs):
    """(108) C_L_star derivate wrt alpha, found using C_L_alpha and C_L_T_alpha"""
    C_L_star_alpha, C_L_alpha, C_L_T_alpha, S, S_T = symbols(
        "C_L_star_alpha, C_L_alpha, C_L_T_alpha, S, S_T")
    expr = solve((C_L_alpha + C_L_T_alpha * S_T / S) - C_L_star_alpha, find)[0]
    return expr.subs(kwargs.items())


'''
4.3 - Tailplane Stability Derivative - Stick Free
'''


def k_bar(find, **kwargs):
    """(114) tailplane lift curve slope with tail effective angle of attack for stick-free condition"""
    k_bar, a_1_bar, A_T, e_T = symbols("k_bar, a_1_bar, A_T, e_T")
    expr = solve(a_1_bar * np.pi * A_T * e_T /
                 (np.pi * A_T * e_T + a_1_bar) - k_bar, find)[0]
    return expr.subs(kwargs.items())


def H_S_stick_free(find, **kwargs):
    """(117) Stability static margin stick free"""
    H_S_stick_free, h, h_0, k_bar, K, C_L_star_alpha, eta_alpha = symbols(
        "H_S_stick_free,h,h_0,k_bar,K,C_L_star_alpha,eta_alpha")
    expr = solve((h_0 - h + K * (k_bar * (1 - eta_alpha) /
                                 C_L_star_alpha)) - H_S_stick_free, find)[0]
    return expr.subs(kwargs.items())


'''
5.1 - Newton’s Law for Steady Pull-up Manoeuvre
'''


def n(find, **kwargs):
    """(123) load factor for a steady pull-up manoeuvre"""
    n, V, R, g = symbols("n,V,R,g")
    expr = solve(V**2 / (g * R) + 1 - n, find)[0]
    return expr.subs(kwargs.items())


def R_found_from_steady_pullup(find, **kwargs):
    """(124) Radius of steady pullup manoeuvre"""
    R_found_from_steady_pullup, V, n, g = symbols(
        "R_found_from_steady_pullup,V,n,g")
    expr = solve(V**2 / (9.81 * (n - 1)) - R_found_from_steady_pullup, find)[0]
    return expr.subs(kwargs.items())


'''
5.2 -  Pitch Rate Induced Change in Tail Angle of Attack

The pitch rate q causes a change in the angle of attack of the tail delta_alpha_T_q
'''


def delta_alpha_T_q_from_V(find, **kwargs):
    """(128) Change in the angle of attack of the tail caused by pitch rate q"""
    delta_alpha_T_q_from_V, n, g, l, V = symbols(
        "delta_alpha_T_q_from_V,n,g,l,V")
    expr = solve((n - 1) * g * l / V**2 - delta_alpha_T_q_from_V, find)[0]
    return expr.subs(kwargs.items())


def delta_alpha_T_q_from_C_W_and_phi(find, **kwargs):
    """(130) Change in the angle of attack of the tail caused by pitch rate q"""
    delta_alpha_T_q, n, C_W, phi = symbols("delta_alpha_T_q, n, C_W, phi")
    expr = solve((n - 1) * C_W * phi - delta_alpha_T_q, find)[0]
    return expr.subs(kwargs.items())


#print(delta_alpha_T_q_from_C_W_and_phi("n", delta_alpha_T_q=2, n=3, C_W=5, phi=10))


def phi(find, **kwargs):
    """(132) mass parameter phi used for (130)"""
    phi, rho, S, l, m = symbols("phi,rho,S,l,m")
    expr = solve(rho * S * l / (2 * m) - phi, find)[0]
    return expr.subs(kwargs.items())


'''
5.3 - Tailplane Effective Angle of Attack for Manoeuvre
'''


def alpha_T_eff_Manoeuvre(find, **kwargs):
    """(133) The effective tail angle of attack for a manoeuvre, this is different to (17)"""
    alpha_T_eff, eta_alpha, alpha, alpha_0, alpha_S, C_L_T, A_T, e_T, delta_alpha_T_q = symbols(
        "alpha_T_eff, eta_alpha, alpha, alpha_0, alpha_S, C_L_T, A_T, e_T, delta_alpha_T_q")
    expr = solve((1 - eta_alpha) * alpha + eta_alpha * alpha_0 + alpha_S -
                 C_L_T / (np.pi * A_T * e_T) + delta_alpha_T_q - alpha_T_eff, find)[0]
    return expr.subs(kwargs.items())


'''
5.4 - Derivative of Pitch Rate Induced Change in TailAngle of Attack with Angle of Attack
'''


def derivative_delta_alpha_T_q(find, **kwargs):
    """(137) Derivative of the pitch rate induced change in angle of attack of tail wrt aoa"""
    derivative_delta_alpha_T_q, C_L_star_alpha, phi = symbols(
        "derivative_delta_alpha_T_q, C_L_star_alpha, phi")
    expr = solve(C_L_star_alpha * phi - derivative_delta_alpha_T_q, find)[0]
    return expr.subs(kwargs.items())


'''
5.5 - Tailplane Stability Derivative - Manoeuvre Stick Fixed
'''


def C_L_T_alpha_stick_fixed(find, **kwargs):
    """(144) C_L_T stick fixed differentiated wrt alpha"""
    C_L_T_alpha_stick_fixed, k, eta_alpha, phi, C_L_alpha, S, S_T = symbols(
        "C_L_T_alpha_stick_fixed, k, eta_alpha, phi, C_L_alpha, S, S_T")
    expr = solve((k * (1 - eta_alpha) + k * phi * C_L_alpha) /
                 (1 - k * phi * S_T / S) - C_L_T_alpha_stick_fixed, find)[0]
    return expr.subs(kwargs.items())


def H_m_stick_fixed(find, **kwargs):
    """(145) Manoeuvre margin stick fixed"""
    H_m_stick_fixed, h, h_0, C_L_T_alpha, C_L_star_alpha, K = symbols(
        "H_m_stick_fixed, h, h_0, C_L_T_alpha, C_L_star_alpha, K")
    expr = solve(h_0 - h + K * C_L_T_alpha /
                 C_L_star_alpha - H_m_stick_fixed, find)[0]
    return expr.subs(kwargs.items())


'''
5.6 - Tailplane Stability Derivative - Manoeuvre Stick Free
'''


def C_L_T_alpha_stick_free_eqn(find, **kwargs):
    """(152) C_L_T stick free differentiated wrt alpha"""
    C_L_T_alpha_stick_free, k_bar, eta_alpha, phi, C_L_alpha, S, S_T = symbols(
        "C_L_T_alpha_stick_free, k_bar, eta_alpha, phi, C_L_alpha, S, S_T")
    expr = solve((k_bar * (1 - eta_alpha) + k_bar * phi * C_L_alpha) /
                 (1 - k_bar * phi * S_T / S) - C_L_T_alpha_stick_free, find)[0]
    return expr.subs(kwargs.items())


def H_m_stick_free_eqn(find, **kwargs):
    """(153) Manoeuvre margin stick free"""
    H_m_stick_free, h, h_0, C_L_T_alpha, C_L_star_alpha, K = symbols(
        "H_m_stick_free, h, h_0, C_L_T_alpha, C_L_star_alpha, K")
    expr = solve(h_0 - h + K * C_L_T_alpha /
                 C_L_star_alpha - H_m_stick_free, find)[0]
    return expr.subs(kwargs.items())
