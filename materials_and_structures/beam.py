from sympy import *
import numpy as np
from matplotlib import pyplot as plt
import numpy.polynomial.polynomial as fitter

from section import *

# Symbols used in formulas
x, R_a, Q_y, M_a, M_z, E_symbol = symbols("x, R_a, Q_y, M_a, M_z, E")


class Beam():
    """Calculate various properties of a beam

    Downwards is force +ve and clockwise is moment +ve"""

    def __init__(self, length, I_zz):
        self.length = length
        self.I_zz = I_zz

        # Actual force values or distributions and their positions
        self.force_terms = []

    def add_point_load(self, size, x_loc):
        """Adds a point load, downwards is positive"""
        self.force_terms.append((size, x_loc))

    def add_distributed_load(self, distribution, start, end):
        self.force_terms.append((distribution, start, end))

    def get_force_terms(self):
        """Returns the distributed force terms summed using heavyside functions"""
        force_eqn = 0
        for term in self.force_terms:
            if len(term) == 2:
                force = term[0]
                x_loc = term[1]
                force_eqn += Heaviside(x - x_loc) * force
            elif len(term) == 3:
                distribution = term[0]
                start = term[1]
                end = term[2]
                force_eqn += Heaviside(x - start) * (distribution * (x - start) * (Heaviside(
                    x - start) - Heaviside(x - end)) + Heaviside(x - end) * distribution * (end - start))

        return force_eqn

    def get_R_a(self):
        resolve_vertical_eqn = R_a - self.get_force_terms()
        resolve_vertical_eqn = resolve_vertical_eqn.subs(x, self.length)

        # Substitute 0 for undefined heavyside values
        resolve_vertical_eqn = resolve_vertical_eqn.subs(Heaviside(0), 0)

        return solve(resolve_vertical_eqn, R_a)[0]

    def shear_force(self):
        """Get the cutface shear force Q_y(x)"""
        return solve(self.get_force_terms() + Q_y - self.get_R_a(), Q_y)[0]

    def get_M_a(self):
        M_a = 0
        for term in self.force_terms:
            if len(term) == 2:
                force = term[0]
                x_loc = term[1]
                M_a -= x_loc * force
            elif len(term) == 3:
                distribution = term[0]
                start = term[1]
                end = term[2]
                M_a -= (end - start) * distribution * \
                    (start + (end - start) / 2)

        return M_a

    def get_moment_terms(self):
        """Returns the distributed moments terms summed using heavyside functions measured from the cut"""
        moment_eqn = 0
        for term in self.force_terms:
            if len(term) == 2:
                force = term[0]
                x_loc = term[1]
                moment_eqn -= Heaviside(x - x_loc) * (x - x_loc) * force
            elif len(term) == 3:
                distribution = term[0]
                start = term[1]
                end = term[2]
                moment_eqn -= Heaviside(x - start) * ((Heaviside(x - start) - Heaviside(x - end)) * (x - start) **
                                                      2 * 0.5 * distribution + Heaviside(x - end) * (end - start) * distribution * (x - (start - (end-start)/2)))

        return moment_eqn

    def bending_moment(self):
        """Get the cutface bending moment M_z(x)"""
        return solve(self.get_R_a() * x + self.get_moment_terms() + self.get_M_a() - M_z, M_z)[0]

    def bending_stress(self, M_y, y, z, I_yy, I_yz, x_value):
        """Get the bending stress σ_xx equation for x position of the cutface"""
        M_z = self.bending_moment()
        M_z = M_z.subs(x, x_value)
        return ((M_z * I_yy - M_y * I_yz) * y + (M_y * self.I_zz - M_z * I_yz)) * z / (I_yy * self.I_zz - I_yz**2)

    def deflection_v(self, E, M_y, I_yy, I_yz, M_z=None, step=0.001, poly_approx=100, should_plot=False):
        if M_z == None:
            M_z = self.bending_moment()
        v_2_d = - (M_z * I_yy - M_y * I_yz) / \
            (E * (I_yy * self.I_zz - I_yz ** 2))
        x_values = np.arange(0, self.length, step)

        def gen_y_values(function):
            y_values = []
            for value in x_values:
                y_values.append(float(function.subs(x, value).evalf()))
            return np.array(y_values)

        y_values = gen_y_values(v_2_d)
        coeffs = np.flip(fitter.polyfit(x_values, y_values, poly_approx))
        poly_v_2_d = Poly(coeffs, x).as_expr()
        poly_v_1_d = integrate(poly_v_2_d, x)
        poly_v_0_d = integrate(poly_v_1_d, x)

        if should_plot:
            plt.plot(x_values, y_values)
            plt.plot(x_values, gen_y_values(poly_v_2_d))
            plt.plot(x_values, gen_y_values(poly_v_1_d))
            plt.plot(x_values, gen_y_values(poly_v_0_d))
            plt.show()

        return poly_v_0_d

    def deflection_w(self, E, M_y, I_yy, I_yz, M_z=None, step=0.001, poly_approx=100, should_plot=False):
        if M_z == None:
            M_z = self.bending_moment()
        w_2_d = - (M_y * self.I_zz - M_z * I_yz) / \
            (E * (I_yy * self.I_zz - I_yz ** 2))
        x_values = np.arange(0, self.length, step)

        def gen_y_values(function):
            y_values = []
            for value in x_values:
                y_values.append(float(function.subs(x, value).evalf()))
            return np.array(y_values)

        y_values = gen_y_values(w_2_d)
        coeffs = np.flip(fitter.polyfit(x_values, y_values, poly_approx))
        poly_w_2_d = Poly(coeffs, x).as_expr()
        poly_w_1_d = integrate(poly_w_2_d, x)
        poly_w_0_d = integrate(poly_w_1_d, x)

        if should_plot:
            plt.plot(x_values, y_values)
            plt.plot(x_values, gen_y_values(poly_w_2_d))
            plt.plot(x_values, gen_y_values(poly_w_1_d))
            plt.plot(x_values, gen_y_values(poly_w_0_d))
            plt.show()

        return poly_w_0_d
    
    def neutral_axis(self, M_y, I_yy, I_yz, I_zz=None, M_z=None, x_value=None):
        """Returns the angle of the neutral axis in degrees"""
        if I_zz == None:
            I_zz = self.I_zz
        if M_z == None:
            M_z = self.bending_moment()
            M_z = M_z.subs(x, x_value)
        return degrees(np.arctan(- (M_y * I_zz - M_z * I_yz) / (M_z*I_yy - M_y*I_yz)))

    def shear_flow(self, y_bar, z_bar, A, Q_z, I_yy, I_yz, I_zz=None, Q_y=None):
        """Returns the shear flow q in the cross section"""
        if I_zz == None:
            I_zz = self.I_zz
        if Q_y == None:
            Q_y = self.shear_force()
        return ((Q_y * I_yy - Q_z * I_yz) * y_bar * A + (Q_z * I_zz - Q_y * I_yz) * z_bar * A) / (I_yy * I_zz - I_yz ** 2)


def bending_shear_force_lecture_1():
    """Example shear force from lecture 1"""
    beam = Beam(1, 24.2 * 10**(-6))
    beam.add_distributed_load(10000, 0, 1)
    beam.add_point_load(5000, 0.5)
    print("R_a", beam.get_R_a())
    print("M_a", beam.get_M_a())
    #print(beam.neutral_axis(0, 6.6 * 10 **-6, 7.2 * 10 **-6))

    bending_stress = beam.bending_stress(0, -0.01, -0.1, 6.6 * 10 **-6, 7.2 * 10 **-6, x_value=0)
    print("bending_stress", bending_stress)

    
    bending_moment = beam.bending_moment()
    shear_force = beam.shear_force()
    x_points = []
    bending_moment_points = []
    shear_force_points = []

    # Enables floating point increments in loops
    def drange(start, stop, step):
        while start < stop:
            yield start
            start += step

    for i in drange(0, 1, 0.01):
        x_points.append(i)
        bending_moment_points.append(bending_moment.subs(x, i))
        shear_force_points.append(shear_force.subs(x, i))

    plt.plot(x_points, bending_moment_points)
    plt.plot(x_points, shear_force_points)

    plt.legend(["Bending Moment", "Shear Force"])
    plt.show()
    


def deflection_lecture_2():
    """Example deflection from lecture 2"""
    shape = [(0, 0), (0.020, 0), (0.020, 0.100),
             (0.200, 0.100), (0.200, 0.120), (0, 0.120)]
    I_yy, I_zz, I_yz = inertia(shape)
    # outline(shape)
    print(I_yy, I_zz, I_yz)

    E = 70 * 10**9
    M_y = 0

    beam = Beam(1, I_zz)
    beam.add_distributed_load(10000, 0, 1)
    v_poly = beam.deflection_v(
        E, M_y, I_yy, I_yz, should_plot=True, step=0.01, poly_approx=10)
    v_tip = v_poly.subs(x, 1) * 1000
    print("The deflection v at the beam tip is", v_tip, "mm")

    w_poly = beam.deflection_w(
        70 * 10**9, 0, I_yy, I_yz, should_plot=True, step=0.01, poly_approx=10)
    w_tip = w_poly.subs(x, 1) * 1000
    print("The deflection w at the beam tip is", w_tip, "mm")


def max_shear_stress_lecture_4_question_1():
    """Example 1 of 2: finding the max bending shear stress in a cross section"""

    # Whole cross section properties
    shape = [(0, 0), (0.020, 0), (0.020, 0.100),
             (0.200, 0.100), (0.200, 0.120), (0, 0.120)]
    #outline(shape)
    A = area(shape)
    y_c, z_c = centroid(shape)
    print("z_c", z_c, "y_c", y_c)
    I_yy, I_zz, I_yz = inertia(shape)
    print("I_yy, I_zz, I_yz")
    print(I_yy, I_zz, I_yz)

    beam = Beam(1, I_zz)
    beam.add_distributed_load(10000, 0, 1)
    #plot(beam.shear_force(), xlim=(0, 1))
    # From graph we can see Q_y max is at x = 0
    Q_z = 0
    Q_y = beam.shear_force().subs(x, 0)

    # Then we start at the free edge and work around the strcture
    S_1, S_2= symbols("S_1, S_2")

    thickness = 0.02

    y_bar_12 = thickness / 2 - y_c
    print("y_bar_12", y_bar_12)
    z_bar_12 = S_1 / 2 - z_c
    print("z_bar_12", z_bar_12)

    A_12 = S_1 * thickness

    q_12 = expand(beam.shear_flow(y_bar_12, z_bar_12, A_12, Q_z, I_yy, I_yz, Q_y=Q_y))
    tau_12 = q_12 / thickness
    print(tau_12)

    # End of section is S_1 = 0.110
    tau_12_max = tau_12.subs(S_1, 0.11)
    print(tau_12_max)

    y_bar_32 = 0.2 - S_2 / 2 - y_c
    print("y_bar_32", y_bar_32)
    z_bar_32 = 0.12 - thickness / 2 - z_c
    print("z_bar_32", z_bar_32)

    A_32 = S_2 * thickness

    q_32 = expand(beam.shear_flow(y_bar_32, z_bar_32, A_32, Q_z, I_yy, I_yz, Q_y=Q_y))
    tau_32 = q_32 / thickness

    # End of section is S_2 = 0.190
    tau_32_max = tau_32.subs(S_2, 0.19)
    print(tau_32_max)

    '''Find min and max bending shear stress using the derivatives for each section
    and see which one is bigger'''


if __name__ == "__main__":
    x = symbols("x")
    bending_shear_force_lecture_1()
    #deflection_lecture_2()
    #max_shear_stress_lecture_4_question_1()
