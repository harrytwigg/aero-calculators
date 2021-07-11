from sympy import *
from sympy.vector import CoordSys3D
from matplotlib import pyplot as plt
from tabulate import tabulate
import string


class Truss():
    """Truss member class right is x positive, upwards is y positive"""
    N = CoordSys3D('N')

    def __init__(self, points, members):
        self.points = points
        self.members = members
        self.point_data = []
        self.symbol_array = []
        self.eqn_array = []
        self.resolved_forces = []
        self.truss_element_table = []
        self.has_generated_eqn_terms = False
        self.applied_point = ''
        self.applied_force = 0

        # Generate equations for all the points
        for point in points:
            terms_at_point = []

            for joint in members:

                # If joint contains the desired point

                # At position 0
                if point[0] == joint[0]:

                    point_a = point
                    point_b = 0

                    # Find the coordinates of point_b
                    for point2 in points:
                        if point2[0] == joint[1]:
                            point_b = point2
                            break

                    terms_at_point.append(
                        self.get_joint_terms(point_a, point_b))

                # Alternatively at position 1
                elif point[0] == joint[1]:

                    point_a = point
                    point_b = 0

                    # Find the coordinates of point_b
                    for point2 in points:
                        if point2[0] == joint[0]:
                            point_b = point2
                            break

                    terms_at_point.append(
                        self.get_joint_terms(point_a, point_b))

            # Check to see if we need a wall term
            if len(point) == 4:
                if point[3] == "W":
                    terms_at_point.append(("wall term", self.add_symbol(
                        "Rx_" + point[0]) * Truss.N.i + self.add_symbol("Ry_" + point[0]) * Truss.N.j))

            self.point_data.append(terms_at_point)

        # Solve the problem
        self.original_eqn_array = []
        for point in self.points:
            self.original_eqn_array.append(
                self.generate_eqn_at_point(point[0]))

    def add_force_at_point(self, point, force):
        """Takes a force vector at a point and inserts it into the equations"""
        self.applied_point = point
        self.applied_force = force
        for i in range(len(self.points)):
            if point == self.points[i][0]:
                self.original_eqn_array[i] = self.original_eqn_array[i] + force
                return

    def solve_equation_system(self):
        """Solves the system of equations"""
        self.resolved_forces.clear()

        self.eqn_array.clear()
        for vector_eqn in self.original_eqn_array:
            self.eqn_array.append(vector_eqn.coeff(Truss.N.i))
            self.eqn_array.append(vector_eqn.coeff(Truss.N.j))
        self.has_generated_eqn_terms = True

        print("""
Resolving Forces:
        """)

        solvable = True
        while solvable:
            solvable = False  # See if we can set to true whilst looping through the equations
            for i in range(len(self.eqn_array)):
                if i % 2 == 0:
                    point_pos = int(i / 2)
                    axis = 'i:'
                else:
                    point_pos = int((i - 1) / 2)
                    axis = 'j:'

                if len(self.eqn_array[i].free_symbols) == 1:
                    solvable = True
                    starting_eqn = self.eqn_array[i]
                    symbol_solved = list(self.eqn_array[i].free_symbols)[0]
                    resolved_force = solve(self.eqn_array[i])[0]
                    self.resolved_forces.append(
                        (symbol_solved, resolved_force))
                    self.sub_all(symbol_solved, resolved_force)
                    print(self.points[point_pos][0], axis, "0 =", starting_eqn,
                          " |  resolving", symbol_solved, "==>", resolved_force)

    def find_truss_elements(self, A, E):
        """Calculate the turn element table under a state of unixial stress"""
        self.truss_element_table.clear()

        def should_add_member(member):
            for row in self.truss_element_table:
                if row[0] == member:
                    return False
            return True

        def get_force(force_symbol):
            for force in self.resolved_forces:
                if force[0] == force_symbol:
                    return force[1]

        for point in self.point_data:
            for term in point:
                if term[0] != "wall term":
                    if string.ascii_uppercase.index(term[2][0]) > string.ascii_uppercase.index(term[3][0]):
                        member_name = term[3][0] + term[2][0]
                    else:
                        member_name = term[2][0] + term[3][0]

                    if should_add_member(member_name):
                        force_symbol = Symbol("F_" + member_name)
                        force_value = get_force(force_symbol)
                        length = term[0]
                        delta = float(force_value * length / (E * A))
                        self.truss_element_table.append(
                            [member_name, force_value, length, f'{delta:g}'])

        table_headers = ["Member", "Force [N]", "Length [m]", "δ [m]"]
        return self.truss_element_table, table_headers

    def principle_virtual_work(self, table=None):
        """Find the deflection u of the original applied force point"""
        if table == None:
            table = self.truss_element_table

        sum_F_delta = 0
        for row in table:
            sum_F_delta += row[1] * float(row[3])

        u = sum_F_delta / self.applied_force.magnitude()

        return u

    @staticmethod
    def displacement_arbitrary_point(real_elements, virtual_elements, virtual_force):
        """Returns the displacement u for the arbitrary point applied to the virtual truss"""
        merged_table = []
        sum_f_internal_delta = 0

        for row_real in real_elements:
            member = row_real[0]
            delta = row_real[3]
            f_internal = 0
            f_internal_delta = 0
            for row_virtual in virtual_elements:
                if row_virtual[0] == member:
                    f_internal = sympify(row_virtual[1])
                    break
            f_internal_delta = sympify(f_internal) * sympify(delta)
            sum_f_internal_delta += f_internal_delta

            merged_table.append([member, delta, f_internal, f_internal_delta])

        table_headers = ["Member", "δ Actual[m]",
                         "F* Virtual [N]", "F* x δ Actual"]

        print("""
Merged Truss Data:
        """)

        print(tabulate(merged_table, headers=table_headers, disable_numparse=True))

        return sum_f_internal_delta / virtual_force.magnitude()

    def print_resolved_forces(self):
        print("""
Resolved Forces:
        """)
        for resolved in self.resolved_forces:
            print(resolved[0], "=", resolved[1])

    def sub_all(self, symbol, value):
        eqn_array_cache = []
        for eqn in self.eqn_array:
            eqn_array_cache.append(eqn.subs(symbol, value))
        self.eqn_array = eqn_array_cache

    def add_symbol(self, name):
        symbol = Symbol(name)

        # Check if symbol already stored
        can_add = True
        for old_symbol in self.symbol_array:
            if old_symbol == symbol:
                can_add = False
                break
        if can_add:
            self.symbol_array.append(symbol)
        return symbol

    def get_joint_terms(self, point_a, point_b):
        """Returns the length and force term at a certian point"""
        vector = (point_b[1] - point_a[1]) * Truss.N.i + \
            (point_b[2] - point_a[2]) * Truss.N.j
        length = vector.magnitude()

        # Alphabetise force names
        underscore = ""
        if string.ascii_uppercase.index(point_a[0]) > string.ascii_uppercase.index(point_b[0]):
            underscore = point_b[0] + point_a[0]
        else:
            underscore = point_a[0] + point_b[0]

        force_symbol = self.add_symbol("F_" + underscore)
        force_term = force_symbol * vector / vector.magnitude()  # Get the unit vector
        return (length, force_term, point_a, point_b)

    def plot(self):
        """Plot truss structure"""

        fig, ax = plt.subplots()

        max_x = "unset"
        max_y = "unset"
        min_x = "unset"
        min_y = "unset"

        for point in self.points:
            if max_x == "unset" or point[1] > max_x:
                max_x = point[1]

            if max_y == "unset" or point[2] > max_y:
                max_y = point[2]

            if min_x == "unset" or point[1] < min_x:
                min_x = point[1]

            if min_y == "unset" or point[2] < min_y:
                min_y = point[2]

        b = 0.1*max(max_x - min_x, max_y - min_y)
        point_size = (max_x - min_x) / 100

        for joint in self.members:
            x_values = []
            y_values = []
            for point in self.points:
                if joint[0] == point[0]:
                    x_values.append(point[1])
                    y_values.append(point[2])
                elif joint[1] == point[0]:
                    x_values.append(point[1])
                    y_values.append(point[2])
            plt.plot(x_values, y_values, color="grey")

        for point in self.points:
            if len(point) == 4:
                if point[3] == "W":
                    ax.add_patch(plt.Circle(
                        (point[1], point[2]), point_size, color="red"))
            else:
                ax.add_patch(plt.Circle(
                    (point[1], point[2]), point_size, color="black"))
            plt.text(point[1] + 2 * point_size,
                     point[2] + 2 * point_size, point[0])

        plt.xlim(xmin=min_y-b, xmax=max_x+b)
        plt.ylim(ymin=min_y-b, ymax=max_y+b)
        plt.show()

    def print_joint_pairs(self):
        for eqn in self.point_data:
            for term in eqn:
                if term[0] != "wall term" and term[0] != "loaded point":
                    print(term[2], term[3])

    def generate_eqn_at_point(self, point):
        point_index = 0
        for i in range(len(self.point_data)):
            if point == self.points[i][0]:
                point_index = i
                break

        x = Symbol("x")
        sum = 0 * Truss.N.i + 0 * Truss.N.j
        for term in self.point_data[point_index]:
            sum += term[1]
        return sum

    def print_all_eqns(self):
        print("""
Equations:
        """)
        if self.has_generated_eqn_terms:
            for i in range(len(self.points)):
                print(self.points[i][0], "i:", "0 =", self.eqn_array[2*i])
                print(self.points[i][0], "j:", "0 =", self.eqn_array[2*i+1])
        else:
            for i in range(len(self.points)):
                print(self.points[i][0], "i:", "0 =",
                      self.original_eqn_array[i].coeff(Truss.N.i))
                print(self.points[i][0], "j:", "0 =",
                      self.original_eqn_array[i].coeff(Truss.N.j))

    def live_eqn_at_point(self, point):
        i = 0
        for current_point in self.points:
            if current_point[0] == point:
                return self.eqn_array[2*i], self.eqn_array[2*i+1]
            i += 1


def pvw_lecture_15():
    """Example question lecture 15, uses real external force"""
    points = [("B", 0, 2, "W"), ("C", 2, 2), ("D", 4, 2),
              ("F", 0, 0, "W"), ("G", 2, 0)]
    members = [("B", "C"), ("C", "D"), ("D", "G"),
               ("F", "G"), ("F", "C"), ("C", "G")]

    N = CoordSys3D('N')
    force_term = - 10000 * N.j

    truss = Truss(points, members)
    truss.add_force_at_point("D", force_term)
    truss.print_all_eqns()
    truss.solve_equation_system()
    truss.print_resolved_forces()

    table, headers = truss.find_truss_elements(1000 * 10**-6, 200 * 10**9)
    print("""
Truss Elements:
    """)

    print(tabulate(table, headers, disable_numparse=True))
    u = truss.principle_virtual_work()

    print("")
    print("Displacement u at the point D is", u, "as float", u.evalf())

    truss.plot()


def virtual_forces_lecture_16():
    """Example question lecture 16, uses a virtual external force"""

    points = [("B", 0, 2, "W"), ("C", 2, 2), ("D", 4, 2),
              ("F", 0, 0, "W"), ("G", 2, 0)]
    members = [("B", "C"), ("C", "D"), ("D", "G"),
               ("F", "G"), ("F", "C"), ("C", "G")]

    A = 1000 * 10**-6
    E = 200 * 10**9

    N = CoordSys3D('N')
    real_force = - 10000 * N.j
    real_point = "D"

    virtual_force = - 1 * N.j
    virtual_point = "G"

    print("""
=====================
    Real Truss
=====================
    """)

    # Apply analysis on real truss in same way as before
    real_truss = Truss(points, members)
    real_truss.add_force_at_point(real_point, real_force)
    real_truss.print_all_eqns()
    real_truss.solve_equation_system()
    real_truss.print_resolved_forces()
    real_elements, headers = real_truss.find_truss_elements(A, E)
    real_truss.plot()
    print("""
Real Truss Elements:
    """)

    print(tabulate(real_elements, headers=headers, disable_numparse=True))

    print("""
=====================
    Virtual Truss
=====================
    """)

    # Make new truss and apply virtual force at point G and solve
    virtual_truss = Truss(points, members)
    virtual_truss.add_force_at_point(virtual_point, virtual_force)
    virtual_truss.print_all_eqns()
    virtual_truss.solve_equation_system()
    virtual_truss.print_resolved_forces()
    virtual_elements = virtual_truss.find_truss_elements(A, E)[0]

    print("""
Virtual Truss Elements:
    """)

    print(tabulate(virtual_elements, headers=headers, disable_numparse=True))

    u = Truss.displacement_arbitrary_point(
        real_elements, virtual_elements, virtual_force)
    print("")
    print("Displacement u at the point G is", u, "as float", u.evalf())


if __name__ == "__main__":
    """Run functions here"""
    # pvw_lecture_15()
    virtual_forces_lecture_16()
