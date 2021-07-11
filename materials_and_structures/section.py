import matplotlib.pyplot as plt
from math import atan2, sin, cos, sqrt, pi, degrees
import numpy as np


def area(points):
    """Area of cross-section using greens theorem

    points - an array of (y, z) point pairs"""
    if points[0] != points[-1]:
        points = points + points[:1]
    y = [point[0] for point in points]
    z = [point[1] for point in points]

    sum = 0

    for i in range(len(points) - 1):
        sum += y[i] * z[i+1] - y[i+1] * z[i]
    return sum / 2


def centroid(points):
    """Location of the centroid

    points - an array of (y, z) point pairs"""

    if points[0] != points[-1]:
        points = points + points[:1]
    y = [point[0] for point in points]
    z = [point[1] for point in points]

    y_sum = z_sum = 0
    A = area(points)

    for i in range(len(points) - 1):
        y_sum += (y[i] + y[i+1]) * (y[i]*z[i+1] - y[i+1] * z[i])
        z_sum += (z[i] + z[i+1]) * (y[i]*z[i+1] - y[i+1] * z[i])

    return y_sum/(6*A), z_sum/(6*A)


def inertia(points):
    """Moments and product of inertia about centroid

    points - an array of (y, z) point pairs"""

    if points[0] != points[-1]:
        points = points + points[:1]
    y = [c[0] for c in points]
    z = [c[1] for c in points]

    syy = szz = syz = 0
    A = area(points)
    y_c, z_c = centroid(points)

    for i in range(len(points) - 1):
        syy += (z[i]**2 + z[i]*z[i+1] + z[i+1]**2)*(y[i]*z[i+1] - y[i+1]*z[i])
        szz += (y[i]**2 + y[i]*y[i+1] + y[i+1]**2)*(y[i]*z[i+1] - y[i+1]*z[i])
        syz += (y[i]*z[i+1] + 2*y[i]*z[i] + 2*y[i+1]*z[i+1] +
                y[i+1]*z[i])*(y[i]*z[i+1] - y[i+1]*z[i])

    return syy/12 - A*z_c**2, szz/12 - A*y_c**2, syz/24 - A*y_c*z_c


def principal(Iyy, Izz, Iyz):
    """Principal moments of inertia and orientation

    points - an arraz of (y, z) point pairs"""

    avg = (Iyy + Izz) / 2
    diff = (Iyy - Izz) / 2
    I1 = avg + sqrt(diff**2 + Iyz**2)
    I2 = avg - sqrt(diff**2 + Iyz**2)
    theta = atan2(-Iyz, diff)/2
    return I1, I2, theta


def summary(points):
    """Text summary of cross-sectional properties

    points - an array of (y, z) point pairs"""

    a = area(points)
    y_c, z_c = centroid(points)
    Iyy, Izz, Iyz = inertia(points)
    I1, I2, theta = principal(Iyy, Izz, Iyz)
    summ = """Area
  A = {}
Centroid
  y_c = {}
  z_c = {}
Moments and product of inertia
  Iyy = {}
  Izz = {}
  Iyz = {}
Principal moments of inertia and direction
  I1 = {}
  I2 = {}
  θ︎ = {}°""".format(a, y_c, z_c, Iyy, Izz, Iyz, I1, I2, degrees(theta))
    return summ


def outline(points, size=(8, 8), dpi=100):
    """Draw an outline of the cross-section with centroid and principal axes

    points - an array of (y, z) point pairs"""

    if points[0] != points[-1]:
        points = points + points[:1]
    y = [c[0] for c in points]
    z = [c[1] for c in points]

    # Get the bounds of the cross-section
    min_y = min(y)
    max_y = max(y)
    min_z = min(z)
    max_z = max(z)

    # Whitespace border is 5% of the larger dimension
    b = .05*max(max_z - min_y, max_y - min_y)

    # Get the properties needed for the centroid and principal axes
    y_c, z_c = centroid(points)
    i = inertia(points)
    p = principal(*i)

    # Principal axes extend 10% of the minimum dimension from the centroid
    length = min(max_z-min_z, max_y-min_y)/10
    a1y = [y_c - length*cos(p[2]), y_c + length*cos(p[2])]
    a1z = [z_c - length*sin(p[2]), z_c + length*sin(p[2])]
    a2y = [y_c - length*cos(p[2] + pi/2), y_c + length*cos(p[2] + pi/2)]
    a2z = [z_c - length*sin(p[2] + pi/2), z_c + length*sin(p[2] + pi/2)]

    fig, ax = plt.subplots(figsize=size)
    ax.plot(z, y, 'k*-', lw=2)
    ax.plot(a1z, a1y, '-', color='#0072B2', lw=2)     # blue
    ax.plot(a2z, a2y, '-', color='#D55E00')           # vermillion
    ax.plot(z_c, y_c, 'ko', mec='k')

    ax.set_xlabel("z")
    ax.set_ylabel("y")
    ax.set_aspect('equal')
    ax.xaxis.tick_top()

    plt.xlim(xmin=min_z-b, xmax=max_z+b)
    plt.ylim(ymin=min_y-b, ymax=max_y+b)

    plt.gca().invert_yaxis()

    plt.show()


if __name__ == "__main__":
    shape = [(0, 0), (120, 0), (120, 120), (110, 120), (110, 10), (0, 10)]

    print(summary(shape))
    outline(shape, size=(8, 6))
