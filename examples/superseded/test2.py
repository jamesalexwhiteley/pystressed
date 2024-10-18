from scipy.stats import norm
from scipy.special import erf, erfinv
import scipy.optimize as opt
from shapely.geometry import LineString, Point, MultiPoint
import cvxpy as cp

import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# ================== Q2 ================== #
def plot_line(start, end, label, alpha, linestyle, extension):
    delta_x = end[0] - start[0]
    delta_y = end[1] - start[1]
    extended_end = (start[0] + delta_x * extension, start[1] + delta_y * extension)
    plt.plot([start[0], extended_end[0]], [start[1], extended_end[1]], 'k', label=label, alpha=alpha, linestyle=linestyle)
    return (start, extended_end)

def get_unique_intersections(lines):
    # find all unique intersections between lines
    intersections = {
        (point.x, point.y)
        for i in range(len(lines))
        for j in range(i + 1, len(lines))
        for intersection in [lines[i].intersection(lines[j])]
        if isinstance(intersection, Point) and not intersection.is_empty
        or isinstance(intersection, MultiPoint) and not intersection.is_empty
        for point in (intersection if isinstance(intersection, MultiPoint) else [intersection])
    }
    return intersections

def plot_magnel(P, Z1, Z2, A, e, line_extension=1.5):
    # plot Magnel diagram        
    Z = [Z1, Z1, Z2, Z2]    
    label = ['c1', 't1', 'c2', 't2'] 
    alpha = [1, .4, 1, .4]
    linestyle = ['-', '-', '--', '--'] 

    plt.figure()
    lines = []
    for i, z in enumerate(Z):
        start, end = plot_line((0, -z / A), (1 / P, e[i]), label[i], alpha[i], linestyle[i], line_extension)
        lines.append(LineString([start, end]))

    print("Intersection Points:")
    for point in get_unique_intersections(lines):
        try:
            p0 = 1/point[0] 
        except:
            p0 = point[0] 
        print(f"P={float(p0)*1e-3:.0f} kN, e={float(point[1]):.0f}")
    
    plt.xlabel("1/P")
    plt.ylabel("e")
    plt.gca().invert_yaxis()
    plt.grid(True)
    plt.legend()
    plt.show()

def stress_inequality(f, Z, M, A, P):
    return -Z/A + f*Z/P + M/P

def calc_stress_inequalitieis(fc, ft, Z1, Z2, Mmax, Mmin, A, P):
    c1 = stress_inequality(fc, Z1, Mmax, A, P)
    t1 = stress_inequality(ft, Z1, Mmin, A, P)
    c2 = stress_inequality(fc, Z2, Mmin, A, P)
    t2 = stress_inequality(ft, Z2, Mmax, A, P)
    # print(f"e TOP (comp): {c1:.0f} mm")
    # print(f"e TOP (tens): {t1:.0f} mm")
    # print(f"e BOT (comp): {c2:.0f} mm")
    # print(f"e BOT (tens): {t2:.0f} mm")
    return (c1, t1, c2, t2)

# parameters 
P = 4800e3 # N 
A = 0.64e6 # mm2
fc = 15 # N/mm2 
ft = 1 # N/mm2
Z1 = -309e6 # mm3
Z2 = 265e6 # mm3
Mmax = 5000e6 # Nmm
Mmin = 1000e6 # Nmm

e = calc_stress_inequalitieis(fc, ft, Z1, Z2, Mmax, Mmin, A, P)
plot_magnel(P, Z1, Z2, A, e)