from shapely.geometry import LineString, Point, MultiPoint
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# Author: James Whiteley (github.com/jamesalexwhiteley)

def plot_line(start, end, label, alpha, linestyle, extension):
    delta_x = end[0] - start[0]
    delta_y = end[1] - start[1]
    extended_end = (start[0] + delta_x * extension, start[1] + delta_y * extension)
    # plt.plot([start[0], extended_end[0]], [start[1], extended_end[1]], 'k', label=label, alpha=alpha, linestyle=linestyle)
    plt.plot([start[1], extended_end[1]], [start[0], extended_end[0]], 'k', label=label, alpha=alpha, linestyle=linestyle)
    return (start, extended_end)

def get_unique_intersections(lines):
    # find all unique intersections between lines
    intersections = set()
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            intersection = lines[i].intersection(lines[j])
            if isinstance(intersection, Point):
                if not intersection.is_empty:
                    intersections.add((intersection.x, intersection.y))
            elif isinstance(intersection, MultiPoint):
                for point in intersection:
                    if not intersection.is_empty:
                        intersections.add((point.x, point.y))
    return intersections

def plot_magnel(P, Z1, Z2, A, e, line_ext=2.0):
    # plot Magnel diagram        
    Z = [Z1, Z1, Z2, Z2]    
    label = ['c1', 't1', 'c2', 't2'] 
    alpha = [1, .4, 1, .4]
    linestyle = ['-', '-', '--', '--'] 

    plt.figure()
    lines = []
    for i, z in enumerate(Z):
        start, end = plot_line((0, -z / A), (1 / P, e[i]), label[i], alpha[i], linestyle[i], line_ext)
        lines.append(LineString([start, end]))
    
    plt.xlabel("1/P")
    plt.ylabel("e")
    plt.gca().invert_yaxis()
    plt.grid(True)
    plt.legend()
    plt.show()

    return get_unique_intersections(lines) 

def e_inequality(f, Z, M, A, P):
    return -Z/A + f*Z/P + M/P

def plot_magnel(state, line_ext=2.0):
    fc, ft, A, Ztop, Zbot, Mmax, Mmin = state.fc, state.ft, state.A, state.Ztop, state.Zbot, state.Mmax, state.Mmin
    P = fc*A/2 # arbitrary, but fc*A/2 suitable  
    c1 = e_inequality(fc, -Ztop, Mmax, A, P)
    t1 = e_inequality(ft, -Ztop, Mmin, A, P)
    c2 = e_inequality(fc, Zbot, Mmin, A, P)
    t2 = e_inequality(ft, Zbot, Mmax, A, P)
    # print(f"e TOP (comp): {c1:.0f} mm")
    # print(f"e TOP (tens): {t1:.0f} mm")
    # print(f"e BOT (comp): {c2:.0f} mm")
    # print(f"e BOT (tens): {t2:.0f} mm")
    e = (c1, t1, c2, t2)
    return plot_magnel(P, -Ztop, Zbot, A, e, line_ext=line_ext)

