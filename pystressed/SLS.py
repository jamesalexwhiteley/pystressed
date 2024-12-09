# from shapely.geometry import LineString, Point, MultiPoint
import matplotlib.pyplot as plt
import cvxpy as cp

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# Author: James Whiteley (github.com/jamesalexwhiteley)

# def get_unique_intersections(lines):
#     # find all unique intersections between lines
#     intersections = {
#         (point.x, point.y)
#         for i in range(len(lines))
#         for j in range(i + 1, len(lines))
#         for intersection in [lines[i].intersection(lines[j])]
#         if isinstance(intersection, Point) and not intersection.is_empty
#         or isinstance(intersection, MultiPoint) and not intersection.is_empty
#         for point in (intersection if isinstance(intersection, MultiPoint) else [intersection])
#     }
#     return intersections

def plot_line(start, end, label, alp, linestyle, extension, color):
    delta_x = end[0] - start[0]
    delta_y = end[1] - start[1]
    extended_end = (start[0] + delta_x * extension, start[1] + delta_y * extension)
    # plt.plot([start[0], extended_end[0]], [start[1], extended_end[1]], label=label, alpha=alp, linestyle=linestyle, color=color)
    plt.plot([start[0], extended_end[0]], [start[1], extended_end[1]], alpha=alp, linestyle=linestyle, color=color)
    return (start, extended_end)

def plot_lines(Z, A, P, e, labels, alp, linestyles, line_ext, color):

    for i, z in enumerate(Z):
        plot_line((0, -z / A), (1 / P, e[i]), labels[i], alp[i], linestyles[i], line_ext, color=color)

    # lines = []
    # for i, z in enumerate(Z):
    #     start, end = plot_line((0, -z / A), (1 / P, e[i]), labels[i], alp[i], linestyles[i], line_ext, color=color)
        # lines.append(LineString([start, end]))

    # if print_intersections: 
    #     print("Intersection Points:")
    #     for point in get_unique_intersections(lines):
    #         if point[0] != 0:
    #             print(f"P={float(1/point[0])*1e-3:.0f} kN, e={float(point[1]):.0f}")

def plot_reference_lines(ax, point, xlim, ylim, color='gray', linestyle='--', linewidth=1, zorder=5):
    """
    Plots dashed lines from the axes to a given point of interest.
    
    Parameters:
        ax (matplotlib.axes.Axes): The axes object to plot on.
        point (tuple): The (x, y) coordinates of the point of interest.
        color (str): Color of the dashed lines.
        linestyle (str): Line style (default is '--').
        linewidth (float): Width of the lines.
        zorder (int): Layer order of the lines.
    """
    # xlim = ax.get_xlim()
    # ylim = ax.get_ylim()
    
    ax.plot([point[0], point[0]], [ylim[0], point[1]], color=color, linestyle=linestyle, linewidth=linewidth, zorder=zorder)
    ax.plot([xlim[0], point[0]], [point[1], point[1]], color=color, linestyle=linestyle, linewidth=linewidth, zorder=zorder)

def _plot_magnel(transfer, P0, e0, service, P1, e1, line_ext, min_max_points):
    
    # plot Magnel diagram      
    Ztop0, Zbot0, A0 = transfer.Ztop, transfer.Zbot, transfer.A
    Z0 = [-Ztop0, -Ztop0, Zbot0, Zbot0]  

    Ztop1, Zbot1, A1 = service.Ztop, service.Zbot, service.A
    Z1 = [-Ztop1, -Ztop1, Zbot1, Zbot1]    
    
    top, bot = 'top', 'bot'
    labels0 = [fr'$f_{{c,{top},tf}}$', fr'$f_{{t,{top},tf}}$', fr'$f_{{c,{bot},tf}}$', fr'$f_{{t,{bot},tf}}$']
    labels1 = [fr'$f_{{c,{top},sv}}$', fr'$f_{{t,{top},sv}}$', fr'$f_{{c,{bot},sv}}$', fr'$f_{{t,{bot},sv}}$'] 
    # alp = [.4, 1, .4, 1]
    # linestyles = [(0, (10, 2, 1, 2)), (0, (10, 2, 1, 2)), '-', '-'] 
    alp = [1, 1, 1, 1]
    linestyles = ['-', '-', '-', '-'] 

    fig, ax = plt.subplots(figsize=(6, 6))

    # account for losses here 
    plot_lines(Z0, A0, P0 / transfer.losses, e0, labels0, alp, linestyles, line_ext, color='k')
    plot_lines(Z1, A1, P1 /  service.losses, e1, labels1, alp, linestyles, line_ext, color='k')

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    if min_max_points != None:
        ax.scatter(min_max_points[0][0], min_max_points[0][1], s=75, color='black', edgecolor='white', zorder=10, label='Minimum prestress')
        plot_reference_lines(ax, min_max_points[0], xlim, ylim)

        ax.scatter(min_max_points[1][0], min_max_points[1][1], s=75, color='white', edgecolor='black', zorder=10, label='Maximum prestress')
        plot_reference_lines(ax, min_max_points[1], xlim, ylim)
    
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel(r"1/P (N$^{-1})$")
    ax.set_ylabel(r"Eccentricity (mm)")
    ax.grid(True, color='grey', linestyle='-', linewidth=0.5, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    plt.savefig('magnel.png', bbox_inches='tight', dpi=600)
    plt.show()

def e_inequality(f, Z, M, A, P):
    return -Z/A + f*Z/P + M/P

def calc_inequalities(state, P):
    
    fc, ft, A, Ztop, Zbot, Mmax, Mmin = state.fc, state.ft, state.A, state.Ztop, state.Zbot, state.Mmax, state.Mmin
    c1 = e_inequality(fc, -Ztop, Mmax, A, P) 
    t1 = e_inequality(ft, -Ztop, Mmin, A, P) 
    c2 = e_inequality(fc, Zbot, Mmin, A, P) 
    t2 = e_inequality(ft, Zbot, Mmax, A, P) 
    
    return c1, t1, c2, t2

def plot_magnel(transfer, service, line_ext=2.0, min_max_points=None):
    """
    Plot the Magnel diagram, loading states: 1) at transfer and 2) in service.

    """

    P0 = transfer.fc * transfer.A/2 # arbitrary, but fc*A/2 suitable
    e0 = calc_inequalities(transfer, P0)
    P1 = service.fc * service.A/2 
    e1 = calc_inequalities(service, P1)

    return _plot_magnel(transfer, P0, e0, service, P1, e1, line_ext=line_ext, min_max_points=min_max_points)

def optimize_magnel(transfer, service, ebounds, mode='min', output=False):
    """
    Solve the linear programming problem, loading states: 1) at transfer and 2) in service.

    """

    A0, fc0, ft0, Ztop0, Zbot0, Mmax0, Mmin0, losses0 = transfer.A, transfer.fc, transfer.ft, transfer.Ztop, transfer.Zbot, transfer.Mmax, transfer.Mmin, transfer.losses
    A1, fc1, ft1, Ztop1, Zbot1, Mmax1, Mmin1, losses1 = service.A, service.fc, service.ft, service.Ztop, service.Zbot, service.Mmax, service.Mmin, service.losses 

    # transfer 
    ct0, ctb0 =  fc0 * -Ztop0 + Mmax0, -Ztop0/A0 # comp at top
    cb0, cbb0 =  fc0 *  Zbot0 + Mmin0,  Zbot0/A0 # comp at bot
    tt0, ttb0 = -ft0 * -Ztop0 + Mmin0, -Ztop0/A0 # tens at top
    tc0, tcb0 = -ft0 *  Zbot0 + Mmax0,  Zbot0/A0 # tens at bot

    # service 
    ct1, ctb1 =  fc1 * -Ztop1 + Mmax1, -Ztop1/A1
    cb1, cbb1 =  fc1 *  Zbot1 + Mmin1,  Zbot1/A1 
    tt1, ttb1 = -ft1 * -Ztop1 + Mmin1, -Ztop1/A1 
    tc1, tcb1 = -ft1 *  Zbot1 + Mmax1,  Zbot1/A1
 
    R = cp.Variable(nonneg=True)
    e = cp.Variable()
    scale = 1e7 # numerical stability 

    if mode == 'max': 
        # maximum P <=> minimum R 
        objective = cp.Minimize(R * scale) 

    if mode == 'min':
        # minimum P <=> maximum R
        objective = cp.Maximize(R * scale)

    # transfer  
    constraints = [
        (ct0 * R*scale)*(1/losses0) - e <= ctb0,   
        (cb0 * R*scale)*(1/losses0) - e >= cbb0, 
        (tt0 * R*scale)*(1/losses0) - e >= ttb0, 
        (tc0 * R*scale)*(1/losses0) - e <= tcb0
    ]

    # service 
    constraints += [
        (ct1 * R*scale)*(1/losses1) - e <= ctb1,   
        (cb1 * R*scale)*(1/losses1) - e >= cbb1, 
        (tt1 * R*scale)*(1/losses1) - e >= ttb1, 
        (tc1 * R*scale)*(1/losses1) - e <= tcb1,
    ]

    # bounds 
    constraints += [
        e >= ebounds[0],
        e <= ebounds[1]
    ]

    problem = cp.Problem(objective, constraints)
    problem.solve()

    if problem.status == cp.OPTIMAL:
        R_opt = R.value*scale
        e_opt = e.value
        if output:
            print(f"Optimal 1/P: {R_opt}")
            print(f"Optimal e  : {e_opt:.0f} mm")
            print(f"Optimal P  : {round(1/R_opt*1e-3):.0f} kN")
        return 1/R_opt, e_opt
    else:
        if output: 
            print(f"Optimization failed. Status: {problem.status}")
        return 0.0, 0.0 
    
def optimize_and_plot_magnel(transfer, service, ebounds, output=False, line_ext=2.0):
    """
    Solve the linear programming problem, loading states: 1) at transfer and 2) in service and plot the associated Magnel diagram. 

    """
        
    min_p, min_e = optimize_magnel(transfer=transfer, service=service, ebounds=ebounds, mode='min', output=output) 
    max_p, max_e = optimize_magnel(transfer=transfer, service=service, ebounds=ebounds, mode='max', output=output) 

    min_max_points = [[1/min_p, min_e], [1/max_p, max_e]]
    
    plot_magnel(transfer=transfer, service=service, line_ext=line_ext, min_max_points=min_max_points) 

