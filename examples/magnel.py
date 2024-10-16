import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linprog
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# Author: James Whiteley (github.com/jamesalexwhiteley)

class SectionForce():
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

def plot_line(start, end, label, alpha, linestyle, extension):
    delta_x = end[0] - start[0]
    delta_y = end[1] - start[1]
    extended_end = (start[0] + delta_x * extension, start[1] + delta_y * extension)
    plt.plot([start[0], extended_end[0]], [start[1], extended_end[1]], 'k', label=label, alpha=alpha, linestyle=linestyle)
    # plt.plot([start[1], extended_end[1]], [start[0], extended_end[0]], 'k', label=label, alpha=alpha, linestyle=linestyle)
    return (start, extended_end)

def plot_lines(Z, A, P, e, labels, alphas, linestyles, line_ext):
    for i, z in enumerate(Z):
        start, end = plot_line((0, -z / A), (1 / P, e[i]), labels[i], alphas[i], linestyles[i], line_ext)

def _plot_magnel(transfer, P0, e0, service, P1, e1, line_ext=2.0):
    
    # plot Magnel diagram      
    Ztop0, Zbot0, A0 = transfer.Ztop, transfer.Zbot, transfer.A
    Z0 = [-Ztop0, -Ztop0, Zbot0, Zbot0]  

    Ztop1, Zbot1, A1 = service.Ztop, service.Zbot, service.A
    Z1 = [-Ztop1, -Ztop1, Zbot1, Zbot1]    
    
    labels = ['c1', 't1', 'c2', 't2'] 
    alphas = [1, .4, 1, .4]
    linestyles = ['-', '-', '--', '--'] 

    plt.figure()

    plot_lines(Z0, A0, P0, e0, labels, alphas, linestyles, line_ext)
    # plot_lines(Z1, A1, P1, e1, labels, alphas, linestyles, line_ext)
    
    plt.xlabel("1/P")
    plt.ylabel("e")
    plt.gca().invert_yaxis()
    plt.grid(True)
    plt.legend()
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

def plot_magnel(transfer, service, ebounds, line_ext=2.0):

    P0 = transfer.fc * transfer.A/2 # arbitrary, but fc*A/2 suitable
    e0 = calc_inequalities(transfer, P0)
    P1 = service.fc * service.A/2 
    e1 = calc_inequalities(service, P1)

    return _plot_magnel(transfer, P0, e0, service, P1, e1, line_ext=line_ext)

# def linear_program(transfer, service, ebounds):
#     fc0, ft0, A0, Ztop, Zbot, Mmax0, Mmin0, alpha0 = transfer.fc, transfer.ft, transfer.A, transfer.Ztop, transfer.Zbot, transfer.Mmax, transfer.Mmin, transfer.alpha
#     fc1, ft1, A1, Ztop, Zbot, Mmax1, Mmin1, alpha1 = service.fc, service.ft, service.A, service.Ztop, service.Zbot, service.Mmax, service.Mmin, transfer.alpha

#     # min c^Tx s.t. Ax<b  
#     A = np.array(
#         [[ -(Mmax0 + Ztop*fc0)*(1/alpha0), 1 ],
#          [ -(Mmax1 + Ztop*fc1)*(1/alpha1), 1 ],
#          [ -(Mmin0 + Ztop*ft0)*(1/alpha0), 1 ],
#          [ -(Mmin1 + Ztop*ft1)*(1/alpha1), 1 ],
#          [ -(Mmin0 + Zbot*fc0)*(1/alpha0), 1 ],
#          [ -(Mmin1 + Zbot*fc1)*(1/alpha1), 1 ],
#          [ -(Mmax0 + Zbot*ft0)*(1/alpha0), 1 ],
#          [ -(Mmax1 + Zbot*ft1)*(1/alpha1), 1 ]])

#     b = np.array([-Ztop/A0, -Ztop/A1, -Ztop/A0, -Ztop/A1, Zbot/A0, Zbot/A1, Zbot/A0, Zbot/A1])  
#     c = np.array([0, 1])

#     # bounds = [(None, None), (ebounds[0], ebounds[1])] 
#     bounds = [(None, None), (None, None)] 
#     res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method='highs') 

#     if res.success:
#         print(f"Optimal solution: x = {res.x}")
#         print(f"Objective value: {res.fun}")
#     else:
#         print("Optimization failed:", res.message)

def _linear_program(transfer, service, ebounds):
    fc, ft, A0, Ztop, Zbot, Mmax, Mmin, _ = transfer.fc, transfer.ft, transfer.A, transfer.Ztop, transfer.Zbot, transfer.Mmax, transfer.Mmin, transfer.alpha

    # min c^Tx s.t. Ax<b  
    A = np.array(
        [[ -(Mmax + Ztop*fc), 1 ],
         [ -(Mmin + Ztop*ft), 1 ],
         [ -(Mmin + Zbot*fc), 1 ],
         [ -(Mmax + Zbot*ft), 1 ]])
    
    b = np.array([-Ztop/A0, -Ztop/A0, Zbot/A0, Zbot/A0])  
    c = np.array([0, 1]) 

    bounds = [(0, None), (0, None)]   
    # bounds = [(0, None), (ebounds[0], ebounds[1])] # 1/P bounds, e bounds 
    res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method='highs')     

    if res.success:
        print(f"Optimal solution: x = {res.x}")
        print(f"Objective value: {res.fun}")
    else:
        print("Optimization failed:", res.message)

# # transfer   
# A = 1.04e6 # mm2
# fc = 20 # N/mm2 
# ft = -1 # N/mm2
# Ztop = 225e6 # mm3  
# Zbot = 225e6 # mm3
# Mmax = 1500e6 # Nmm
# Mmin = 1300e6 # Nmm
# alpha = 0.95
# transfer = SectionForce(A=A, fc=fc, ft=ft, Ztop=Ztop, Zbot=Zbot, Mmax=Mmax, Mmin=Mmin, alpha=alpha)

# service 
A = 1.04e6 # mm2
fc = 30 # N/mm2 
ft = 0 # N/mm2
Ztop = 225e6 # mm3  
Zbot = 225e6 # mm3
Mmax = 2300e6 # Nmm
Mmin = 1500e6 # Nmm
alpha = 0.85
service = SectionForce(A=A, fc=fc, ft=ft, Ztop=Ztop, Zbot=Zbot, Mmax=Mmax, Mmin=Mmin, alpha=alpha)

# plot_magnel(transfer=transfer, service=service, ebounds=[0, 550], line_ext=1.5)

# 4d7
P = 1280e3 # N 
A = 0.16e6 # mm2
fc = 16 # N/mm2 
ft = 0 # N/mm2
Ztop = 34.4e6 # mm3
Zbot = 29.3e6 # mm3
Mmax = 530e6 # Nmm
Mmin = 108e6 # Nmm  
alpha = 1.0

transfer = SectionForce(A=A, fc=fc, ft=ft, Ztop=Ztop, Zbot=Zbot, Mmax=Mmax, Mmin=Mmin, alpha=alpha)

_linear_program(transfer=transfer, service=service, ebounds=[0, 550])
# plot_magnel(transfer=transfer, service=service, ebounds=[0, 550], line_ext=1.5)