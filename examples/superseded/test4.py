import cvxpy as cp
import numpy as np

def optimize(A, fc, ft, Ztop, Zbot, Mmax, Mmin, ebounds, mode='min'):

    # transfer 
    ct, ctb =  fc * -Ztop + Mmax, -Ztop/A # comp@top
    cb, cbb =  fc *  Zbot + Mmin,  Zbot/A # comp@bot
    tt, ttb = -ft * -Ztop + Mmin, -Ztop/A # tens@top
    tc, tcb = -ft *  Zbot + Mmax,  Zbot/A # tens@bot

    print(cb)
    print(cbb)

    R = cp.Variable(nonneg=True)
    e = cp.Variable()

    scale = 1e7

    if mode == 'max': 
        # maximum P <=> minimum R 
        objective = cp.Minimize(R * scale) 

    if mode == 'min':
        # minimum P <=> maximum R
        objective = cp.Maximize(R * scale)

    constraints = [
        (ct * R*scale) - e <= ctb,   
        (cb * R*scale) - e >= cbb, 
        (tt * R*scale) - e >= ttb, 
        (tc * R*scale) - e <= tcb
    ]

    constraints += [
        e >= ebounds[0],
        e <= ebounds[1]
    ]

    # if mode == 'max':
    #     # maximum P <=> minimum R
    #     objective = cp.Minimize(R) 
    #     constraints = [
    #         (ct * R) - e <= ctb,   
    #         (cb * R) - e >= cbb, 
    #         (tt * R) - e >= ttb, 
    #         (tc * R) - e <= tcb, 
    #         # e >= ebounds[0],
    #         # e <= ebounds[1]
    #     ] 

    # if mode == 'min':
    #     # minimum P <=> maximum R
    #     objective = cp.Maximize(R)
    #     constraints = [
    #         (ct * R) - e <= ctb,   
    #         (cb * R) - e >= cbb, 
    #         (tt * R) - e >= ttb, 
    #         (tc * R) - e <= tcb, 
    #         # e >= ebounds[0], 
    #         # e <= ebounds[1] 
    #     ] 

    problem = cp.Problem(objective, constraints)
    problem.solve()

    if problem.status == cp.OPTIMAL:
        R_opt = R.value*scale
        e_opt = e.value
        print(f"Optimal R: {R_opt}")
        print(f"Optimal e: {e_opt}")
        print(f"Optimal P: {1/R_opt*1e-3}")
    else:
        print(f"Optimization failed. Status: {problem.status}")

# Q1
A = 0.64e6 # mm2
fc = 15 # N/mm2 
ft = 1 # N/mm2
Ztop = 309e6 # mm3
Zbot = 265e6 # mm3
Mmax = 5000e6 # Nmm
Mmin = 1000e6 # Nmm
ebounds = [0, 1000]

# # Q2
# # parameters 
# A = 0.16e6 # mm2
# fc = 16 # N/mm2 
# ft = 0 # N/mm2
# Ztop = 34.4e6 # mm3
# Zbot = 29.3e6 # mm3
# Mmax = 530e6 # Nmm
# Mmin = 108e6 # Nmm
# ebounds = [0, np.inf]

optimize(A, fc, ft, Ztop, Zbot, Mmax, Mmin, ebounds, mode='min')