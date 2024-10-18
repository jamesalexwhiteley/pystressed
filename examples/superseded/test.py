import numpy as np
from scipy.optimize import linprog

np.random.seed(1)

# min c^Tx s.t. Ax<b 
A = np.random.randn(8, 2) 
b = np.random.randn(8)    
c = np.array([0, 1])

emin, emax = 2, 3
bounds = [(None, None), (emin, emax)] 
# bounds = [(None, None), (None, None)]
res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method='highs')

if res.success:
    print(f"Optimal solution: x = {res.x}")
    print(f"Objective value: {res.fun}")
else:
    print("Optimization failed:", res.message)

