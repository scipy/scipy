"""
Solving a Transportation Problem with linprog
=============================================

This example shows how to solve a classic transportation problem
using `scipy.optimize.linprog`.

A company has 3 warehouses and 4 stores. Each warehouse has supply,
each store has demand. We minimize total shipping cost.

Supply: W0=20, W1=30, W2=25
Demand: S0=15, S1=25, S2=20, S3=15
Costs:
  [[8,  6, 10, 9],
   [9, 12, 13, 7],
   [14, 9, 16, 5]]
"""

import numpy as np
from scipy.optimize import linprog

supply = [20, 30, 25]
demand = [15, 25, 20, 15]
costs = [[8, 6, 10, 9],
         [9, 12, 13, 7],
         [14, 9, 16, 5]]

nw = len(supply)
ns = len(demand)
c = np.array(costs).flatten()

# Supply constraints
A_ub = np.zeros((nw, nw * ns))
for i in range(nw):
    A_ub[i, i*ns:(i+1)*ns] = 1
b_ub = supply

# Demand constraints
A_eq = np.zeros((ns, nw * ns))
for j in range(ns):
    A_eq[j, j::ns] = 1
b_eq = demand

bounds = [(0, None)] * len(c)
res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
              bounds=bounds, method='highs')

print("Solution found:", res.success)
print("Total cost:", res.fun)
x = res.x.reshape((nw, ns))
print("\nShipping plan:")
for i in range(nw):
    for j in range(ns):
        if x[i, j] > 1e-6:
            print(f"  W{i} -> S{j}: {x[i, j]:.1f}")
