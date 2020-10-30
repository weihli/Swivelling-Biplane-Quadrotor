import numpy as np
from RK4_matrix import RK4_matrix as RK4_matrix
from scipy.integrate import quad


def f(t, R):
    return np.array([[1.0, 0.0],[0.0,0.0]])



R0 = np.eye(2)
t0 = 0.0
tf = 10.0
h = 0.0001
n = 2
m = 2

print(RK4_matrix(f,R0, t0, tf, h, n, m)[0])

# print(0.01*np.eye(3))
