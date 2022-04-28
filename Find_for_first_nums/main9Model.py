import numpy as np
from scipy.sparse import dia_matrix
from numpy import linalg as LA

beta = 0.25
alpha = 1
x0, xN = 0, 1
h = 0.001
L = int(xN / h)  # число отрезков

a_arr = [1 for i in range(1, L)] + [0] + [1000]
c_arr = [1000] + [2] + [1 for i in range(1, L)]


def calcDet(lam):
    b_arr = [lam * beta * h ** 2 / alpha - 2 for i in range(L)] + [1]
    return LA.det(dia_matrix((np.array([b_arr, a_arr, c_arr]), [0, -1, 1]), shape=(L + 1, L + 1)).toarray())


# Локализуем корни:
n = 10
z1, z2 = 0, 500
value = [calcDet(z1 + i * (z2 - z1) / n) for i in range(n + 1)]
roots = []
for i in range(n):
    if value[i] * value[i + 1] < 0:
        x1 = z1 + i * (z2 - z1) / n
        x2 = z1 + (i + 1) * (z2 - z1) / n
        roots.append((x1, x2))
print('\nЛокализуем корни:')
for i in range(len(roots)): print(roots[i])

# Уточняем значение корня методом половинного деления:
eps = 1e-4
Y = []
for i in range(len(roots)):
    a = roots[i][0]
    b = roots[i][1]
    while (b - a > eps):
        c = (a + b) / 2
        if calcDet(c) * calcDet(a) < 0:
            b = c
        else:
            a = c
    Y.append((a + b) / 2)
print('\nУточняем значение корней:')
for i in range(len(Y)): print(Y[i])

# Локализуем корни:
# (0.0, 50.0)
# (50.0, 100.0)
# (200.0, 250.0)
# (450.0, 500.0)
# Уточняем значение корней:
# 9.869623184204102
# 88.82641792297363
# 246.74010276794434
# 483.6103916168213