import numpy as np
from scipy.sparse import dia_matrix
from numpy import linalg as LA


def q(x):
    if 0.5 >= x >= 0:
        return 2
    elif 0.5 <= x <= 1:
        return 1


def alpha(x):
    return 1


x0, xN = 0, 1
h = 0.0002
L = int(xN / h)  # число отрезков
x_arr = [x0 + i * (xN - x0) / L for i in range(L + 1)]
beta = [q(x) for x in x_arr]
#alpha = [1 for x in x_arr]
#alpha_minus_one = 2 - (-h)
#alpha = 1


a_arr = [1 for i in range(L - 1)] + [0] + [1000]
c_arr = [1000] + [2] + [1 for i in range(2, L + 1)]


def calcDet(lam):
    b_arr = [lam * beta[i] * h * h - 2 for i in range(L)] + [1]
    return LA.det(dia_matrix((np.array([b_arr, a_arr, c_arr]), [0, -1, 1]), shape=(L + 1, L + 1)).toarray())


# Локализуем корни:
n = 10
z1, z2 = 0, 250
value = [calcDet(z1 + i * (z2 - z1) / n) for i in range(n + 1)]
roots = []
for i in range(n):
    if value[i] * value[i + 1] < 0:
        x1 = z1 + i * (z2 - z1) / n
        x2 = z1 + (i + 1) * (z2 - z1) / n
        roots.append((x1, x2))
print('\nЛокализуем корни:')
for i in range(len(roots)):
    print(roots[i])

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
