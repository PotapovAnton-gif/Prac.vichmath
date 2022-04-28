# Задание 4

from math import *

# Аналитическое решение
def an_sol(x):
    x0 = 0.525
    ka = x0 + 1
    kb = x
    qa = x0 * x0
    qb = exp(-x0 * x0)
    fa = 1
    fb = cos(x0)
    la = sqrt(qa / ka)
    lb = sqrt(qb / kb)
    ma = fa / qa
    mb = fb / qb
    u0 = 1
    u1 = 0
    u = 0

    a11 = exp(-la * x0) - exp(la * x0)
    a12 = exp(lb * (2 - x0)) - exp(lb * x0)
    a21 = ka * la * (exp(la * x0) + exp(-la * x0))
    a22 = kb * lb * (exp(lb * (2 - x0)) + exp(lb * x0))
    b1 = mb - ma + (ma - u0) * exp(la * x0) - (mb - u1) * exp(lb * (1 - x0))
    b2 = ka * la * (u0 - ma) * exp(la * x0) + kb * lb * (u1 - mb) * exp(lb * (1 - x0))
    c1 = (((u0 - ma) * a11 - b1) * a22 - ((u0 - ma) * a21 - b2) * a12) / (a11 * a22 - a12 * a21)
    c2 = (b1 * a22 - b2 * a12) / (a11 * a22 - a12 * a21)
    c3 = (b2 * a11 - b1 * a21) / (a11 * a22 - a12 * a21)
    c4 = (u1 - mb) * exp(lb) - c3 * exp(2 * lb)

    if x <= x0:
        u = c1 * exp(la * x) + c2 * exp(-la * x) + ma
    else:
        u = c3 * exp(lb * x) + c4 * exp(-lb * x) + mb
    return u

# Численное решение модельной задачи
def num_sol_const(N):
    x0 = 1 / sqrt(2)
    ka = exp(-x0)
    kb = 1
    qa = x0 * x0
    qb = exp(-x0 * x0)
    fa = 1
    fb = cos(x0)
    u0 = 1
    u1 = 0

    L = N - 1
    h = 1 / L
    ia = floor(x0 / h)
    ib = ia + 1

    al = [0] + [ka for i in range(1, ia)] + [0, 0] + [kb for i in range(ib + 1, L)] + [0]
    # print(len(al))
    bl = [0] + [-2 * ka - qa * h * h for i in range(1, ia)] + [0, 0] + [-2 * kb - qb * h * h for i in
                                                                        range(ib + 1, L)] + [0]
    cl = [0] + [ka for i in range(1, ia)] + [0, 0] + [kb for i in range(ib + 1, L)] + [0]
    dl = [0] + [-fa * h * h for i in range(1, ia)] + [0, 0] + [-fb * h * h for i in range(ib + 1, L)] + [0]

    alpl = [0, -al[1] / bl[1]]
    for i in range(2, ia):
        alpl.append(-al[i] / (bl[i] + cl[i] * alpl[i - 1]))
    alpl.append(0)
    alpl.append(0)
    temp1 = [-cl[L - 1] / bl[L - 1]]
    for i in range(L - 2, ib, -1):
        temp1.append(-cl[i] / (bl[i] + al[i] * temp1[-1]))
    alpl = alpl + temp1[::-1]
    alpl.append(0)

    betl = [0, (dl[1] - cl[1] * u0) / bl[1]]
    for i in range(2, ia):
        betl.append((dl[i] - cl[i] * betl[i - 1]) / (bl[i] + cl[i] * alpl[i - 1]))
    betl.append(0)
    betl.append(0)
    temp2 = [(dl[L - 1] - cl[L - 1] * u1) / bl[L - 1]]
    for i in range(L - 2, ib, -1):
        temp2.append((dl[i] - al[i] * temp2[-1]) / (bl[i] + al[i] * alpl[i + 1]))
    betl = betl + temp2[::-1]
    betl.append(0)

    u = [u0] + [0 for i in range(1, N - 1)] + [u1]
    u[ia] = (ka * betl[ia - 1] + kb * betl[ib + 1]) / (ka * (1 - alpl[ia - 1]) + kb * (1 - alpl[ib + 1]))
    u[ib] = u[ia]
    # print(ia)
    u[ia - 1] = alpl[ia - 1] * u[ia] + betl[ia - 1]
    u[ib + 1] = alpl[ib + 1] * u[ib] + betl[ib + 1]

    for i in range(ia - 2, 0, -1):
        u[i] = alpl[i] * u[i + 1] + betl[i]
    for i in range(ib + 2, L):
        u[i] = alpl[i] * u[i - 1] + betl[i]
    return u

# Численное решение задачи с перем. к-тами
def num_sol_variable(N):
    x0 = 1 / sqrt(2)

    def ka(x):
        return exp(-x)

    def kb(x):
        return 1

    def qa(x):
        return x * x

    def qb(x):
        return exp(-x * x)

    def fa(x):
        return 1

    def fb(x):
        return cos(x)

    u0 = 1
    u1 = 0

    L = N - 1
    h = 1 / L
    ia = floor(x0 / h)
    ib = ia + 1

    al = [0] + [ka(i * h + h / 2) for i in range(1, ia)] + [0, 0] + [kb(i * h + h / 2) for i in range(ib + 1, L)] + [0]

    bl = [0] + [-(ka(i * h + h / 2) + ka(i * h - h / 2) + qa(i * h) * h * h) for i in range(1, ia)] + [0, 0] + [
        -(kb(i * h + h / 2) + kb(i * h - h / 2) + qb(i * h) * h * h) for i in range(ib + 1, L)] + [0]

    cl = [0] + [ka(i * h - h / 2) for i in range(1, ia)] + [0, 0] + [kb(i * h - h / 2) for i in range(ib + 1, L)] + [0]

    dl = [0] + [-fa(i * h) * h * h for i in range(1, ia)] + [0, 0] + [-fb(i * h) * h * h for i in range(ib + 1, L)] + [0]

    alpl = [0, -al[1] / bl[1]]
    for i in range(2, ia):
        alpl.append(-al[i] / (bl[i] + cl[i] * alpl[i - 1]))
    alpl.append(0)
    alpl.append(0)
    temp1 = [-cl[L - 1] / bl[L - 1]]
    for i in range(L - 2, ib, -1):
        temp1.append(-cl[i] / (bl[i] + al[i] * temp1[-1]))
    alpl = alpl + temp1[::-1]
    alpl.append(0)

    betl = [0, (dl[1] - cl[1] * u0) / bl[1]]
    for i in range(2, ia):
        betl.append((dl[i] - cl[i] * betl[i - 1]) / (bl[i] + cl[i] * alpl[i - 1]))
    betl.append(0)
    betl.append(0)
    temp2 = [(dl[L - 1] - cl[L - 1] * u1) / bl[L - 1]]
    for i in range(L - 2, ib, -1):
        temp2.append((dl[i] - al[i] * temp2[-1]) / (bl[i] + al[i] * alpl[i + 1]))
    betl = betl + temp2[::-1]
    betl.append(0)

    u = [u0] + [0 for i in range(1, N - 1)] + [u1]
    u[ia] = (ka(h * ia) * betl[ia - 1] + kb(h * ia) * betl[ib + 1]) / (
            ka(h * ia) * (1 - alpl[ia - 1]) + kb(h * ia) * (1 - alpl[ib + 1]))
    u[ib] = u[ia]
    u[ia - 1] = alpl[ia - 1] * u[ia] + betl[ia - 1]
    u[ib + 1] = alpl[ib + 1] * u[ib] + betl[ib + 1]

    for i in range(ia - 2, 0, -1):
        u[i] = alpl[i] * u[i + 1] + betl[i]
    for i in range(ib + 2, L):
        u[i] = alpl[i] * u[i - 1] + betl[i]
    return u


N = 81921
n = int(N / 10)
X = [i / 10 for i in range(11)]
An_sol = [an_sol(x) for x in X]
Num_sol_const = num_sol_const(N)[::n]
Diff = [abs(An_sol[i] - Num_sol_const[i]) for i in range(11)]
MaxDiff = [max(Diff)]
Num_sol_variable = num_sol_variable(N)[::n]

# print(X)
# print(An_sol)
# print(Num_sol_const)
# print(Diff)
# print(MaxDiff)
# print(Num_sol_variable)

tabledata = [["Икс"] + X, ["Аналит. решение"] + An_sol, ["Числен. решение"] + Num_sol_const, ["Модуль разности"] + Diff,
             ["Мах модуля разности"] + MaxDiff, ["Реш. с перем. к-тами"] + Num_sol_variable]
from tabulate import tabulate

print(tabulate(tabledata))
print("Число точек: ", N)
