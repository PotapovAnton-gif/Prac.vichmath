# Работа 6
# Вариант 1
# Задание 9

from math import *
from random import randint
from matplotlib import pyplot as plt
import numpy as np
import copy as copy
from scipy import rand

nu = 0
mu = 1/3
C_0 = 15
С_1 = 1
eps = 1e-4
T = 1


def anal_sol(t, r):
    return ((С_1 + r) ** (2 / mu)) / (C_0 - 2 * (mu + 2) * t / mu) ** (1 / mu)


def solve_triagonal_slae(a:float, b:float, c:float, d:float) -> float:
    # a - нижняя диагональ
    # b - диагональ
    # c - верхняя диагональ
    N = len(a) - 1

    # Вычисляем p1 и q1
    p = [0] + [-c[0] / b[0]]
    q = [0] + [d[0] / b[0]]

    for i in range(1, N):
        znam = a[i] * p[i] + b[i]
        p.append(-c[i] / znam)
        q.append((d[i] - a[i] * q[i]) / znam)

    # Обратный ход:
    yN = (d[N] - a[N] * q[N]) / (a[N] * p[N] + b[N])
    result = [0 for _ in range(N)] + [yN]
    for i in range(N - 1, -1, -1):
        result[i] = p[i + 1] * result[i + 1] + q[i + 1]
    return result


def num_sol(L, N):
    t_0, t_end = 0, T
    r_0, r_end = 0, 1

    # N = 10  # число отрезков по t
    tau = (t_end - t_0) / N
    # L = 10  # число отрезков по r
    h = (r_end - r_0) / L

    # Начальное условие
    def u_0(x):
        return ((1 + x)**6)/3375

    # Левое граничное условие
    def u_left(t):
        return 1/((15-14*t)**3)

    # Правое граничное условие
    def u_right(t):
        return 64/((15-14*t)**3)

    def r_func(l):
        return r_0 + (r_end - r_0) * l / L

    result = [u_0(r_func(l)) for l in range(L + 1)]
    for n in range(N):
        u_last = copy.copy(result)

        index = 0
        while (True):
            A = [0] + [(u_last[l] ** mu + u_last[l - 1] ** mu)*tau / (2 * h * h)
                       for l in range(1, L)]
            B = [0] + [(u_last[l+1] ** mu + u_last[l] ** mu)*tau / (2 * h * h)
                       for l in range(1, L)]

            a = [0] + [-B[l] for l in range(1, L)] + [0]
            b = [1] + [1 + A[l] + B[l] for l in range(1, L)] + [1]
            c = [0] + [-A[l] for l in range(1, L)]
            d = [u_left(tau * (n + 1))] + [result[l] for l in range(1, L)] + [u_right(tau * (n + 1))]

            u_cur = solve_triagonal_slae(a, b, c, d)

            error = max(
                [(u_cur[l] - u_last[l]) / u_cur[l] if u_cur[l] != 0 else (u_cur[l] - u_last[l]) for l in range(L + 1)])

            u_last = copy.copy(u_cur)
            index += 1

            if error <= eps:
                break
        # print(index)

        result = copy.copy(u_last)

    return result
    #2
#0,69
#0,3






# L - число точек по x
# N - число точек по t
L = 21
N = 121



n_1 = 1
# N = ceil(int(L / 10) * 10 / CFL + 1)
n = int(L / 10)
X = [i / 10 for i in range(11)]
An_sol = [anal_sol(T, x) for x in X]
Num_sol =  num_sol(L - 1, N - 1)[::n]





Diff = [abs((An_sol[i] - Num_sol[i])/10)/n_1 for i in range(11)]
MaxDiff = [max(Diff)]
print()
# print(X)
# print(An_sol)
# print(Num_sol_const)
# print(Diff)
# print(MaxDiff)
# print(Num_sol_variable)

tabledata = [["x"] + X, ["Аналит. решение"] + An_sol, ["Числен. решение"] + Num_sol, ["Модуль разности"] + Diff,
             ["Мах модуля разности"] + MaxDiff]
from tabulate import tabulate

print(tabulate(tabledata))
print("Число точек по r: ", L)
print("Число точек по t: ", N)
#