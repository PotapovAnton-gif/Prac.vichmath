# Работа 5
# Вариант 3
# Задание 3

from math import *


# Аналитическое решение
def an_sol(x, t):
    return (x + t)**2 + cos(x)


# Численное решение
def num_sol(L, CFL):
    h = 1 / (L - 1)
    tau = CFL * h
    t = 0
    x = [h * i for i in range(L)]

    u = [x[i] for i in range(L)]
    u_next = [0 for i in range(L)]

    eps = 1e-10
    while t < 1 - eps:
        if abs(t + tau - 1) >= eps and t + tau > 1:
            tau = 1 - t
        u0 = -(1 + t + tau)**2 + cos(1)
        u_next[L - 1] = u0
        u_next[L - 2] = u0 - h
        u_next[L - 3] = u0 - 2 * h
        for l in range(L - 3):
            u_next[l] = u[l] + tau/(6*h) * (2*u[l+3] - 9*u[l+2] + 18*u[l+1] - 11*u[l]) + tau**2/(2*h**2) * (
                                -u[l+3] + 4*u[l+2] - 5*u[l+1] + 2*u[l]) + tau*sin(h*l) + 1/2*tau**2*cos(h*l) - tau**3/(
                                    6*h**3) * (u[l+3] - 3*u[l+2] + 3*u[l+1] - u[l]) - tau**3/6*sin(h*l)
        t += tau
        #print(t)
        for i in range(L):
            u[i] = u_next[i]
    return u


# L - число точек по x
# N - число точек по t
L = 11
CFL = 0.25
N = ceil(int(L / 10) * 10 / CFL + 1)
n = int(L / 10)
X = [i / 10 for i in range(11)]
An_sol = [an_sol(x, 1) for x in X]
Num_sol = num_sol(L, CFL)[::n]
Diff = [abs(An_sol[i] - Num_sol[i]) for i in range(11)]
MaxDiff = [max(Diff)]

# print(X)
# print(An_sol)
# print(Num_sol_const)
# print(Diff)
# print(MaxDiff)
# print(Num_sol_variable)

tabledata = [["Икс"] + X, ["Аналит. решение"] + An_sol, ["Числен. решение"] + Num_sol, ["Модуль разности"] + Diff,
             ["Мах модуля разности"] + MaxDiff]
from tabulate import tabulate

print(tabulate(tabledata))
print("Число точек по x: ", L)
print("Число точек по t: ", N)
