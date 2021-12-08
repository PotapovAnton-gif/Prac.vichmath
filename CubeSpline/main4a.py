import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def return_cofficient():    # Вариант 1
    # Начальные данные:
    x = [0.17453, 0.5236, 0.87267, 1.22173, 1.5708, 1.91986, 2.26893]
    y = [0.00162, 0.00252, 0.00498, 0.0129, 0.03964, 0.1207, 0.34097]

    #x = [-2, -1, 0, 1, 2]
    #y = [17, 5, 5, 5, 17]

    N = len(y) - 1  # Порядок полинома

    # Вычисляю интерполянт методом Ньютона:
    t = x
    difference = [y]
    for i in range(N):
        temp = []
        for j in range(N - i):
            temp.append((difference[i][j] - difference[i][j + 1]) / (t[j] - t[j + i + 1]))
        difference.append(temp)
    odds = [difference[i][0] for i in range(N + 1)]  # массив искомых коэффициентов
    # Коэффициенты при степенях x:
    a0 = odds[0] + odds[1] * (-t[0]) + odds[2] * (t[0] * t[1]) + odds[3] * (-t[0] * t[1] * t[2]) + odds[4] * (
            t[0] * t[1] * t[2] * t[3])
    a1 = odds[1] + odds[2] * (-t[1] - t[0]) + odds[3] * (t[0] * t[1] + t[1] * t[2] + t[0] * t[2]) + odds[4] * (
            -t[0] * t[1] * t[2] - t[0] * t[1] * t[3] - t[1] * t[2] * t[3] - t[0] * t[2] * t[3])
    a2 = odds[2] + odds[3] * (-t[1] - t[0] - t[2]) + odds[4] * (
            t[0] * t[1] + t[1] * t[2] + t[0] * t[2] + t[1] * t[3] + t[0] * t[3] + t[2] * t[3])
    a3 = odds[3] + odds[4] * (-t[1] - t[0] - t[2] - t[3])
    a4 = odds[4]
    print(a0, a1, a2, a3, a4)
    z1 = 2
    print(a0 + a1 * z1 + a2 * z1 * z1 + a3 * z1 ** 3 + a4 * z1 ** 4)

    # Строю график интерполянта:
    fig = plt.figure()
    x_ = np.arange(0.5, 2.5, 0.05)
    Interpolant = np.zeros(40)
    for j in range(N + 1):
        p = np.ones(40)
        for i in range(j):
            p = p * (x_ - x[i])
        Interpolant = Interpolant + odds[j] * p
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(x,y,'bo')
    ax.scatter(x, y)
    ax.plot(x_, Interpolant)
    ax.set_xticks(x)
    ax.grid()
    plt.show()

    P = [a1 + 2 * a2 * u + 3 * a3 * u * u + 4 * a4 * u * u * u for u in x]
    a = [[], [], [], [], [], []]
    for i in range(N):
        a[i].append(
            (-1 * P[i + 1] * x[i] * x[i] * x[i + 1] * (x[i + 1] - x[i]) + y[i + 1] * x[i] * x[i] * (
                    3 * x[i + 1] - x[i])) / (
                    (x[i + 1] - x[i]) ** 3) + (
                    y[i] * x[i + 1] * x[i + 1] * (x[i + 1] - 3 * x[i]) - P[i] * x[i] * x[i + 1] * x[i + 1] * (
                    x[i + 1] - x[i])) / ((x[i + 1] - x[i]) ** 3))
        a[i].append(
            (P[i + 1] * x[i] * (2 * x[i + 1] + x[i]) * (x[i + 1] - x[i]) - 6 * (y[i + 1] - y[i]) * x[i] * x[i + 1]) / (
                    (x[i + 1] - x[i]) ** 3) + (P[i] * x[i + 1] * (x[i + 1] + 2 * x[i]) * (x[i + 1] - x[i])) / (
                    (x[i + 1] - x[i]) ** 3))
        a[i].append(
            (-1 * P[i + 1] * (x[i + 1] - x[i]) * (x[i + 1] + 2 * x[i]) + 3 * (y[i + 1] - y[i]) * (x[i + 1] + x[i])) / (
                    (x[i + 1] - x[i]) ** 3) - (P[i] * (x[i + 1] - x[i]) * (x[i] + 2 * x[i + 1])) / ((x[i + 1] - x[i]) ** 3))
        a[i].append(
            (P[i + 1] * (x[i + 1] - x[i]) - 2 * (y[i + 1] - y[i]) + P[i] * (x[i + 1] - x[i])) / ((x[i + 1] - x[i]) ** 3))

    
    data = pd.DataFrame(a)
    return data
    z2 = 2
    k = 0
    '''
    while x[k + 1] < z2:
        k = k + 1
    print(a[k][0] + a[k][1] * z2 + a[k][2] * z2 * z2 + a[k][3] * z2 ** 3)
    '''
print(return_cofficient())