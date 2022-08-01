# Вариант 6
import math as m
import numpy as np
import copy


# Правая часть уравнения
def f(x, y):
    return -25 * m.pi * m.pi * m.sin(3 * m.pi * x) * m.sin(4 * m.pi * y)


# Начальное условие
def u0(x, y):
    return 0


# Граничное условие
def phi(x, y):
    return 0


# Аналитическое решение
def anSol(x, y):
    return m.sin(3 * m.pi * x) * m.sin(4 * m.pi * y)


def firstNormMatrix(matrix):
    # суммируем внутри строки
    return max(np.sum(np.abs(matrix), axis=1))


def alternativeNorm(matrix):
    return np.max(np.abs(matrix))


def isEven(x):
    return (x % 2) == 0


coef = 8
M = 41  # число отрезков по икс
x0, xEnd = 0, 1
hX = (xEnd - x0) / M
L = 41  # число отрезков по игрек
y0, yEnd = 0, 1
hY = (yEnd - y0) / L
# eps = 1e-3
eps = 1e-4

# Аналитическое решение
anSolMatrix = np.zeros((6, 6))
for i in range(1, 5):
    for j in range(1, 5):
        anSolMatrix[i][j] = anSol(i * 0.2, j * 0.2)
# anSolMatrix = np.fromfunction(lambda i, j: anSol(i * hX, j * hY), (M + 1, L + 1))

# Начальное условие
uLast = np.zeros((M + 1, L + 1))
for i in range(1, M):
    for j in range(1, L):
        uLast[i][j] = u0(i * hX, j * hY)
# uLast = np.fromfunction(lambda i, j: u0(i * hX, j * hY), (M + 1, L + 1))


muMin = 2 * m.pi * m.pi
muMax = 4 * (L * L + M * M)

# Число шагов по времени:
#N = m.ceil(m.log(2 / eps) / m.log((m.sqrt(muMax) + m.sqrt(muMin)) / (m.sqrt(muMax) - m.sqrt(muMin))))
N = 512
print("Число шагов по времени:", N)

indexesForAlgorithm = [N]  # все индексы
indexesForAlgorithmWithLine = set()  # индексы с подчёркиванием

# Вычисляем индексы
while indexesForAlgorithm[-1] != 1:
    if isEven(indexesForAlgorithm[-1]):
        indexesForAlgorithm.append(int(indexesForAlgorithm[-1] / 2))
    else:
        indexesForAlgorithm.append(indexesForAlgorithm[-1] - 1)

# Вычисялем индексы с подчеркиванием
for i in range(1, len(indexesForAlgorithm) - 1):
    if isEven(indexesForAlgorithm[i]) and not isEven(indexesForAlgorithm[i + 1]) and not isEven(
            indexesForAlgorithm[i - 1]):
        indexesForAlgorithmWithLine.add(indexesForAlgorithm[i])
for i in range(len(indexesForAlgorithm) - 1):
    if not isEven(indexesForAlgorithm[i]):
        indexesForAlgorithmWithLine.add(indexesForAlgorithm[i + 1])
        break

indexesForAlgorithm.reverse()

# Вычисляем тэты
teta = {}
teta[1] = [1]
for i in range(len(indexesForAlgorithm) - 1):
    if 2 * indexesForAlgorithm[i] == indexesForAlgorithm[i + 1]:
        if indexesForAlgorithm[i + 1] not in indexesForAlgorithmWithLine:
            K = int(indexesForAlgorithm[i + 1] / 2)
            array1 = []
            for j in range(K):
                array1.append(teta.get(K)[j])
                array1.append(4 * K - array1[-1])
            teta[indexesForAlgorithm[i + 1]] = array1
        else:
            K = int(indexesForAlgorithm[i + 1] / 2)
            array2 = []
            for j in range(K):
                array2.append(teta.get(K)[j])
                array2.append(4 * K + 2 - array2[-1])
            teta[indexesForAlgorithm[i + 1]] = array2
    else:
        if indexesForAlgorithm[i] not in indexesForAlgorithmWithLine:
            array3 = copy.copy(teta.get(indexesForAlgorithm[i]))
            array3.append(2 * indexesForAlgorithm[i] + 1)
            teta[indexesForAlgorithm[i + 1]] = array3
        else:
            array4 = copy.copy(teta.get(indexesForAlgorithm[i]))
            array4.append(indexesForAlgorithm[i + 1])
            teta[indexesForAlgorithm[i + 1]] = array4

# Очёредность взятия шагов по времени:
tauIndexes = [int((element + 1) / 2) for element in teta.get(indexesForAlgorithm[-1])]
print(teta.get(indexesForAlgorithm[-1]))
print(tauIndexes)
print()

# Вычисляем величину шагов по времени:
tauArray = []
for n in range(1, N + 1):
    tauArray.append(2 / (muMin + muMax + (muMax - muMin) * m.cos(m.pi * (2 * n - 1) / 2 / N)))

# Процесс подсчёта
u = np.zeros((M + 1, L + 1))
for n in range(N):
    tau = tauArray[tauIndexes[n] - 1]
    for i in range(1, M):
        for j in range(1, L):
            u[i][j] = uLast[i][j] + tau / hX / hX * (
                    uLast[i + 1][j] - 2 * uLast[i][j] + uLast[i - 1][j]) + tau / hY / hY * (
                              uLast[i][j + 1] - 2 * uLast[i][j] + uLast[i][j - 1]) - tau * f(i * hX, j * hY)

    uLast = copy.copy(u)

print("Аналитическое решение")
print(anSolMatrix)
print()

m = int(M / 5)
result = np.zeros((6, 6))
for i in range(6):
    string = []
    for j in range(6):
        string.append(uLast[i * m, j * m])
    result[i, :] = np.array(string)
print("Численное решение")
print(result)
print()

print("Разность аналитического и численного решения")
error = np.abs(result - anSolMatrix)
print(error)
print()

print("Число шагов по времени:", N)
print("Число шагов по пространству:", M)
print("Первая норма разности")
print(firstNormMatrix(error))
