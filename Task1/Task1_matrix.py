import math
import numpy as np
from numpy import linalg

N = 50
a = math.cos(2)
b = 5
c = 3.141592653589793238462643383279

matrix = np.zeros((N + 1, N + 1))
matrix[0][0] = 1
matrix[1][1] = 1
for i in range(2, N + 1):
    matrix[i][i - 2] = c
    matrix[i][i - 1] = b
    matrix[i][i] = a

reversed_matrix = linalg.inv(matrix)
print(max(np.sum(np.absolute(reversed_matrix), axis=1)))



# расчет задачи численно
def y_next(y0, y1):
    return (2.718281828459045235360287471352 - 5 * y1 - 3.141592653589793238462643383279 * y0) / math.cos(2)


y0 = 2.718281828459045235360287471352 / (math.cos(2) + 5 + 3.141592653589793238462643383279)
y1 = y0
N = 50
i = 1
while N > 0:
    a = y1
    y1 = y_next(y0, y1)
    y0 = a
    N -= 1
    print(i, ": ", y1)
    i += 1
print()
# расчет задачи аналитически
lambda_1 = (-5 + (25 - 4 * 3.141592653589793238462643383279 * math.cos(2)) ** 0.5) / (2 * math.cos(2))
lambda_2 = (-5 - (25 - 4 * 3.141592653589793238462643383279 * math.cos(2)) ** 0.5) / (2 * math.cos(2))
c_2 = (1 - lambda_1) / (lambda_2 - lambda_1)
c_1 = (1 - lambda_2) / (lambda_1 - lambda_2)

y_error = 10 ** -6
y_10 = c_1 * (lambda_1 ** 10) * y_error + c_2 * (lambda_2 ** 10) * y_error + \
       (2.718281828459045235360287471352 / (math.cos(2)) + 5 + 3.141592653589793238462643383279)
print("lambda_1:", lambda_1)
print("lambda_2:", lambda_2)
print("const_1:", c_1)
print("const_2:", c_2)
print("Аналитическая ошибка, при n=10, y_error=10^(-6): ", y_10)
