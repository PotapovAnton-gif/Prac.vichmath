
import numpy as np
import matplotlib.pyplot as plt

print(1)
dots = np.array([1, 2, 4, 5])

f_x_es = np.array([1, 3, 1, 3])

def divided_differences(x, f):
    n = len(x);
    F = np.empty((n, n))
    F[:, 0] = f
    for k in range(1, n):
        F[0:n-k, k] = (F[1:n-k+1, k-1] - F[0:n-k, k-1]) / (x[k:] - x[:-k])
    return F # F[i, k] = f(x_i, x_{i+1}, ..., x_{i+k})


F = divided_differences(dots, f_x_es)
print(F)

def evaluate(x, F, x0):
    n = len(x);
    P = 0;
    xprod = 1.0 # (x - x1) (x - x2) ... (x - xi)
    for i in range(n):
        P += F[0, i] * xprod
        xprod *= (x0 - x[i])
    return P

X = np.linspace(0.5, 5.5, 1000)
plt.plot(X, evaluate(dots, F, X), 'r', lw=2, label='$P(x)$')
plt.plot(dots, f_x_es, 'g.', ms=15, label='$f(x_i)$')
plt.legend(loc='center left', bbox_to_anchor=(1, .5)); plt.show()


