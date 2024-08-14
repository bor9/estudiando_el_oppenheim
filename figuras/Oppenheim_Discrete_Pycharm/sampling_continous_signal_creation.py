import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

__author__ = 'ernesto'

# se crea una señal suave definiendo algunos puntos xn en los instantes n
# y luego empleando interpolación cúbica.

xn = [3.5, 3, 2.2, 2.3, 2, 4, 5, 3.5, 2, 3, 4, 1.5, 2, 2.2, 3, 4.5, 3.8, 2.7, 2.5]
N = len(xn)
n = np.arange(N)

f = interpolate.interp1d(n, xn, kind='cubic')

ni = np.linspace(0, N-1, 4000)
xni = f(ni)
print(n[-1])
print(ni[-1])

plt.plot(n, xn, 'ks')
plt.plot(ni, xni)
plt.show()

