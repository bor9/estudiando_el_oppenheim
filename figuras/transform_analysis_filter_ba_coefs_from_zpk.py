import numpy as np
from scipy import signal

# polos
c_k = 0.95 * np.exp(1j * (0.15 * np.pi + 0.02 * np.pi * np.arange(1, 5)))
# construcción del vector de ceros y polos y cálculo de la gananacia: (z, p, k)
# ceros dobles en 1/c_k, 1/c^*_k. polos dobles en c_k, c^*_k.
p = np.concatenate((c_k, c_k, np.conjugate(c_k), np.conjugate(c_k)))
z = np.concatenate((1 / c_k, 1 / c_k, 1 / np.conjugate(c_k), 1 / np.conjugate(c_k)))
k = np.abs(np.prod(c_k)) ** 4
# cálculo de los coeficientes (b, a) a partir de (z, p, k)
b, a = signal.zpk2tf(z, p, k)
