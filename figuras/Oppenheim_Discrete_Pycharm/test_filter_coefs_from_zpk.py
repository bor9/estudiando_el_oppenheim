import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from numpy.polynomial import polynomial

__author__ = 'ernesto'

# Construcción del filtro del ejemplo de la sección 5.1.2 del libro
#  A. V. Oppenheim and R. W. Schafer, Discrete-Time Signal Processing. Prentice Hall, 3rd ed., 2009.

# Polos
c_k = 0.95 * np.exp(1j * (0.15 * np.pi + 0.02 * np.pi * np.arange(1, 5)))

# Construcción del vector de polos y ceros
# Ceros dobles: 1/c_k, 1/c^*_k. Polos dobles: c_k, c^*_k.
p = np.concatenate((c_k, c_k, np.conjugate(c_k), np.conjugate(c_k)))
z = np.concatenate((1 / c_k, 1 / c_k, 1 / np.conjugate(c_k), 1 / np.conjugate(c_k)))
k = np.abs(np.prod(c_k)) ** 4 # el 4 es porque los polos son dobles.

b, a = signal.zpk2tf(z, p, k)

# Lo mismo como multiplicación de polinomios
aa = [1]
bb = [1]
for c in c_k:
    aa = polynomial.polymul(aa, [-1, c])
    aa = polynomial.polymul(aa, [-1, c])
    aa = polynomial.polymul(aa, [-1, np.conjugate(c)])
    aa = polynomial.polymul(aa, [-1, np.conjugate(c)])
    bb = polynomial.polymul(bb, [-c, 1])
    bb = polynomial.polymul(bb, [-c, 1])
    bb = polynomial.polymul(bb, [-np.conjugate(c), 1])
    bb = polynomial.polymul(bb, [-np.conjugate(c), 1])

aa = aa.real
bb = bb.real

nw = 1024

w, H1 = signal.freqz(b, a, nw)
w, H2 = signal.freqz(bb, aa, nw)

_, grdH = signal.group_delay((b, a), nw)

# Cálculo manual del retardo de grupo
grdH_alt = -np.diff(np.unwrap(np.angle(H2))) * nw / np.pi

print(b-bb)
print(a-aa)

fig = plt.figure(0)
plt.subplot(311)
plt.plot(w, np.abs(H2))
plt.subplot(312)
plt.plot(w, np.unwrap(np.angle(H2)))
plt.subplot(313)
plt.plot(grdH_alt)
plt.show()
