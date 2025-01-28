import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import polynomial as P
from scipy import signal

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# parámetros
a = 0.1
b = 1
Td = 1

# Filtro analógico H_c(jw)
ac = [1, 2 * a, a ** 2 + b ** 2]
bc = [1, a ** 2]
Wmax = np.pi / Td
W = np.linspace(0, Wmax, 500)

_, Hc = signal.freqs(bc, ac, W)

# Filtro digital mediante invarianza al impulso
a1 = [1, -2 * np.exp(-a * Td) * np.cos(b * Td), np.exp(-2 * a * Td)]
b1 = [1, -np.exp(-a * Td) * np.cos(b * Td)]

w = np.linspace(0, np.pi, 500)
_, H1 = signal.freqz(b1, a1, w)

# Filtro digital mediante invarianza al escalón
a2 = a1
c = np.exp(-a * Td)
b2 = np.array([0, a * (1 - c * np.cos(b * Td)) + b * c * np.sin(b * Td),
               c * (a * (c - np.cos(b * Td)) - b * np.sin(b * Td))]) / (a ** 2 + b ** 2)

_, H2 = signal.freqz(b2, a2, w)

# Comprobación
ejw = np.exp(-1j * w)
H2_alt = a / (a ** 2 + b ** 2) \
         - P.polyval(ejw, [1, -1]) / (2 * (a - 1j * b) * P.polyval(ejw, [1, -np.exp(-(a - 1j * b) * Td)]))\
         - P.polyval(ejw, [1, -1]) / (2 * (a + 1j * b) * P.polyval(ejw, [1, -np.exp(-(a + 1j * b) * Td)]))

print('Diferencia: {}'.format(np.sum(np.abs(H2 - H2_alt))))

fs = 11

fig = plt.figure(0, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=2)
plt.xlim(w[0], w[-1])
plt.plot(W * Td, np.absolute(Hc) / Td, label=r'$\dfrac{1}{T}H_c(j\Omega)$')
plt.plot(w, np.absolute(H1), label=r'$H_1(e^{j\omega})$')
plt.plot(w, np.absolute(H2) / Td, label=r'$\dfrac{1}{T}H_2(e^{j\omega})$')

plt.legend(loc='upper right', fontsize=fs, frameon=False, framealpha=1)
plt.xlabel(r'$\omega\,(\textrm{rad}),\,\Omega=\omega/T\,(\textrm{rad/s})$', fontsize=fs)
plt.ylabel(r'$\textrm{Magnitud}$', fontsize=fs)

ymax = ax.get_ylim()[1]
plt.ylim(0, ymax)
plt.text(0.2, ymax - 0.3, '$T=1$', ha='left', va='top', fontsize=fs)


ax = plt.subplot2grid((4, 4), (0, 2), rowspan=4, colspan=2)
plt.xlim(w[0], w[-1])
plt.plot(w, np.angle(H1))
plt.plot(w, np.angle(H2))
plt.plot(W * Td, np.angle(Hc))

ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
plt.xlabel(r'$\omega\,(\textrm{rad}),\,\Omega=\omega/T\,(\textrm{rad/s})$', fontsize=fs)
plt.ylabel(r'$\textrm{Fase (rad)}$', fontsize=fs)

# save as pdf image
plt.savefig('filter_design_problem_07_01.pdf', bbox_inches='tight')

plt.show()
