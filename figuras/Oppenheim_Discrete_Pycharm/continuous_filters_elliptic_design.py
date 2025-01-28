import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import polynomial as P
from scipy import special
import math
from matplotlib.patches import Polygon

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

"""
Procedimiento de diseño de un filtro elíptico usando la integral elíptica completa 
y el seno elíptico de Jacobi.
"""

Amax = 0.1
Amin = 40
FB = 20
FH = 26

ee = 10 ** (Amax / 10) - 1
L = np.sqrt((10 ** (Amin / 10) - 1) / ee)
xL = FH / FB

kkl = special.ellipk((1 / L) ** 2)
kk1l = special.ellipk((1 - (1 / L) ** 2) ** 2)
kk1 = special.ellipk((1 - (1 / xL) ** 2) ** 2)
kk = special.ellipk((1 / xL) ** 2)

N = math.ceil(kk1l * kk / (kkl * kk1))

Ni = 1 if N % 2 == 0 else 0
Nh = N // 2

xxz = np.zeros((Nh, ))
xx = np.zeros((Nh, ))
for i in np.arange(Nh):
    u = (2 * (i + 1) - Ni) * kk / N
    s, _, _, _ = special.ellipj(u, (1 / xL) ** 2)
    xxz[i] = s
    xx[i] = xL / s

print("polos: {}".format(FB * xx))

# Cálculo de Rn
C = np.prod((1 - xx ** 2)) / np.prod((1 - xxz ** 2))

Rn_num = 1 if N % 2 == 0 else [0, 1]
Rn_den = 1
for i in np.arange(Nh):
    Rn_num = P.polymul(Rn_num, [-(xxz[i] ** 2), 0, 1])
    Rn_den = P.polymul(Rn_den, [-(xx[i] ** 2), 0, 1])

w = np.linspace(0, 100, 1000)
ww = w / FB
Rn = C * P.polyval(ww, Rn_num) / P.polyval(ww, Rn_den)

H = 1 / (np.sqrt(1 + ee * (Rn ** 2)))

# gráficas

delta_p = 1 - 10 ** (-Amax / 20)
delta_s = 10 ** (-Amin / 20)
ws = FH
wp = FB

### Parámetros de la gráfica
fontsize = 12
# x y ticks labels margin
xtm = -0.07
ytm = -0.04
display_length = 6

xmin_ax = 0
xmax_ax = w[-1]
ymin_ax = 0
ymax_ax = 1.1

grey = [0.9, 0.9, 0.9]

fig = plt.figure(0, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# magnitud de las respuestas en frecuencia
plt.plot(w, H)
# mascara
plt.plot([0, ws], [1, 1], 'k-')
plt.plot([0, wp], [1 - delta_p, 1 - delta_p], 'k-')
plt.plot([ws, w[-1]], [delta_s, delta_s], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p], 'k-')
plt.plot([ws, ws], [delta_s, 1], 'k-')
# región pintada
vert = np.vstack(([0, wp, wp, 0], [0, 0, 1 - delta_p, 1 - delta_p]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
mask_max = ymax_ax
vert = np.vstack(([0, ws, ws, 0], [1, 1, mask_max, mask_max]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax_ax, xmax_ax, ws], [delta_s, delta_s, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)
plt.title(r'$\textrm{{Filtro el\'iptico (N={})}}$'.format(N))
plt.xlabel(r'$\Omega\,(\textrm{rad/s})$')
plt.ylabel(r'$|H_c(j\Omega)|$')

plt.show()
