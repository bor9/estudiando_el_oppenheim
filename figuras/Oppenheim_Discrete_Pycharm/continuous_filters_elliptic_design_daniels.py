import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import polynomial as P
import math
from matplotlib.patches import Polygon

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

"""
Procedimiento de diseño de un filtro elíptico.

Traducción al lenguaje Python del código en lenguaje TELCOMP 2 del capítulo 5 del libro

R. W. Daniels, Approximation methods for electronic filter design: with applications to passive, active,
and digital networks. McGraw-Hill, 1974.
"""


# complete elliptic integral
def complete_elliptic_integral(k):
    a = np.arctan2(k, np.sqrt(1 - k ** 2))
    o = np.pi / 2
    p = 1
    while True:
        x = 2 / (1 + np.sin(a)) - 1
        y = np.sin(a) * np.sin(o)
        a = np.arctan2(np.sqrt(1 - x ** 2), x)
        o = (o + np.arctan2(y, np.sqrt(1 - y ** 2))) / 2
        e = 1 - 2 * a / np.pi
        p *= (1 + np.cos(a))
        if e < 10e-14:
            break
    x = np.pi / 4 + o / 2
    return p * np.log(np.sin(x) / np.cos(x))


# seno elíptico de Jacobi
def sn(u, k):
    q = np.exp(-np.pi * kk1 / kk)
    v = np.pi / 2 * u / kk
    sn = 0
    j = 0
    while True:
        w = q ** (j + 0.5)
        sn += w * np.sin((2 * j + 1) * v) / (1 - w ** 2)
        j += 1
        if w < 10e-14:
            break
    return sn * 2 * np.pi / (k * kk)


Amax = 0.1
Amin = 40
FB = 20
FH = 26

ee = 10 ** (Amax / 10) - 1
L = np.sqrt((10 ** (Amin / 10) - 1) / ee)
xL = FH / FB

kkl = complete_elliptic_integral(1 / L)
kk1l = complete_elliptic_integral(np.sqrt(1 - (1 / L) ** 2))
kk1 = complete_elliptic_integral(np.sqrt(1 - (1 / xL) ** 2))
kk = complete_elliptic_integral(1 / xL)

N = math.ceil(kk1l * kk / (kkl * kk1))

Ni = 1 if N % 2 == 0 else 0
Nh = N // 2

xxz = np.zeros((Nh, ))
xx = np.zeros((Nh, ))
xxz_alt = np.zeros((Nh, ))
for i in np.arange(Nh):
    u = (2 * (i + 1) - Ni) * kk / N
    s = sn(u, 1 / xL)
    xxz[i] = s
    xx[i] = xL / s

print("polos: {}".format(FB * xx))

# Cálculo de Rn
C = np.prod(1 - xx ** 2) / np.prod(1 - xxz ** 2)
Rn_num = 1 if N % 2 == 0 else [0, 1]
Rn_den = 1
for i in np.arange(Nh):
    Rn_num = P.polymul(Rn_num, [-(xxz[i] ** 2), 0, 1])
    Rn_den = P.polymul(Rn_den, [-(xx[i] ** 2), 0, 1])
w = np.linspace(0, 100, 1000)
ww = w / FB
Rn = C * P.polyval(ww, Rn_num) / P.polyval(ww, Rn_den)

H = 1 / np.sqrt(1 + ee * (Rn ** 2))

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
