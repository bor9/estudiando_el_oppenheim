import matplotlib.pyplot as plt
import numpy as np

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Parametros del filtro Butterworth en tiempo continuo
N = 5
Omega_c = 1

#
# Diagrama de polos y ceros
#

# Polos de H(s)H(-s)
ks = np.arange(2 * N)
sk = np.exp(1j * np.pi * (2 * ks + N + 1) / (2 * N)) * Omega_c

# círculo donde están los polos
M = 200
thetas = np.linspace(0, 2 * np.pi, M)
xc = Omega_c * np.cos(thetas)
yc = Omega_c * np.sin(thetas)

# valores maximos y minimos de los ejes
max_ax = 1.5
xmin = -max_ax
xmax = max_ax
ymin = -max_ax
ymax = max_ax
# axis parameters
xmin_ax = xmin
xmax_ax = xmax
ymin_ax = ymin
ymax_ax = ymax

display_length = 6

# x ticks labels margin
xtm = -0.23
ytm = -0.1
# font size
fontsize = 14
fontsize2 = 10

fig = plt.figure(1, figsize=(4, 4), frameon=False)
ax = fig.add_subplot(111)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.gca().set_aspect('equal', adjustable='box')

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

# círculo unidad
plt.plot(xc, yc, 'k--', lw=1)
# ceros
zeros_marker_style = dict(marker='o', linestyle='', markersize=8, markerfacecolor='w', markeredgecolor='k',
                          markeredgewidth=1.5)
# polos
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgecolor='k', markeredgewidth=1.5)
plt.plot(sk.real, sk.imag, **polos_marker_style)

# angulo
k = 7
plt.plot([0, sk[k].real], [0, sk[k].imag], 'k--', lw=1)
plt.plot([0, sk[k + 1].real], [0, sk[k + 1].imag], 'k--', lw=1)
r = 0.6
nn = 20
ths = np.linspace(np.angle(sk[k]), np.angle(sk[k + 1]), nn)
xa = r * np.cos(ths)
ya = r * np.sin(ths)
plt.plot(xa, ya, 'k-', lw=1)
plt.text(xa[nn // 2] + 0.14, ya[nn // 2] + 0.04, '$\dfrac{{\pi}}{{{}}}$'.format(N),
         fontsize=fontsize, ha='center', va='center')

# longitud
k = 3
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(sk[k].real, sk[k].imag), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=3, headlength=8, facecolor='black', shrink=0.002))
plt.text(sk[k].real / 2 + 0.05, sk[k].imag / 2, r'$\Omega_c$', fontsize=fontsize, ha='left', va='top')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, r'$\textrm{Re}(s)$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, r'$\textrm{Im}(s)$', fontsize=fontsize, ha='left', va='center')

plt.axis('off')

# save as pdf image
plt.savefig('continuous_filters_butterworth_poles.pdf', bbox_inches='tight')

plt.show()
