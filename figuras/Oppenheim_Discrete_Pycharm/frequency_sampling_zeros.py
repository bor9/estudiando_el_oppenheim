import matplotlib.pyplot as plt
import numpy as np
import math

__author__ = 'ernesto'

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


# parámetros
Ne = 6
No = Ne + 1

ke = np.arange(Ne)
ze1 = np.exp(1j * 2 * np.pi * ke / Ne)
ze2 = np.exp(1j * 2 * np.pi * (ke + 0.5)/ Ne)
ko = np.arange(No)
zo1 = np.exp(1j * 2 * np.pi * ko / No)
zo2 = np.exp(1j * 2 * np.pi * (ko + 0.5)/ No)


# circulo |z|=1
M = 200
thetas = np.linspace(0, 2 * math.pi, M)
# circulo unidad
xd = np.cos(thetas)
yd = np.sin(thetas)

# valores maximos y minimos de los ejes
max_ax = 1.5
xmin_ax = -max_ax
xmax_ax = max_ax
ymin_ax = -max_ax
ymax_ax = max_ax

# x ticks labels margin
xtm = -0.18
ytm = -0.1
# font size
fontsize = 14

# estilo de las marcas de los ceros y los polos
zeros_marker_style = dict(marker='o', linestyle='', markersize=8, markerfacecolor='w',
                          markeredgewidth=1.5)

fig = plt.figure(0, figsize=(8, 4), frameon=False)
ax = fig.add_subplot(121)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.gca().set_aspect('equal', adjustable='box')

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

# círculo unidad
plt.plot(xd, yd, 'k', lw=2)
# ceros

plt.plot(ze1.real, ze1.imag, **zeros_marker_style, markeredgecolor='b', zorder=10, label=r'$\alpha=0$')
plt.plot(ze2.real, ze2.imag, **zeros_marker_style, markeredgecolor='r', zorder=10, label=r'$\alpha=1/2$')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')

#leg = plt.legend(loc=1, bbox_to_anchor=(0.4, 1.05), frameon=False, fontsize=fontsize, framealpha=1)
plt.text(0, ymax_ax + 0.2, r'$N\textrm{{ par }}(N={})$'.format(Ne), ha='center', va='baseline', fontsize=fontsize)

plt.axis('off')

ax = fig.add_subplot(122)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.gca().set_aspect('equal', adjustable='box')

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

# círculo unidad
plt.plot(xd, yd, 'k', lw=2)
# ceros

plt.plot(zo1.real, zo1.imag, **zeros_marker_style, markeredgecolor='b', zorder=10, label=r'$\alpha=0$')
plt.plot(zo2.real, zo2.imag, **zeros_marker_style, markeredgecolor='r', zorder=10, label=r'$\alpha=1/2$')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')

leg = plt.legend(loc=1, bbox_to_anchor=(0.15, 1.05), frameon=False, fontsize=fontsize, framealpha=1)
plt.text(0, ymax_ax + 0.2, r'$N\textrm{{ impar }}(N={})$'.format(No), ha='center', va='baseline', fontsize=fontsize)

plt.axis('off')

# save as pdf image
plt.savefig('frequency_sampling_zeros.pdf', bbox_inches='tight')

plt.show()
