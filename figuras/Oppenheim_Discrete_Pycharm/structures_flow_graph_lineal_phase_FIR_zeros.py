import matplotlib.pyplot as plt
import numpy as np
import math

__author__ = 'ernesto'

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Ceros y polos de H(z)


# z1: cero complejo generico
z1 = 0.7 * np.exp(1j * np.pi / 3)
# z2: cero real
z2 = 0.6
# z3: cero complejo sobre el círculo unidad
z3 = np.exp(1j * 4 * np.pi / 5)
# z4: cero real sobre el círculo unidad
z4 = -1

# circulo |z|=1
M = 200
thetas = np.linspace(0, 2 * math.pi, M)
# circulo unidad
xc = np.cos(thetas)
yc = np.sin(thetas)

# axis parameters
xmin_ax = -1.5
xmax_ax = 2
ymin_ax = -1.5
ymax_ax = 1.5

# x ticks labels margin
xtm = -0.22
ytm = -0.1
# font size
fontsize = 14

# estilo de las marcas de los ceros y los polos
zeros_marker_style = dict(marker='o', linestyle='', markersize=8, markerfacecolor='w',
                          markeredgewidth=1.5)
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgewidth=2)

fig = plt.figure(0, figsize=(5, 4), frameon=False)
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
plt.plot(xc, yc, 'k', lw=2)
dt = 0.12
# z1
plt.plot(z1.real, z1.imag, **zeros_marker_style, markeredgecolor='r', zorder=10)
plt.text(z1.real + dt, z1.imag, '$z_1$', fontsize=fontsize, ha='left', va='center')
plt.plot((1 / z1).real, -(1 / z1).imag, **zeros_marker_style, markeredgecolor='r', zorder=10)
plt.text((1 / z1).real + dt, -(1 / z1).imag, r'$\dfrac{1}{z_1^*}$', fontsize=fontsize, ha='left', va='center')
plt.plot([0, (1 / z1).real], [0, -(1 / z1).imag], 'k--', zorder=1, lw=1)
plt.plot(z1.real, -z1.imag, **zeros_marker_style, markeredgecolor='r', zorder=10)
plt.text(z1.real + dt, -z1.imag, '$z_1^*$', fontsize=fontsize, ha='left', va='center')
plt.plot((1 / z1).real, (1 / z1).imag, **zeros_marker_style, markeredgecolor='r', zorder=10)
plt.text((1 / z1).real + dt, (1 / z1).imag, r'$\dfrac{1}{z_1}$', fontsize=fontsize, ha='left', va='center')
plt.plot([0, (1 / z1).real], [0, (1 / z1).imag], 'k--', zorder=1, lw=1)
# z2
plt.plot(z2.real, z2.imag, **zeros_marker_style, markeredgecolor='b', zorder=10)
plt.text(z2, xtm, '$z_2$', fontsize=fontsize, ha='center', va='baseline')
plt.plot((1/z2).real, (1/z2).imag, **zeros_marker_style, markeredgecolor='b', zorder=10)
plt.text(1/z2, 0.3, r'$\dfrac{1}{z_2}$', fontsize=fontsize, ha='center', va='baseline')
# z3
plt.plot(z3.real, z3.imag, **zeros_marker_style, markeredgecolor='g', zorder=10)
plt.text(z3.real - dt, z3.imag, '$z_3$', fontsize=fontsize, ha='right', va='bottom')
plt.plot(z3.real, -z3.imag, **zeros_marker_style, markeredgecolor='g', zorder=10)
plt.text(z3.real - dt, -z3.imag, '$z_3^*$', fontsize=fontsize, ha='right', va='top')
# z4
plt.plot(z4.real, z4.imag, **zeros_marker_style, markeredgecolor='m', zorder=10)
dt = 0.08
plt.text(z4.real - dt, z4.imag + dt, '$z_4$', fontsize=fontsize, ha='right', va='bottom')

plt.plot(xc, yc, 'k--', lw=1)

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')

plt.axis('off')

# save as pdf image
plt.savefig('structures_flow_graph_lineal_phase_FIR_zeros.pdf', bbox_inches='tight')

plt.show()
