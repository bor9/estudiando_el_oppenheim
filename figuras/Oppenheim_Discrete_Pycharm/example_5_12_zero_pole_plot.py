import matplotlib.pyplot as plt
import numpy as np
import math

__author__ = 'ernesto'

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Ceros y polos de H(z)
z1 = (3 / 2) * np.exp(-1j * 3 * np.pi / 4)
z2 = np.conjugate(z1)
p = 1 / 3

# circulo |z|=1
M = 200
thetas = np.linspace(0, 2 * math.pi, M)
# circulo unidad
xc = np.cos(thetas)
yc = np.sin(thetas)

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

# length of the ticks for all subplot (6 pixels)
display_length = 6  # in pixels
# x ticks labels margin
xtm = -0.2
ytm = -0.1
# font size
fontsize = 14
fontsize2 = 10

# estilo de las marcas de los ceros y los polos
zeros_marker_style = dict(marker='o', linestyle='', markersize=8, markerfacecolor='w',
                          markeredgewidth=1.5)
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgewidth=2)

fig = plt.figure(0, figsize=(10, 4), frameon=False)
ax = plt.subplot2grid((4, 12), (0, 0), rowspan=4, colspan=4)

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
# ceros
plt.plot(z1.real, z1.imag, **zeros_marker_style, markeredgecolor='b', zorder=10)
plt.plot(z2.real, z2.imag, **zeros_marker_style, markeredgecolor='r', zorder=10)
plt.plot([0, z1.real], [0, z1.imag], 'k--', zorder=1, lw=1)
plt.plot([0, z2.real], [0, z2.imag], 'k--', zorder=1, lw=1)

# polos
plt.plot(p.real, p.imag, **polos_marker_style, markeredgecolor='k')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')
# etiqueta del polo
plt.text(p.real, -0.3, r'$\frac{1}{3}$', fontsize=fontsize, ha='center', va='baseline')
# titulo
plt.text(xmax_ax, ymax_ax, r'$H(z)$', fontsize=fontsize, ha='right', va='baseline')


plt.axis('off')

###### Sistema de fase minima
ax = plt.subplot2grid((4, 12), (0, 4), rowspan=4, colspan=4)

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
# ceros
plt.plot((1/np.conjugate(z1)).real, (1/np.conjugate(z1)).imag, **zeros_marker_style, markeredgecolor='b', zorder=10)
plt.plot((1/np.conjugate(z2)).real, (1/np.conjugate(z2)).imag, **zeros_marker_style, markeredgecolor='r', zorder=10)
plt.plot([0, (1/np.conjugate(z1)).real], [0, (1/np.conjugate(z1)).imag], 'k--', zorder=1, lw=1)
plt.plot([0, (1/np.conjugate(z2)).real], [0, (1/np.conjugate(z2)).imag], 'k--', zorder=1, lw=1)

# polos
plt.plot(p.real, p.imag, **polos_marker_style, markeredgecolor='k')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')
# etiqueta del polo
plt.text(p.real, -0.3, r'$\frac{1}{3}$', fontsize=fontsize, ha='center', va='baseline')
# titulo
plt.text(xmax_ax, ymax_ax, r'$H_\textrm{min}(z)$', fontsize=fontsize, ha='right', va='baseline')

plt.axis('off')

###### Sistema pasatodo
ax = plt.subplot2grid((4, 12), (0, 8), rowspan=4, colspan=4)

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
# ceros
plt.plot(z1.real, z1.imag, **zeros_marker_style, markeredgecolor='b', zorder=10)
plt.plot(z2.real, z2.imag, **zeros_marker_style, markeredgecolor='r', zorder=10)
plt.plot([0, z1.real], [0, z1.imag], 'k--', zorder=1, lw=1)
plt.plot([0, z2.real], [0, z2.imag], 'k--', zorder=1, lw=1)
# polos
plt.plot((1/np.conjugate(z1)).real, (1/np.conjugate(z1)).imag, **polos_marker_style, markeredgecolor='b', zorder=10)
plt.plot((1/np.conjugate(z2)).real, (1/np.conjugate(z2)).imag, **polos_marker_style, markeredgecolor='r', zorder=10)
plt.plot(xc, yc, 'k--', lw=1)

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')
# titulo
plt.text(xmax_ax, ymax_ax, r'$H_\textrm{ap}(z)$', fontsize=fontsize, ha='right', va='baseline')

plt.axis('off')


# save as pdf image
plt.savefig('example_5_12_zero_pole_plot.pdf', bbox_inches='tight')

plt.show()
