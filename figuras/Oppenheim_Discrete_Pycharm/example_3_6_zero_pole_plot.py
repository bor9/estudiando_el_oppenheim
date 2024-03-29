import matplotlib.pyplot as plt
import numpy as np
import math

__author__ = 'ernesto'

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


# auxiliar function for plot ticks of equal length in x and y axis despite its scales.
def convert_display_to_data_coordinates(transData, length=10):
    # create a transform which will take from display to data coordinates
    inv = transData.inverted()
    # transform from display coordinates to data coordinates in x axis
    data_coords = inv.transform([(0, 0), (length, 0)])
    # get the length of the segment in data units
    yticks_len = data_coords[1, 0] - data_coords[0, 0]
    # transform from display coordinates to data coordinates in y axis
    data_coords = inv.transform([(0, 0), (0, length)])
    # get the length of the segment in data units
    xticks_len = data_coords[1, 1] - data_coords[0, 1]
    return xticks_len, yticks_len


# parámetros
a = 0.8
N = 16

# circulo |z|=1
M = 200
thetas = np.linspace(0, 2 * math.pi, M)
# circulo unidad
xd = np.cos(thetas)
yd = np.sin(thetas)
# ceros
k = np.arange(1, N)
z = a * np.exp(1j * 2 * np.pi * k / N)
# círculo donde están los ceros
xz = a * np.cos(thetas)
yz = a * np.sin(thetas)

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
xtm = -0.18
ytm = -0.1
# font size
fontsize = 14
fontsize2 = 10

fig = plt.figure(0, figsize=(4, 4), frameon=False)
ax = fig.add_subplot(111)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.gca().set_aspect('equal', adjustable='box')

# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

# círculo unidad
plt.plot(xd, yd, 'k', lw=2)
# ceros
zeros_marker_style = dict(marker='o', linestyle='', markersize=8, markerfacecolor='w', markeredgecolor='k',
                          markeredgewidth=1.5)
plt.plot(z.real, z.imag, **zeros_marker_style, zorder=10)
plt.plot(xz, yz, 'k--', lw=1)
# polos
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgecolor='k', markeredgewidth=2)
plt.plot(0, 0, **polos_marker_style)
# multiplicidad del polo
plt.text(-0.07, 0.07, '{:d}'.format(N - 1), fontsize=fontsize2, ha='right', va='baseline')

# puntos a y 1
ms = 7
# 1
plt.plot(1, 0, 'k.', ms=ms)
# a
plt.plot(a, 0, 'k.', ms=ms)

# angulo
plt.plot([0, z[0].real], [0, z[0].imag], 'k--', lw=1)
r = 0.5
nn = 20
ths = np.linspace(0, np.angle(z[0]), nn)
xa = r * np.cos(ths)
ya = r * np.sin(ths)
plt.plot(xa, ya, 'k-', lw=1)
plt.annotate('$\dfrac{\pi}{8}$', xytext=(1.2, 0.4), xycoords='data', xy=(xa[int(nn/2)], ya[int(nn/2)]),
             textcoords='data', fontsize=fontsize, va="center", ha="left",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0, 0.4),
                             patchA=None, patchB=None, shrinkA=4, shrinkB=0,
                             connectionstyle="angle3,angleA={:d},angleB={:d}".format(-145, 10)))

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')

# circle labels
plt.text(1.05, xtm, '$1$', fontsize=fontsize, ha='left', va='baseline')
plt.text(a + 0.04, xtm, '$a$', fontsize=fontsize, ha='left', va='baseline')

plt.axis('off')

# save as pdf image
plt.savefig('example_3_6_zero_pole_plot.pdf', bbox_inches='tight')

plt.show()
