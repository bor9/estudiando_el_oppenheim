import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.patches import Polygon

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
# polos
ps = [0.4, 0.7, 1.3]
pn = ['a', 'b', 'c']

# circulo |z|=1
M = 200
thetas = np.linspace(0, 2 * math.pi, M)
# circulo unidad
xd = np.cos(thetas)
yd = np.sin(thetas)

# valores maximos y minimos de los ejes
max_ax = 1.8
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
xtm = -0.25
ytm = -0.1
# font size
fontsize = 14
fontsize2 = 10

# estilo de las marcas de los ceros y los polos
zeros_marker_style = dict(marker='o', linestyle='', markersize=8, markerfacecolor='w', markeredgecolor='k',
                          markeredgewidth=1.5)
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgecolor='k', markeredgewidth=2)

fig = plt.figure(0, figsize=(11, 8), frameon=False)
ax = plt.subplot2grid((8, 12), (2, 0), rowspan=4, colspan=4)

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
# polos
for i in np.arange(len(ps)):
    plt.plot(ps[i], 0, **polos_marker_style)
    plt.text(ps[i], xtm, '${}$'.format(pn[i]), fontsize=fontsize, ha='center', va='baseline')

# punto 1
ms = 7
# 1
plt.plot(1, 0, 'k.', ms=ms)

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')
# circle labels
plt.text(1.05, xtm, '$1$', fontsize=fontsize, ha='left', va='baseline')
# titulo
dx = 0.2
plt.text(0, ymax_ax + dx, '${\\rm Diagrama\;de\;polos\;y\;ceros}$', fontsize=fontsize, ha='center', va='bottom',
         ma='center')

plt.axis('off')

####### Secuencia hacia adelante
ax = plt.subplot2grid((8, 12), (0, 4), rowspan=4, colspan=4)

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
plt.plot(xd, yd, 'k', lw=2, zorder=3)
# frontera de la ROC
ind = 2
xr = ps[ind] * np.cos(thetas)
yr = ps[ind] * np.sin(thetas)
plt.plot(xr, yr, 'k--', lw=1, zorder=3)

# polos
for i in np.arange(len(ps)):
    plt.plot(ps[i], 0, **polos_marker_style, zorder=3)
    if i != ind:
        plt.text(ps[i], xtm, '${}$'.format(pn[i]), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(ps[i]+0.05, xtm, '${}$'.format(pn[i]), fontsize=fontsize, ha='left', va='baseline')

# punto 1
ms = 7
# 1
plt.plot(1, 0, 'k.', ms=ms, zorder=3)

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')
# circle labels
plt.text(1.05, xtm, '$1$', fontsize=fontsize, ha='left', va='baseline')

# relleno de la ROC
v = xmax - 0.2
vert = np.vstack(([v, -v, -v, v], [v, v, -v, -v]))
p2 = Polygon(vert.T, fc=(0.9, 0.9, 0.9), alpha=1, zorder=1)
ax.add_artist(p2)
vert = np.vstack((xr, yr))
p1 = Polygon(vert.T, fc='w', alpha=1, zorder=2)
ax.add_artist(p1)
# titulo
plt.text(0, ymax_ax + dx, '${\\rm Secuencia\;hacia\;adelante}$', fontsize=fontsize, ha='center', va='bottom',
         ma='center')
plt.axis('off')

####### Secuencia hacia atras
ax = plt.subplot2grid((8, 12), (0, 8), rowspan=4, colspan=4)

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
plt.plot(xd, yd, 'k', lw=2, zorder=3)
# frontera de la ROC
ind = 0
xr = ps[ind] * np.cos(thetas)
yr = ps[ind] * np.sin(thetas)
plt.plot(xr, yr, 'k--', lw=1, zorder=3)

# polos
for i in np.arange(len(ps)):
    plt.plot(ps[i], 0, **polos_marker_style, zorder=3)
    if i != ind:
        plt.text(ps[i], xtm, '${}$'.format(pn[i]), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(ps[i]+0.05, xtm, '${}$'.format(pn[i]), fontsize=fontsize, ha='left', va='baseline')

# punto 1
ms = 7
# 1
plt.plot(1, 0, 'k.', ms=ms, zorder=3)

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')
# circle labels
plt.text(1.05, xtm, '$1$', fontsize=fontsize, ha='left', va='baseline')

# relleno de la ROC
vert = np.vstack((xr, yr))
p1 = Polygon(vert.T, fc=(0.9, 0.9, 0.9), alpha=1, zorder=2)
ax.add_artist(p1)
# titulo
plt.text(0, ymax_ax + dx, '${\\rm Secuencia\;hacia\;atr\\acute{a}s}$', fontsize=fontsize, ha='center', va='bottom',
         ma='center')
plt.axis('off')

####### Secuencia hacia ambos lados
ax = plt.subplot2grid((8, 12), (4, 4), rowspan=4, colspan=4)

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
plt.plot(xd, yd, 'k', lw=2, zorder=3)
# frontera de la ROC

xr1 = ps[0] * np.cos(thetas)
yr1 = ps[0] * np.sin(thetas)
plt.plot(xr1, yr1, 'k--', lw=1, zorder=3)
xr2 = ps[1] * np.cos(thetas)
yr2 = ps[1] * np.sin(thetas)
plt.plot(xr2, yr2, 'k--', lw=1, zorder=3)

ind = 2
# polos
for i in np.arange(len(ps)):
    plt.plot(ps[i], 0, **polos_marker_style, zorder=3)
    if i == ind:
        plt.text(ps[i], xtm, '${}$'.format(pn[i]), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(ps[i]+0.05, xtm, '${}$'.format(pn[i]), fontsize=fontsize, ha='left', va='baseline')

# punto 1
ms = 7
# 1
plt.plot(1, 0, 'k.', ms=ms, zorder=3)

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')
# circle labels
plt.text(1.05, xtm, '$1$', fontsize=fontsize, ha='left', va='baseline')

# relleno de la ROC
vert = np.vstack((xr2, yr2))
p2 = Polygon(vert.T, fc=(0.9, 0.9, 0.9), alpha=1, zorder=1)
ax.add_artist(p2)
vert = np.vstack((xr1, yr1))
p2 = Polygon(vert.T, fc='w', alpha=1, zorder=2)
ax.add_artist(p2)


# titulo
plt.text(0, ymin_ax - dx, '${\\rm Secuencia\;hacia\;ambos\;lados}$', fontsize=fontsize, ha='center', va='top',
         ma='center')
plt.axis('off')

####### Secuencia hacia ambos lados
ax = plt.subplot2grid((8, 12), (4, 8), rowspan=4, colspan=4)

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
plt.plot(xd, yd, 'k', lw=2, zorder=3)
# frontera de la ROC

xr2 = ps[1] * np.cos(thetas)
yr2 = ps[1] * np.sin(thetas)
plt.plot(xr2, yr2, 'k--', lw=1, zorder=3)
xr3 = ps[2] * np.cos(thetas)
yr3 = ps[2] * np.sin(thetas)
plt.plot(xr3, yr3, 'k--', lw=1, zorder=3)

ind = 0
# polos
for i in np.arange(len(ps)):
    plt.plot(ps[i], 0, **polos_marker_style, zorder=3)
    if i == ind:
        plt.text(ps[i], xtm, '${}$'.format(pn[i]), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(ps[i]+0.05, xtm, '${}$'.format(pn[i]), fontsize=fontsize, ha='left', va='baseline')

# punto 1
ms = 7
# 1
plt.plot(1, 0, 'k.', ms=ms, zorder=3)

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')
# circle labels
plt.text(1.05, xtm, '$1$', fontsize=fontsize, ha='left', va='baseline')

# relleno de la ROC
vert = np.vstack((xr3, yr3))
p2 = Polygon(vert.T, fc=(0.9, 0.9, 0.9), alpha=1, zorder=1)
ax.add_artist(p2)
vert = np.vstack((xr2, yr2))
p2 = Polygon(vert.T, fc='w', alpha=1, zorder=2)
ax.add_artist(p2)


# titulo
plt.text(0, ymin_ax - dx, '${\\rm Secuencia\;hacia\;ambos\;lados}$', fontsize=fontsize, ha='center', va='top',
         ma='center')
plt.axis('off')
















# save as pdf image
plt.savefig('z_transform_roc_posibilities.pdf', bbox_inches='tight')

plt.show()
