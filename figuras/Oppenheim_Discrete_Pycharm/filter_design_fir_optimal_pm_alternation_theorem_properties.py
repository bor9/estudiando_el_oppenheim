import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from matplotlib.patches import Polygon
from numpy.polynomial import polynomial as P
from scipy.signal import argrelextrema

__author__ = 'ernesto'

# use latex
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


# Especificaciones del filtro
delta_1 = 0.05
delta_2 = 0.1
wp = 0.2765 * np.pi
ws = 0.4 * np.pi


# largo del filtro
M = 15

wb = [0, wp, ws, np.pi]
K = delta_1 / delta_2

h = signal.remez(M, wb, [1, 0], [1 / K, 1], fs=2*np.pi)

# Cálculo de la respuesta en frecuencia
nw = 1024
w = np.linspace(0, np.pi, nw)
_, H = signal.freqz(h, 1, worN=w, whole=False, plot=None, fs=2*np.pi)
# eliminación del componente de fase lineal
H = np.real(H * np.exp(1j * w * (M - 1) / 2))

# se redefine wp a ser un poco menor
wp = 0.26 * np.pi

# cálculo del error de aproximación para obtener los indices de las frecuencias de alternancia
# respuesta deseada
Hd = np.zeros(H.shape)
Hd[w <= wp] = 1
# función de ponderación
W = np.ones(H.shape)
W[w <= wp] = 1 / K
# error de aproximación
E = W * (Hd - H)
E[np.logical_and(w > wp, w < ws)] = 'nan'
# indices de las frecuencias de alternancia (maximos del valor absoluto del error)
alt_idx = argrelextrema(np.abs(E), np.greater)[0]
#alt_idx = np.concatenate([[0], alt_idx, [nw - 1]])

# cálculo de delta para el filtro obtenido
delta = np.max(np.abs(H[w >= ws]))
# se redefinen delta 1 y delta 2
delta_2 = delta
delta_1 = K * delta


### Parámetros de la gráfica
fontsize = 12
# x y ticks labels margin
xtm = -0.065
ytm = -0.04
display_length = 6

xmin = 0
xmax = np.pi
dx = 0.2
xmax_ax = xmax + dx
xmin_ax = xmin - dx
ymin_ax = -0.2
ymax_ax = 1.3

mask_amp = 0.08

grey = [0.9, 0.9, 0.9]

alt_marker_style = dict(marker='o', linestyle='', markersize=5, markerfacecolor='w',
                        markeredgewidth=1.5, markeredgecolor='r')

fig = plt.figure(0, figsize=(8, 5), frameon=False)
ax = plt.subplot2grid((32, 4), (0, 0), rowspan=22, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
# filtro
plt.plot(w, H, 'r-', lw=1.5, zorder=10)
# mascara
plt.plot([0, wp], [1 + delta_1, 1 + delta_1], 'k-', lw=1)
plt.plot([0, wp], [1 - delta_1, 1 - delta_1], 'k-', lw=1)
plt.plot([ws, np.pi], [delta_2, delta_2], 'k-', lw=1)
plt.plot([ws, np.pi], [-delta_2, -delta_2], 'k-', lw=1)

plt.plot([w[alt_idx], w[alt_idx]], [np.zeros(alt_idx.shape), H[alt_idx]], 'k--', lw=1)
plt.plot(w[alt_idx], H[alt_idx], **alt_marker_style, zorder=15)
plt.plot(w[0], H[0], **alt_marker_style, zorder=15)
plt.plot(w[-1], H[-1], **alt_marker_style, zorder=15)
for k, ai in enumerate(alt_idx):
    if H[ai] > 0:
        plt.text(w[ai], xtm, '$\omega_{}$'.format(k + 1), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(w[ai], 0.03, '$\omega_{}$'.format(k + 1), fontsize=fontsize, ha='center', va='baseline')

# bandas
plt.plot([wp, wp], [0, 1 + delta_1], 'k--', lw=1)
plt.text(wp, xtm, '$\omega_p$', fontsize=fontsize, ha='center', va='baseline')
plt.plot([ws, ws], [-delta_2, delta_2], 'k--', lw=1)
plt.text(ws - 0.02, xtm, '$\omega_s$', fontsize=fontsize, ha='right', va='baseline')
plt.plot([np.pi, np.pi], [-delta_2, delta_2], 'k--', lw=1)

# región pintada
vert = np.vstack(([0, wp, wp, 0], [1 - delta_1 - mask_amp, 1 - delta_1 - mask_amp, 1 - delta_1, 1 - delta_1]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
vert = np.vstack(([0, wp, wp, 0], [1 + delta_1, 1 + delta_1, 1 + delta_1 + mask_amp, 1 + delta_1 + mask_amp]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [delta_2, delta_2, delta_2 + mask_amp, delta_2 + mask_amp]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)
vert = np.vstack(([ws, xmax, xmax, ws], [-delta_2, -delta_2, -delta_2 - mask_amp, -delta_2 - mask_amp]))
p4 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p4)

# etiquetas
# eje x
plt.text(np.pi + 0.03, xtm, '$\pi$', fontsize=fontsize, ha='left', va='baseline')
# eje y
plt.text(ytm, 1 - delta_1, '$1-\delta_1$', fontsize=fontsize, ha='right', va='top')
plt.text(ytm, 1 + delta_1, '$1+\delta_1$', fontsize=fontsize, ha='right', va='bottom')
plt.plot([0, ytl], [1, 1], 'k-', lw=1)
plt.text(ytm, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, ytl], [delta_2, delta_2], 'k-', lw=1)
plt.text(ytm, delta_2, '$\delta_2$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, ytl], [-delta_2, -delta_2], 'k-', lw=1)
plt.text(ytm, -delta_2, '$-\delta_2$', fontsize=fontsize, ha='right', va='center')
# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')

plt.text(0.06, ymax_ax, '$A_e(e^{j\omega})$', fontsize=fontsize, ha='left', va='center')

plt.axis('off')

#####################################################################
#####################################################################
#####################################################################


ymax_ax = 3 * delta_2
ymin_ax = -2.4 * delta_2

ax = plt.subplot2grid((32, 4), (23, 0), rowspan=8, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(w, E, 'r-', lw=1.5, zorder=10)
plt.plot(w[alt_idx], E[alt_idx], **alt_marker_style, zorder=15)
plt.plot(w[0], E[0], **alt_marker_style, zorder=15)
plt.plot(w[-1], E[-1], **alt_marker_style, zorder=15)

# mascara
plt.plot([0, wp], [delta_2, delta_2], 'k-', lw=1)
plt.plot([0, wp], [-delta_2, -delta_2], 'k-', lw=1)
plt.plot([ws, np.pi], [delta_2, delta_2], 'k-', lw=1)
plt.plot([ws, np.pi], [-delta_2, -delta_2], 'k-', lw=1)
# región pintada
vert = np.vstack(([0, wp, wp, 0], [-delta_2 - mask_amp, -delta_2 - mask_amp, -delta_2, -delta_2]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
vert = np.vstack(([0, wp, wp, 0], [delta_2, delta_2, delta_2 + mask_amp, delta_2 + mask_amp]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [delta_2, delta_2, delta_2 + mask_amp, delta_2 + mask_amp]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)
vert = np.vstack(([ws, xmax, xmax, ws], [-delta_2, -delta_2, -delta_2 - mask_amp, -delta_2 - mask_amp]))
p4 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p4)

# etiquetas
# bandas
plt.plot([wp, wp], [-delta_2, delta_2], 'k--', lw=1)
plt.text(wp + 0.02, xtm, '$\omega_p$', fontsize=fontsize, ha='left', va='baseline')
plt.plot([ws, ws], [-delta_2, delta_2], 'k--', lw=1)
plt.text(ws - 0.02, xtm, '$\omega_s$', fontsize=fontsize, ha='right', va='baseline')
plt.plot([np.pi, np.pi], [-delta_2, delta_2], 'k--', lw=1)
plt.text(np.pi + 0.03, xtm, '$\pi$', fontsize=fontsize, ha='left', va='baseline')

# eje y
plt.text(ytm, delta_2+0.01, '$\delta_2$', fontsize=fontsize, ha='right', va='bottom')
plt.text(ytm, -delta_2-0.02, '$-\delta_2$', fontsize=fontsize, ha='right', va='top')
# etiquetas de los ejes
plt.text(0.04, ymax_ax-0.02, '$E(\omega)$', fontsize=fontsize, ha='left', va='center')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.axis('off')

plt.savefig('filter_design_fir_optimal_pm_alternation_theorem_properties.pdf', bbox_inches='tight')

plt.show()
