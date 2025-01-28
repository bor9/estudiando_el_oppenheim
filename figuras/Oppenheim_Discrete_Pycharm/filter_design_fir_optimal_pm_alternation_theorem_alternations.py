import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from matplotlib.patches import Polygon
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
K = delta_1 / delta_2

# largo del filtro
M = 15

# Para el cálculo de la respuesta en frecuencia
nw = 2 * 1024
w = np.linspace(0, np.pi, nw)

### Parámetros de la gráfica
fontsize = 11
# x y ticks labels margin
xtm = -0.08
ytm = -0.04
display_length = 6

xmin = 0
xmax = np.pi
dx = 0.2
xmax_ax = xmax + 1.5 * dx
xmin_ax = xmin - dx
ymin_ax = -0.2
ymax_ax = 1.3

mask_amp = 0.08

grey = [0.9, 0.9, 0.9]

alt_marker_style = dict(marker='o', linestyle='', markersize=5, markerfacecolor='w',
                        markeredgewidth=1.5, markeredgecolor='r')

fig = plt.figure(0, figsize=(8.3, 5.8), frameon=False)

#############
# Filtro 1: #
#############

wp = 0.2765 * np.pi
ws = 0.4 * np.pi

wb = [0, wp, ws, np.pi]

h = signal.remez(M, wb, [1, 0], [1 / K, 1], fs=2*np.pi)

# Cálculo de la respuesta en frecuencia
_, H = signal.freqz(h, 1, worN=w, whole=False, plot=None, fs=2*np.pi)
# eliminación del componente de fase lineal
H = np.real(H * np.exp(1j * w * (M - 1) / 2))

# cálculo del error de aproximación para obtener los indices de las frecuencias de alternancia
# respuesta deseada
Hd = np.zeros(H.shape)
Hd[w <= wp] = 1
# función de ponderación
W = np.ones(H.shape)
W[w <= wp] = 1 / K
# error de aproximación
E = W * (Hd - H)
E[np.logical_and(w > wp, w < ws)] = 0
# indices de las frecuencias de alternancia (maximos del valor absoluto del error)
alt_idx = argrelextrema(np.abs(E), np.greater)[0]

# cálculo de delta para el filtro obtenido
delta = np.max(np.abs(H[w >= ws]))
# se redefinen delta 1 y delta 2
delta_2 = delta
delta_1 = K * delta

###################################################################

ax = plt.subplot2grid((4, 4), (0, 0), rowspan=2, colspan=2)
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
plt.plot([wp, wp], [0, 1 - delta_1], 'k--', lw=1)
plt.plot([ws, ws], [-delta_2, delta_2], 'k--', lw=1)
plt.plot([np.pi, np.pi], [-delta_2, delta_2], 'k--', lw=1)

plt.plot(w[alt_idx], H[alt_idx], **alt_marker_style, zorder=15)
plt.plot(w[0], H[0], **alt_marker_style, zorder=15)
plt.plot(w[-1], H[-1], **alt_marker_style, zorder=15)

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
plt.text(wp, xtm, '$\omega_p$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ws - 0.02, xtm, '$\omega_s$', fontsize=fontsize, ha='right', va='baseline')
plt.text(np.pi + 0.04, xtm, '$\pi$', fontsize=fontsize, ha='left', va='baseline')
# eje y
plt.plot([0, ytl], [1, 1], 'k-', lw=1)
plt.plot([0, ytl], [delta_2, delta_2], 'k-', lw=1)
plt.plot([0, ytl], [-delta_2, -delta_2], 'k-', lw=1)
# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')

plt.text(0.1, ymax_ax, '$A_e(e^{j\omega})$', fontsize=fontsize, ha='left', va='center')

plt.axis('off')

#############
# Filtro 2: #
#############

wp = 0.27 * np.pi
ws = 0.4 * np.pi

wb = [0, wp, ws, np.pi]

h = signal.remez(M, wb, [1, 0], [1 / K, 1], fs=2*np.pi)

# Cálculo de la respuesta en frecuencia
_, H = signal.freqz(h, 1, worN=w, whole=False, plot=None, fs=2*np.pi)
# eliminación del componente de fase lineal
H = np.real(H * np.exp(1j * w * (M - 1) / 2))

# cálculo del error de aproximación para obtener los indices de las frecuencias de alternancia
# respuesta deseada
Hd = np.zeros(H.shape)
Hd[w <= wp] = 1
# función de ponderación
W = np.ones(H.shape)
W[w <= wp] = 1 / K
# error de aproximación
E = W * (Hd - H)
E[np.logical_and(w > wp, w < ws)] = 0
# indices de las frecuencias de alternancia (maximos del valor absoluto del error)
alt_idx = argrelextrema(np.abs(E), np.greater)[0]

# cálculo de delta para el filtro obtenido
delta = np.max(np.abs(H[w >= ws]))
# se redefinen delta 1 y delta 2
delta_2 = delta
delta_1 = K * delta


###################################################################

ax = plt.subplot2grid((4, 4), (0, 2), rowspan=2, colspan=2)
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
plt.plot([wp, wp], [0, 1 - delta_1], 'k--', lw=1)
plt.plot([ws, ws], [-delta_2, delta_2], 'k--', lw=1)
plt.plot([np.pi, np.pi], [-delta_2, delta_2], 'k--', lw=1)

plt.plot(w[alt_idx], H[alt_idx], **alt_marker_style, zorder=15)
plt.plot(w[-1], H[-1], **alt_marker_style, zorder=15)

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
plt.text(wp, xtm, '$\omega_p$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ws - 0.02, xtm, '$\omega_s$', fontsize=fontsize, ha='right', va='baseline')
plt.text(np.pi + 0.04, xtm, '$\pi$', fontsize=fontsize, ha='left', va='baseline')
# eje y
plt.plot([0, ytl], [1, 1], 'k-', lw=1)
plt.plot([0, ytl], [delta_2, delta_2], 'k-', lw=1)
plt.plot([0, ytl], [-delta_2, -delta_2], 'k-', lw=1)
# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')

plt.text(0.1, ymax_ax, '$A_e(e^{j\omega})$', fontsize=fontsize, ha='left', va='center')

plt.axis('off')

#############
# Filtro 3: #
#############

wp = 0.28 * np.pi
ws = 0.402 * np.pi

wb = [0, wp, ws, np.pi]

h = signal.remez(M, wb, [1, 0], [1 / K, 1], fs=2*np.pi)

# Cálculo de la respuesta en frecuencia
_, H = signal.freqz(h, 1, worN=w, whole=False, plot=None, fs=2*np.pi)
# eliminación del componente de fase lineal
H = np.real(H * np.exp(1j * w * (M - 1) / 2))

# cálculo del error de aproximación para obtener los indices de las frecuencias de alternancia
# respuesta deseada
Hd = np.zeros(H.shape)
Hd[w <= wp] = 1
# función de ponderación
W = np.ones(H.shape)
W[w <= wp] = 1 / K
# error de aproximación
E = W * (Hd - H)
E[np.logical_and(w > wp, w < ws)] = 0
# indices de las frecuencias de alternancia (maximos del valor absoluto del error)
alt_idx = argrelextrema(np.abs(E), np.greater)[0]

# cálculo de delta para el filtro obtenido
delta = np.max(np.abs(H[w >= ws]))
# se redefinen delta 1 y delta 2
delta_2 = delta
delta_1 = K * delta

###################################################################

ax = plt.subplot2grid((4, 4), (2, 0), rowspan=2, colspan=2)
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
plt.plot([wp, wp], [0, 1 - delta_1], 'k--', lw=1)
plt.plot([ws, ws], [-delta_2, delta_2], 'k--', lw=1)
plt.plot([np.pi, np.pi], [-delta_2, delta_2], 'k--', lw=1)

plt.plot(w[alt_idx], H[alt_idx], **alt_marker_style, zorder=15)
plt.plot(w[0], H[0], **alt_marker_style, zorder=15)

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
plt.text(wp, xtm, '$\omega_p$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ws - 0.02, xtm, '$\omega_s$', fontsize=fontsize, ha='right', va='baseline')
plt.text(np.pi + 0.04, xtm, '$\pi$', fontsize=fontsize, ha='left', va='baseline')
# eje y
plt.plot([0, ytl], [1, 1], 'k-', lw=1)
plt.plot([0, ytl], [delta_2, delta_2], 'k-', lw=1)
plt.plot([0, ytl], [-delta_2, -delta_2], 'k-', lw=1)
# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')

plt.text(0.1, ymax_ax, '$A_e(e^{j\omega})$', fontsize=fontsize, ha='left', va='center')

plt.axis('off')

#############
# Filtro 4: #
#############

wp = 0.3 * np.pi
ws = 0.43 * np.pi

wb = [0, wp, ws, np.pi]

h = signal.remez(M, wb, [1, 0], [1 / K, 1], fs=2*np.pi)

# Cálculo de la respuesta en frecuencia
_, H = signal.freqz(h, 1, worN=w, whole=False, plot=None, fs=2*np.pi)
# eliminación del componente de fase lineal
H = np.real(H * np.exp(1j * w * (M - 1) / 2))

# cálculo del error de aproximación para obtener los indices de las frecuencias de alternancia
# respuesta deseada
Hd = np.zeros(H.shape)
Hd[w <= wp] = 1
# función de ponderación
W = np.ones(H.shape)
W[w <= wp] = 1 / K
# error de aproximación
E = W * (Hd - H)
E[np.logical_and(w > wp, w < ws)] = 0
# indices de las frecuencias de alternancia (maximos del valor absoluto del error)
alt_idx = argrelextrema(np.abs(E), np.greater)[0]

# cálculo de delta para el filtro obtenido
delta = np.max(np.abs(H[w >= ws]))
# se redefinen delta 1 y delta 2
delta_2 = delta
delta_1 = K * delta

###################################################################

ax = plt.subplot2grid((4, 4), (2, 2), rowspan=2, colspan=2)
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
plt.plot([wp, wp], [0, 1 - delta_1], 'k--', lw=1)
plt.plot([ws, ws], [-delta_2, delta_2], 'k--', lw=1)
plt.plot([np.pi, np.pi], [-delta_2, delta_2], 'k--', lw=1)

plt.plot(w[alt_idx], H[alt_idx], **alt_marker_style, zorder=15)
plt.plot(w[0], H[0], **alt_marker_style, zorder=15)
plt.plot(w[-1], H[-1], **alt_marker_style, zorder=15)

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
plt.text(wp, xtm, '$\omega_p$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ws - 0.02, xtm, '$\omega_s$', fontsize=fontsize, ha='right', va='baseline')
plt.text(np.pi + 0.04, xtm, '$\pi$', fontsize=fontsize, ha='left', va='baseline')
# eje y
plt.plot([0, ytl], [1, 1], 'k-', lw=1)
plt.plot([0, ytl], [delta_2, delta_2], 'k-', lw=1)
plt.plot([0, ytl], [-delta_2, -delta_2], 'k-', lw=1)
# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')

plt.text(0.1, ymax_ax, '$A_e(e^{j\omega})$', fontsize=fontsize, ha='left', va='center')

plt.axis('off')

plt.savefig('filter_design_fir_optimal_pm_alternation_theorem_alternations.pdf', bbox_inches='tight')

plt.show()
