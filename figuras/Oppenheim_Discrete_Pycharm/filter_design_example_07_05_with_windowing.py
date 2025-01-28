import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from matplotlib.patches import Polygon
import math

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
wp = 0.4 * np.pi
ws = 0.6 * np.pi
delta_p1 = 0.01
delta_p2 = 0.01
delta_s = 0.001

# Parámetros del filtro FIR
# ripple permitido
delta = np.min([delta_p1, delta_p2, delta_s])
delta_dB = 20 * np.log10(delta)
# frecuencia de corte
wc = (wp + ws) / 2
# ancho de la banda de transición
Dw = ws - wp

print('delta_dB = {:.6f}'.format(delta_dB))
print('Dw = {:.6f}'.format(Dw))

# delta_dB = -60 dB. Se necesita ventana de Blackman.
# Para ventana de Blackman: Dw = 12 * pi / M
Mb1 = 12 * np.pi / Dw
print('M = {:.6f}'.format(Mb1))
Mb1 = math.ceil(Mb1)
# Otro tamaño de ventana menor que cumple las especificaciones.
Mb2 = Mb1 - 9
Mh = 60

# Diseño de filtro FIR por enventanado
# Construcción de la respuesta al impulso. Se incluye ventana de Hamming para comparar.
wb1 = signal.windows.get_window('blackman', Mb1, fftbins=False)
wb2 = signal.windows.get_window('blackman', Mb2, fftbins=False)
wh = signal.windows.get_window('hamming', Mh, fftbins=False)
n1 = np.arange(Mb1)
n2 = np.arange(Mb2)
n3 = np.arange(Mh)
alpha1 = (Mb1 - 1) / 2
alpha2 = (Mb2 - 1) / 2
alpha3 = (Mh - 1) / 2
hb1 = np.sin(wc * (n1 - alpha1)) / (np.pi * (n1 - alpha1)) * wb1
hb2 = np.sin(wc * (n2 - alpha2)) / (np.pi * (n2 - alpha2)) * wb2
hh = np.sin(wc * (n3 - alpha3)) / (np.pi * (n3 - alpha3)) * wh
if Mb1 % 2 != 0:
    hb1[(Mb1 - 1) // 2] = wc / np.pi
if Mb2 % 2 != 0:
    hb2[(Mb2 - 1) // 2] = wc / np.pi
if Mh % 2 != 0:
    hh[(Mh - 1) // 2] = wc / np.pi

# Cálculo de las respuestas en frecuencia
nw = 1024
w = np.linspace(0, np.pi, nw)
_, Hb1 = signal.freqz(hb1, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
_, Hb2 = signal.freqz(hb2, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
_, Hh = signal.freqz(hh, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)

### Parámetros de la gráfica
fontsize = 10
fontsize2 = 11
# x y ticks labels margin
xtm = -0.07
ytm = -0.04
display_length = 6

xmin = 0
xmax = np.pi
dx = 0.15
xmax_ax = xmax + dx
xmin_ax = xmin - dx
ymin_ax = -0.1
ymax_ax = 1.2

grey = [0.9, 0.9, 0.9]

fig = plt.figure(0, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hb1), label=r'$\textrm{{Blackman}}\;(M={})$'.format(Mb1), zorder=10)
plt.plot(w, np.absolute(Hb2), label=r'$\textrm{{Blackman}}\;(M={})$'.format(Mb2))
plt.plot(w, np.absolute(Hh), label=r'$\textrm{{Hamming}}\;(M={})$'.format(Mh))
# mascara
plt.plot([0, ws], [1 + delta_p1, 1 + delta_p1], 'k-')
plt.plot([0, wp], [1 - delta_p2, 1 - delta_p2], 'k-')
plt.plot([ws, np.pi], [delta_s, delta_s], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p2], 'k-')
plt.plot([ws, ws], [delta_s, 1 + delta_p1], 'k-')
# región pintada
vert = np.vstack(([0, wp, wp, 0], [0, 0, 1 - delta_p2, 1 - delta_p2]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [1 + delta_p1, 1 + delta_p1, mask_max, mask_max]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [delta_s, delta_s, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)

# etiquetas
# eje x
xtm2 = -0.11
plt.plot([wp, wp], [0, xtl], 'k-', lw=1)
plt.text(wp, xtm, '${:.1f}\pi$'.format(wp / np.pi), fontsize=fontsize, ha='center', va='baseline')
plt.plot([ws, ws], [0, xtl], 'k-', lw=1)
plt.text(ws, xtm, '${:.1f}\pi$'.format(ws / np.pi), fontsize=fontsize, ha='center', va='baseline')
plt.plot([np.pi, np.pi], [0, xtl], 'k-', lw=1)
plt.text(np.pi, xtm, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm, xtm, '$0$', fontsize=fontsize, ha='right', va='baseline')
# eje y
plt.plot([0, xtl], [1 - delta_p2, 1 - delta_p2], 'k-', lw=1)
plt.text(ytm, 1 - delta_p2, '${:.2f}$'.format(1 - delta_p2), fontsize=fontsize, ha='right', va='top')
plt.plot([0, xtl], [1 + delta_p1, 1 + delta_p1], 'k-', lw=1)
plt.text(ytm, 1 + delta_p1, '${:.2f}$'.format(1 + delta_p1), fontsize=fontsize, ha='right', va='bottom')
plt.plot([0, xtl], [delta_s, delta_s], 'k-', lw=1)
plt.text(ytm, delta_s, '${:.3f}$'.format(delta_s), fontsize=fontsize, ha='right', va='bottom')
# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize2, ha='center', va='baseline')
plt.text(0.06, ymax_ax, '$|H(e^{j\omega})|$', fontsize=fontsize2, ha='left', va='center')

plt.legend(bbox_to_anchor=(0.035, 0.3), loc='lower left', frameon=False, fontsize=fontsize2, framealpha=1)
plt.axis('off')
plt.savefig('filter_design_example_07_05_with_windowing.pdf', bbox_inches='tight')

###############################################

# x y ticks labels margin
display_length = 6

xmin = 0
dw = 0.15
xmax = wp + dw
xmax_ax = xmax
xmin_ax = xmin
ymin_ax = 0.96
ymax_ax = 1.02

xticks = np.arange(0, xmax_ax / np.pi, 0.1)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[0] = '$0$'
xticks *= np.pi

fig = plt.figure(1, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hb1), label=r'$\textrm{{Blackman}}\;(M={})$'.format(Mb1), zorder=10)
plt.plot(w, np.absolute(Hb2), label=r'$\textrm{{Blackman}}\;(M={})$'.format(Mb2))
plt.plot(w, np.absolute(Hh), label=r'$\textrm{{Hamming}}\;(M={})$'.format(Mh))
# mascara
plt.plot([0, ws], [1 + delta_p1, 1 + delta_p1], 'k-')
plt.plot([0, wp], [1 - delta_p2, 1 - delta_p2], 'k-')
plt.plot([ws, np.pi], [delta_s, delta_s], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p2], 'k-')
plt.plot([ws, ws], [delta_s, 1 + delta_p1], 'k-')
# región pintada
vert = np.vstack(([0, wp, wp, 0], [0, 0, 1 - delta_p2, 1 - delta_p2]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [1 + delta_p1, 1 + delta_p1, mask_max, mask_max]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [delta_s, delta_s, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)

# etiquetas
plt.xticks(xticks, xticks_labels, usetex=True)
# etiquetas de los ejes
plt.xlabel('$\omega$', fontsize=fontsize2)
plt.ylabel('$|H(e^{j\omega})|$', fontsize=fontsize2)
plt.legend(bbox_to_anchor=(0.035, 0.06), loc='lower left', frameon=False, fontsize=fontsize2, framealpha=1)

xmin = ws - dw
xmax = np.pi
xmax_ax = xmax
xmin_ax = xmin
ymin_ax = 0
ymax_ax = 0.003

xticks = np.arange(round(ws / np.pi, 1), 1.01, 0.1)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[-1] = '$\pi$'
xticks *= np.pi
yticks = np.linspace(ymin_ax, ymax_ax, int(ymax_ax/0.001 + 1))
yticks_labels = ['${:.3f}$'.format(yt) for yt in yticks]

ax = plt.subplot2grid((4, 4), (0, 2), rowspan=4, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hb1), zorder=10)
plt.plot(w, np.absolute(Hb2))
plt.plot(w, np.absolute(Hh))
# mascara
plt.plot([0, ws], [1 + delta_p1, 1 + delta_p1], 'k-')
plt.plot([0, wp], [1 - delta_p2, 1 - delta_p2], 'k-')
plt.plot([ws, np.pi], [delta_s, delta_s], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p2], 'k-')
plt.plot([ws, ws], [delta_s, 1 + delta_p1], 'k-')
# región pintada
vert = np.vstack(([0, wp, wp, 0], [0, 0, 1 - delta_p2, 1 - delta_p2]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [1 + delta_p1, 1 + delta_p1, mask_max, mask_max]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [delta_s, delta_s, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)

# etiquetas
plt.xticks(xticks, xticks_labels, usetex=True)
plt.yticks(yticks, yticks_labels, usetex=True)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
# etiquetas de los ejes
plt.xlabel('$\omega$', fontsize=fontsize2)
plt.ylabel('$|H(e^{j\omega})|$', fontsize=fontsize2)

plt.savefig('filter_design_example_07_05_with_windowing_zoom.pdf', bbox_inches='tight')
plt.show()

