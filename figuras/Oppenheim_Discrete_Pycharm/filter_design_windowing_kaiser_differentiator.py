import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from itertools import cycle

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


# parámetros de la ventana
Ms = [10, 9]
beta = 2.4

nw = 2 * 1024
w = np.linspace(0, np.pi, nw)

# frecuencia maxima hasta donde calcular el error
wm = 0.8 * np.pi

Nwins = len(Ms)
ns = []
hks = []
Hs = np.zeros((nw, Nwins))
Ediffs = np.zeros((nw, Nwins))

for i in np.arange(Nwins):
    M = Ms[i]
    # construcción de la respuesta al impulso
    wk = signal.windows.kaiser(M + 1, beta)
    n = np.arange(M + 1)
    alpha = M / 2
    hk = (np.cos(np.pi * (n - alpha)) / (n - alpha) - np.sin(np.pi * (n - alpha)) /
           (np.pi * ((n - alpha) ** 2))) * wk
    if M % 2 == 0:
        hk[M // 2] = 0
    # cálculo de la respuesta en frecuencia
    _, H = signal.freqz(hk, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
    # cálculo del error de aproximación
    # eliminación del componente de fase lineal
    Ao = np.real(H * np.exp(1j * (w * alpha - np.pi / 2)))
    Ediff = w - Ao
    Ediff[w >= wm] = 'nan'
    # se almacenan los datos
    ns.append(n)
    hks.append(hk)
    Hs[:, i] = np.abs(H)
    Ediffs[:, i] = Ediff


### Parámetros de la gráfica
fontsize = 10
fontsize2 = 11
# x y ticks labels margin
xtm = -0.07
ytm = -0.04
display_length = 6

# xticks
xticks = np.linspace(0, 1, 6)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[0] = '$0$'
xticks_labels[-1] = '$\pi$'
xticks *= np.pi

xmin = 0
xmax = np.pi
xmax_ax = xmax
xmin_ax = xmin
ymin_ax = 0
ymax_ax = np.pi

grey = 0.8 * np.ones((3, ))

leg = ['$|H(e^{{j\omega}})|\;(M={})$'.format(M) for M in Ms]
leg.append(r'$|H_\textrm{diff}(e^{{j\omega}})|\;(\textrm{ideal})$')

fig = plt.figure(0, figsize=(8, 5), frameon=False)
ax = plt.subplot2grid((16, 4), (0, 0), rowspan=9, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# magnitud de las respuestas en frecuencia
plt.plot(w, Hs, lw=1.5, zorder=2)
plt.plot([0, np.pi], [0, np.pi], lw=1.5, zorder=1)

plt.xticks(xticks, [], usetex=True)
plt.yticks(xticks, xticks_labels, usetex=True)
plt.legend(leg, loc='upper left', frameon=False, fontsize=fontsize2, framealpha=1)

#
# gráfica del error de aproximación
#
ymax_ax = 0.12
ymin_ax = -0.12
leg = ['$M={}$'.format(M) for M in Ms]

ax = plt.subplot2grid((16, 4), (10, 0), rowspan=6, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.plot(w, Ediffs, lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel('$\omega$', fontsize=fontsize2)
plt.ylabel(r'$E_\textrm{diff}(\omega)$', fontsize=fontsize2)
plt.legend(leg, loc='upper right', frameon=False, fontsize=fontsize2, framealpha=1)

plt.savefig('filter_design_windowing_kaiser_differentiator_approx_error.pdf', bbox_inches='tight')

#
# grafica de la respuesta al impulso
#
dM = 2
xmax_ax = max(Ms) + dM
xmin_ax = 0 - dM
ymax_ax = 1.5
ymin_ax = -1.5

# ciclco de colores
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = cycle(prop_cycle.by_key()['color'])
colors_srt = [next(colors) for i in np.arange(Nwins)]

fig = plt.figure(1, figsize=(8, 3), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
i = 0
(markers, stemlines, bl) = plt.stem(ns[i], hks[i], linefmt=colors_srt[i], markerfmt='s', use_line_collection=True)
plt.setp(markers, markersize=3.6, markeredgecolor=colors_srt[i], markerfacecolor=colors_srt[i])
plt.setp(bl, visible=False)
plt.plot([xmin_ax, xmax_ax], [0, 0], color=colors_srt[i], lw=1, zorder=-1)
plt.text(xmin_ax + 0.7, ymax_ax - 0.2, '$M={}$'.format(Ms[i]), fontsize=fontsize2, ha='left', va='top')

plt.xlabel('$n$')

ax = plt.subplot2grid((4, 4), (0, 2), rowspan=4, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
i = 1
(markers, stemlines, bl) = plt.stem(ns[i], hks[i], linefmt=colors_srt[i], markerfmt='s', use_line_collection=True)
plt.setp(markers, markersize=3.6, markeredgecolor=colors_srt[i], markerfacecolor=colors_srt[i])
plt.setp(bl, visible=False)
plt.plot([xmin_ax, xmax_ax], [0, 0], color=colors_srt[i], lw=1, zorder=-1)
ax.set_yticklabels([])
plt.text(xmin_ax + 0.7, ymax_ax - 0.2, '$M={}$'.format(Ms[i]), fontsize=fontsize2, ha='left', va='top')

plt.xlabel('$n$')

plt.savefig('filter_design_windowing_kaiser_differentiator_impulse_response.pdf', bbox_inches='tight')


plt.show()

