import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Diseño de filtro FIR óptimo multibanda

wb = [0, 0.3 * np.pi, 0.35 * np.pi, 0.6 * np.pi, 0.7 * np.pi, np.pi]
Kd = [1, 1, 0.2]
Gd = [0, 1, 0]
M = 74

h = signal.remez(M + 1, wb, Gd, Kd, fs=2*np.pi)

# Cálculo de la respuesta en frecuencia
nw = 8 * 1024
w = np.linspace(0, np.pi, nw)
_, H = signal.freqz(h, 1, worN=w, whole=False, plot=None, fs=2*np.pi)
# eliminación del componente de fase lineal
H = np.real(H * np.exp(1j * w * M / 2))

# calculo del error de aproximación
# respuesta deseada
Hd = np.zeros(H.shape)
Hd[np.logical_and(w >= wb[2], w <= wb[3])] = 1
# función de ponderación
W = np.ones(H.shape)
W[w <= wb[4]] = Kd[2]
# error de aproximación
E = Hd - H
E[np.logical_and(w > wb[1], w < wb[2])] = 0
E[np.logical_and(w > wb[3], w < wb[4])] = 0
delta = np.max(np.abs(E[w <= wb[3]]))
# indices de las frecuencias de alternancia (máximos del valor absoluto del error)
alt_idx = signal.argrelextrema(np.abs(E), np.greater)[0]
alt_idx = np.concatenate([[0], alt_idx, [nw - 1]])
E[np.logical_and(w > wb[1], w < wb[2])] = 'nan'
E[np.logical_and(w > wb[3], w < wb[4])] = 'nan'

# ripple en cada banda
deltas = delta / Kd
print('delta = {}'.format(delta))
print('deltas = {}'.format(deltas))

### Parámetros de la gráfica
fontsize = 10
fontsize2 = 11
# x y ticks labels margin
xtm = -0.07
ytm = -0.04
display_length = 6

xmin = 0
xmax = np.pi
xmax_ax = xmax
xmin_ax = xmin
ymin_ax = -90
ymax_ax = 10

grey = [0.9, 0.9, 0.9]
xticks = np.linspace(0, 1, 11)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[0] = '$0$'
xticks_labels[-1] = '$\pi$'
xticks *= np.pi

fig = plt.figure(0, figsize=(8, 5), frameon=False)
ax = plt.subplot2grid((11, 1), (0, 0), rowspan=6, colspan=1)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# magnitud de las respuestas en frecuencia
plt.plot(w, 20 * np.log10(np.abs(H)), 'k-', lw=1.5, zorder=10)
# mascara
plt.plot([wb[0], wb[1]], [20 * np.log10(Gd[0] + deltas[0]), 20 * np.log10(Gd[0] + deltas[0])], 'r--', lw=1)
plt.plot([wb[4], wb[5]], [20 * np.log10(Gd[2] + deltas[2]), 20 * np.log10(Gd[2] + deltas[2])], 'r--', lw=1)
plt.text(0.99, 0.9, '$|H(e^{j\omega})|$', fontsize=fontsize2, ha='right', va='baseline', transform=ax.transAxes)
plt.xlabel('$\omega$', fontsize=fontsize2)
plt.ylabel(r'$\textrm{dB}$', fontsize=fontsize2)
plt.xticks(xticks, xticks_labels, usetex=True)

dM = 0
xmax_ax = M + dM
xmin_ax = 0 - dM
ymax_ax = 0.6
ymin_ax = -0.2

n = np.arange(M + 1)

ax = plt.subplot2grid((11, 1), (7, 0), rowspan=5, colspan=1)
plt.xlim(xmin_ax, xmax_ax)
(markers, stemlines, bl) = plt.stem(n, h, linefmt='k', markerfmt='s', use_line_collection=True)
plt.setp(markers, markersize=3.5, markeredgecolor='k', markerfacecolor='k')
plt.setp(bl, visible=False)
plt.plot([xmin_ax, xmax_ax], [0, 0], 'k-', lw=1, zorder=-1)

plt.text(0.99, 0.85, '$h[n]$', fontsize=fontsize2, ha='right', va='baseline', transform=ax.transAxes)
plt.xlabel('$n$')
plt.ylabel(r'$\textrm{Amplitud}$', fontsize=fontsize2)

plt.savefig('filter_design_fir_optimal_pm_multiband_example_absH_h.pdf', bbox_inches='tight')

#
# gráfica del error de aproximación
#
xmax_ax = np.pi
xmin_ax = 0
f = 1.2
ymax_ax = np.max(deltas) * f
ymin_ax = -ymax_ax

alt_marker_style = dict(marker='o', linestyle='', markersize=5, markerfacecolor='w',
                        markeredgewidth=1.5, markeredgecolor='r')

fig = plt.figure(2, figsize=(8, 3), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.plot(w, E, 'k-', lw=1.5, label='$E_a(\omega)$')
plt.plot(w[alt_idx], E[alt_idx], **alt_marker_style, zorder=15, label=r'$\textrm{Alternancias}$')
plt.plot([wb[0], wb[1]], [deltas[0], deltas[0]], 'r--', lw=1, label=r'$1\pm K_i\delta$')
plt.plot([wb[0], wb[1]], [-deltas[0], -deltas[0]], 'r--', lw=1)
plt.plot([wb[2], wb[3]], [deltas[1], deltas[1]], 'r--', lw=1)
plt.plot([wb[2], wb[3]], [-deltas[1], -deltas[1]], 'r--', lw=1)
plt.plot([wb[4], wb[5]], [deltas[2], deltas[2]], 'r--', lw=1)
plt.plot([wb[4], wb[5]], [-deltas[2], -deltas[2]], 'r--', lw=1)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel('$\omega$', fontsize=fontsize2)
plt.ylabel('$E_A(\omega)$', fontsize=fontsize2)
plt.legend(loc='upper left', fontsize=fontsize2, frameon=False, framealpha=1)

plt.savefig('filter_design_fir_optimal_pm_multiband_example_aprox_error.pdf', bbox_inches='tight')

plt.show()
