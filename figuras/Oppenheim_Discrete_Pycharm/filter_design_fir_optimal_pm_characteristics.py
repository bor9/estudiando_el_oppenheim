import matplotlib.pyplot as plt
import numpy as np
from scipy import signal, fft

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

#
# parámetros de los filtros pasabajos
#
# ancho de la banda de transición
dw = 0.2 * np.pi
# ponderacion de las bandas
K = 1
# Rango de parámetros M
Ms = np.arange(8, 16)
# Rango de frecuencias de corte
wps = np.linspace(0, np.pi - dw, 600)
# almacenamiento de valores de delta
nM = len(Ms)
nw = len(wps)
deltas = np.zeros((nw, nM))

# Cálculo de la respuesta en frecuencia
w = np.linspace(0, np.pi, 1024)

for j in np.arange(nM):
    taplenght = Ms[j] + 1
    for i in np.arange(nw):
        # diseño del filtro
        wp = wps[i]
        ws = wp + dw
        wb = [0, wp, ws, np.pi]
        h = signal.remez(taplenght, wb, [1, 0], [1 / K, 1], fs=2 * np.pi)
        # calculo de la respuesta en frecuencia
        _, H = signal.freqz(h, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
        # eliminación del componente de fase lineal
        H = np.real(H * np.exp(1j * w * (taplenght - 1) / 2))
        H[w <= wp] = 1 - H[w <= wp]
        H[np.logical_and(w > wp, w < ws)] = 0
        deltas[i, j] = np.max(np.abs(H))

#
# Grafica
#
fs = 11
leg = ['$M={}$'.format(Mj) for Mj in Ms]

xmin_ax = wps[0]
xmax_ax = wps[-1]
ymin_ax = 0
ymax_ax = 0.14

xticks = np.linspace(0, np.round((np.pi - dw) / np.pi, decimals=1), 5)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[0] = '$0$'
xticks *= np.pi

fig = plt.figure(0, figsize=(7, 5), frameon=False)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.plot(wps, deltas)
plt.legend(leg, fontsize=fs, ncol=4, frameon=False, framealpha=1, loc='center', bbox_to_anchor=(0.5, 0.05))
plt.xlabel(r'$\omega_p\textrm{ (rad)}$', fontsize=fs)
plt.ylabel(r'$\delta\textrm{ (ripple en la banda pasante o suprimida)}$', fontsize=fs)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.savefig('filter_design_fir_optimal_pm_characteristics.pdf', bbox_inches='tight')
plt.show()
