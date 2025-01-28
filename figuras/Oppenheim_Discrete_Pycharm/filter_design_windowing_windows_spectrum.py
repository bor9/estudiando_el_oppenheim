import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from itertools import cycle

__author__ = 'ernesto'

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# parámetros
# largo de las ventanas
M = 31 # este valor es M+1 en Oppenheim

n = np.arange(M)

# construcción de las ventanas
# matriz con las ventnas en las columnas
Nwins = 5
wins = np.zeros((M, Nwins))
wins[:, 0] = signal.windows.boxcar(M)
wins[:, 1] = signal.windows.bartlett(M)
wins[:, 2] = signal.windows.hann(M)
wins[:, 3] = signal.windows.hamming(M)
wins[:, 4] = signal.windows.blackman(M)

# espectros
nw = 1024
w = np.linspace(0, np.pi, nw)
Ws = np.zeros((nw, Nwins))
for i in np.arange(Nwins):
    _, W = signal.freqz(wins[:, i], 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
    Ws[:, i] = 20 * np.log10(np.abs(W))
# para ganacia 0 dB en omega = 0
Ws -= Ws[0, :]

# Gráfica

# ciclco de colores
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = cycle(prop_cycle.by_key()['color'])
colors_srt = [next(colors) for i in np.arange(Nwins)]

# graficas de las ventanas en el tiempo
# valores maximos y minimos de los ejes
dM = 4
xmax_ax = M + dM
xmin_ax = 0 - dM
ymax_ax = 1.2
ymin_ax = -0.2

fig = plt.figure(0, figsize=(8, 5), frameon=False)

# ventana rectangular
for i in np.arange(Nwins):

    ax = plt.subplot2grid((5, 11), (i, 0), rowspan=1, colspan=5)

    plt.xlim(xmin_ax, xmax_ax)
    plt.ylim(ymin_ax, ymax_ax)

    (markers, stemlines, bl) = plt.stem(n, wins[:, i], linefmt=colors_srt[i], markerfmt='s', use_line_collection=True)
    plt.setp(markers, markersize=3, markeredgecolor=colors_srt[i], markerfacecolor=colors_srt[i])
    plt.setp(bl, visible=False)
    plt.plot([xmin_ax, xmax_ax], [0, 0], 'k-', lw=1, zorder=-1)
    if i != Nwins-1:
        ax.xaxis.set_ticklabels([])

plt.xlabel('$n$')

# graficas de las ventanas en la frecuencia
xmax_ax = np.pi
xmin_ax = 0
ymax_ax = 10
ymin_ax = -120
xticks = [0, np.pi / 2, np.pi]
xticks_labels = ['$0$', '$\dfrac{\pi}{2}$', '$\pi$']
leg = [r'$\textrm{Rectangular}$', r'$\textrm{Bartlett}$', r'$\textrm{Hann}$', r'$\textrm{Hamming}$',
       r'$\textrm{Blackman}$']

ax = plt.subplot2grid((5, 11), (0, 6), rowspan=5, colspan=5)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.plot(w, Ws, lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel('$\omega$')
plt.ylabel('$20\log_{10}|W(e^{j\omega})|$')
ax.yaxis.set_label_coords(-0.14, 0.5)
plt.legend(leg, loc='lower left', bbox_to_anchor=(0.5, 0.7), frameon=False, framealpha=1)

plt.savefig('filter_design_windowing_windows_spectrum.pdf', bbox_inches='tight')

plt.show()





