import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from itertools import cycle
from scipy.special import i0

__author__ = 'ernesto'

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

#
# Ventanas de Kaiser variando beta
#
# largo de las ventanas
M = 25 # este valor es M+1 en Oppenheim

n = np.arange(M)
# construcción de las ventanas
# matrices con las ventanas en las columnas
betas = [0, 3, 6, 9]
Nwins = len(betas)
wins = np.zeros((M, Nwins))

# espectros
nw = 1024
w = np.linspace(0, np.pi, nw)
Ws = np.zeros((nw, Nwins))
for i in np.arange(Nwins):
    wins[:, i] = signal.windows.kaiser(M, betas[i])
    _, W = signal.freqz(wins[:, i], 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
    Ws[:, i] = 20 * np.log10(np.abs(W))
# para ganancia de 0dB en omega=0
Ws = Ws - Ws[0, :]

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

for i in np.arange(Nwins):

    ax = plt.subplot2grid((Nwins, 11), (i, 0), rowspan=1, colspan=5)

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
ymin_ax = -140
xticks = [0, np.pi / 2, np.pi]
xticks_labels = ['$0$', '$\dfrac{\pi}{2}$', '$\pi$']
leg = [r'$\beta={}$'.format(beta) for beta in betas]

ax = plt.subplot2grid((Nwins, 11), (0, 6), rowspan=Nwins, colspan=5)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.plot(w, Ws, lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel('$\omega$')
plt.ylabel('$20\log_{10}|W(e^{j\omega})|$')
ax.yaxis.set_label_coords(-0.14, 0.5)
plt.legend(leg, loc='upper right', frameon=False, framealpha=1)

plt.savefig('filter_design_windowing_kaiser_windows_and_spectrums_beta.pdf', bbox_inches='tight')

#
# Ventanas de Kaiser variando M
#
beta = 6
Ms = [10, 25, 40]
# construcción de las ventanas
# matrices con las ventanas en las columnas
Nwins = len(Ms)
# espectros
nw = 2048
w = np.linspace(0, np.pi, nw)
Ws = np.zeros((nw, Nwins))
for i in np.arange(Nwins):
    win = signal.windows.kaiser(Ms[i], beta)
    _, W = signal.freqz(win, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
    Ws[:, i] = 20 * np.log10(np.abs(W))
# para ganancia de 0dB en omega=0
Ws = Ws - Ws[0, :]

# cálculo de la amplitud del primer lobulo secundario
# eliminación del componente de fase lineal paara encontrar el primer cruce por cero
alpha = (Ms[i] - 1) / 2
W = np.real(W * np.exp(1j * w * alpha))
# indice del primer cruce por cero de W
izc = np.argmax(W < 0)
# amplitud del primer lobulo secundario
# es el valor maximo del espectro luegodel primer curce por cero
lobAmp = np.max(Ws[izc:-1, i])

ymax_ax = 10
ymin_ax = -130
leg = [r'$M={}$'.format(M) for M in Ms]

fig = plt.figure(1, figsize=(7, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.plot(w, Ws, lw=1.5)
plt.plot([xmin_ax, xmax_ax], [lobAmp, lobAmp], ls='--', lw=1, zorder=-1, color = 0.7 * np.ones((3, )))
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel('$\omega$')
plt.ylabel('$20\log_{10}|W(e^{j\omega})|$')
plt.legend(leg, loc='upper right', frameon=False, framealpha=1)

plt.savefig('filter_design_windowing_kaiser_windows_and_spectrums_M.pdf', bbox_inches='tight')

plt.show()

# Prueba: construcción de una ventana de Kaiser manualmente
# a partir de la función modificada de Bessel de primer tipo.
M = 31
beta = 6
alpha = (M - 1) / 2
n = np.arange(M)
wm = i0(beta * np.sqrt(1 - ((n - alpha) / alpha) ** 2)) / i0(beta)
w = signal.windows.kaiser(M, beta)

print('sum |wm - w| = {}'.format(np.sum(np.abs(wm-w))))

