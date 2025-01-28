import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from itertools import cycle

__author__ = 'ernesto'

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# parámetros
# frecuencia de corte del pasabajos
wc = np.pi / 2
# largo de la ventana
M = 33
# respuesta al impulso del pasabajos ideal
n = np.arange(M)
alpha = (M - 1) / 2
h_lp = np.sin(wc * (n - alpha)) / (np.pi * (n - alpha))
if M % 2 != 0:
    h_lp[(M - 1) // 2] = wc / np.pi
# espectro
nw = 8 * 1024
w = np.linspace(0, np.pi, nw)
# ventanas
win_names = ['boxcar', 'cosine', 'hann', 'hamming', 'blackman', 'blackmanharris']

# cálculo del ancho de banda y el error de aproximación
Nwins = len(win_names)
Dw = np.zeros((Nwins, ))
Pr = np.zeros((Nwins, ))

for j, wname in enumerate(win_names):
    window = signal.windows.get_window(wname, M, fftbins=False)
    # respuesta al impulso de pasabajos enventanado
    h = h_lp * window
    _, W = signal.freqz(h, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
    # eliminación del componente de fase lineal
    W = np.real(W * np.exp(1j * w * alpha))
    ## cálculo de A
    # indice del primer cruce por cero de W
    izc = np.argmax(W < 0)
    W_dB = 20 * np.log10(np.abs(W))
    i_max = izc + np.argmax(W_dB[izc:])
    A = W_dB[i_max]
    # cálculo de Dw
    delta = 10 ** (A / 20)
    i1 = np.argmax(W < 1 - delta)
    i2 = np.argmax(W < delta)
    Dw[j] = w[i2] - w[i1]
    Pr[j] = A

# lo mismo para ventanas de Kaiser
betas = np.arange(1, 13)
# cálculo del ancho de banda y el error de aproximación
NwinsK = len(betas)
DwK = np.zeros((NwinsK, ))
PrK = np.zeros((NwinsK, ))

for j, beta in enumerate(betas):
    window = signal.windows.kaiser(M, beta, sym=True)
    # respuesta al impulso de pasabajos enventanado
    h = h_lp * window
    _, W = signal.freqz(h, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
    # eliminación del componente de fase lineal
    W = np.real(W * np.exp(1j * w * alpha))
    ## cálculo de A
    # indice del primer cruce por cero de W
    izc = np.argmax(W < 0)
    W_dB = 20 * np.log10(np.abs(W))
    i_max = izc + np.argmax(W_dB[izc:])
    A = W_dB[i_max]
    # cálculo de Dw
    delta = 10 ** (A / 20)
    i1 = np.argmax(W < 1 - delta)
    i2 = np.argmax(W < delta)
    DwK[j] = w[i2] - w[i1]
    PrK[j] = A


# Valores teóricos
Dws = np.linspace(0, np.pi / 2, 100)
As = -((M - 1) * 2.285 * Dws + 7.95)

# verificación del cálculo del ancho de banda
fig = plt.figure(0, figsize=(7, 5), frameon=False)
plt.plot(w, W)
plt.plot(w[i1], W[i1], 'r.')
plt.plot(w[i2], W[i2], 'r.')

# Gráfica
leg = [r'$\textrm{{{}}}$'.format(win.capitalize()) for win in win_names]
# se modifican manualmente algunos nombres
leg[0] = r'$\textrm{Rectangular}$'
leg[1] = r'$\textrm{Coseno}$'
leg[-1] = r'$\textrm{Blackman Harris}$'

# ciclco de colores
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = cycle(prop_cycle.by_key()['color'])
color1 = next(colors)
color2 = next(colors)

fs = 11
fs2 = 10

fig = plt.figure(1, figsize=(7, 4), frameon=False)

# Gráfica Dw vs A

# valores maximos y minimos de los ejes
xmax_ax = np.pi / 2
xmin_ax = 0
ymax_ax = -15
ymin_ax = -130
xticks = np.linspace(0.0, 0.5, 6)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[0] = '$0$'
xticks *= np.pi

ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.plot(DwK, PrK, color=color1, ls='', marker='s', markersize=5)
plt.plot(Dw, Pr, color=color2, ls='', marker='o', markersize=5)
plt.plot(Dws, As, 'k--', lw=1, zorder=0)
for i in np.arange(NwinsK):
    plt.text(DwK[i]-0.02, PrK[i], r'$\textrm{{Kaiser}}({})$'.format(betas[i]), fontsize=fs2, ha='right', va='top')
for i in np.arange(Nwins - 1):
    plt.text(Dw[i] + 0.02, Pr[i], leg[i], fontsize=fs2, ha='left', va='bottom')
i += 1
plt.text(Dw[i] + 0.04, Pr[i] + 7, r'$\textrm{Blackman}$', fontsize=fs2, ha='right', va='bottom')
plt.text(Dw[i] + 0.04, Pr[i] + 2, r'$\textrm{Harris}$', fontsize=fs2, ha='right', va='bottom')

plt.xlabel('$\Delta\omega$', fontsize=fs)
plt.ylabel('$20\log_{10}\delta$', fontsize=fs)

plt.xticks(xticks, xticks_labels, usetex=True)

plt.savefig('filter_design_windowing_kaiser_bandwidth_and_ripple.pdf', bbox_inches='tight')

plt.show()

