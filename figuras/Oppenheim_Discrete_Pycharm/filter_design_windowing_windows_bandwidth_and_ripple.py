import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.optimize import curve_fit

__author__ = 'ernesto'

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


# función para ajuste mediante linea recta: y=f(x)
def f(x, A, B):
    return A * x + B

# parámetros
# frecuencia de corte del pasabajos
wc = np.pi / 2
# largo de las ventanas
Ms = np.arange(22, 100)
# espectros
nw = 8 * 1024
w = np.linspace(0, np.pi, nw)
# ventanas
win_names = ['boxcar', 'hann', 'hamming', 'blackman', 'blackmanharris', 'cosine']

Nwins = len(win_names)
NM = len(Ms)
Dw = np.zeros((NM, Nwins))
Pr = np.zeros((NM, Nwins))

for j, wname in enumerate(win_names):
    for i, M in enumerate(Ms):
        window = signal.windows.get_window(wname, M, fftbins=False)
        # espectro de la ventana para calcular el ancho del lóbulo principal
        _, W = signal.freqz(window, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
        # eliminación del componente de fase lineal
        alpha = (M - 1) / 2
        W = np.real(W * np.exp(1j * w * alpha))
        # indice del primer cruce por cero de W
        izc = np.argmax(W < 0)
        Dw[i, j] = 2 * w[izc - 1]

        # respuesta al impulso de pasabajos enventanado
        n = np.arange(M)
        h = np.sin(wc * (n - alpha)) / (np.pi * (n - alpha)) * window
        if M % 2 != 0:
            h[(M - 1) // 2] = wc / np.pi
        _, W = signal.freqz(h, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
        # eliminación del componente de fase lineal
        W = np.real(W * np.exp(1j * w * alpha))
        # indice del primer cruce por cero de W
        izc = np.argmax(W < 0)
        W = 20 * np.log10(np.abs(W))
        i_max = izc + np.argmax(W[izc:])
        Pr[i, j] = W[i_max]

# Ajuste mediante una recta de 2 * pi / Dw
Dw_adj = np.zeros((2, Nwins))
for i in np.arange(Nwins):
    popt, _= curve_fit(f, Ms, 2 * np.pi / Dw[:, i])
    Dw_adj[:, i] = popt

print('Dw/(2pi)={}'.format(1 / Dw_adj[0, :]))
print('delta (dB)={}'.format(np.median(Pr, axis=0)))

# Gráfica
leg = [r'$\textrm{{{}}}$'.format(win.capitalize()) for win in win_names]
# se modifican manualmente algunos nombres
leg[0] = r'$\textrm{Rectangular}$'
leg[-1] = r'$\textrm{Coseno}$'
leg[-2] = r'$\textrm{Blackman Harris}$'

fs = 11

fig = plt.figure(0, figsize=(7, 5), frameon=False)

# Gráfica de 2 * pi / Dw
# valores maximos y minimos de los ejes
xmax_ax = Ms[-1]
xmin_ax = Ms[0]
ymax_ax = 55
ymin_ax = 0

ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.plot(Ms, 2 * np.pi / Dw, ls='-', marker='s', markersize=2, lw=1)
for i in np.arange(Nwins):
    plt.plot(Ms, Ms * Dw_adj[0, i] + Dw_adj[1, i], 'k', lw=1, zorder=-1)
plt.xlabel('$M$', fontsize=fs)
plt.ylabel('$\dfrac{2\pi}{\Delta\omega}$', fontsize=fs)

# grafica de delta
ymax_ax = -10
ymin_ax = -120
plt.legend(leg, loc='upper left', frameon=False, framealpha=1, fontsize=fs)

ax = plt.subplot2grid((4, 4), (0, 2), rowspan=4, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.plot(Ms, Pr, ls='-', marker='s', markersize=2, lw=1)
plt.xlabel('$M$', fontsize=fs)
plt.ylabel('$20\log_{10}\delta$', fontsize=fs)
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")

plt.savefig('filter_design_windowing_windows_bandwidth_and_ripple.pdf', bbox_inches='tight')

plt.show()

