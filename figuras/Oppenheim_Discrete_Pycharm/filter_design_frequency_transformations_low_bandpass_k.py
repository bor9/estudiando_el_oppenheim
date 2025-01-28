import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import numpy.polynomial.polynomial as poly

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


def filter_transformation(b, a, gnum, gden):
    """
    Transformación en frecuencia de un filtro pasabajos.
    Cálculo de los coefcientes del numerador y denominador del filtro
    transformado a partir de los coeficientes del numerador y denominador
    del filtro pasabajos prototipo y la transformación.

    Parámetros
    ----------
    b, a : vectores de dimensión (N + 1, )
        Coeficientes del numerador y denominador del filtro prototipo
    gnum, gden : vectores de una dimensión de tamaño arbitrario
        Coeficientes del numerador y denominador de la transformación

    Retorna
    -------
    bt, at : vectores de una dimensión
        Coeficientes del numerador y denominador del filtro transformado
    """
    # Orden del filtro a transformar
    N = len(b) - 1
    # Lista de polinomios que multiplican a los coeficientes b[k], a[k].
    # Son de la forma (gden ** (N - k)) * (gnum ** k)
    pfactors = [1] * (N + 1)
    k = N
    for i in np.arange(N + 1):
        for j in np.arange(k):
            pfactors[i] = poly.polymul(pfactors[i], gden)
        for j in np.arange(N - k):
            pfactors[i] = poly.polymul(pfactors[i], gnum)
        k -= 1
    # Cálculo de los coeficientes del filtro transformado
    bt = 0
    at = 0
    for i in np.arange(N + 1):
        bt = poly.polyadd(bt, b[i] * pfactors[i])
        at = poly.polyadd(at, a[i] * pfactors[i])
    return bt, at


# pasabajos a pasabajos
# frecuencia original
theta_p = np.pi / 2
# parámetros del pasabanda
# alpha = 0 (w1 + w2 = pi)
w_sum = np.pi
# frecuencia w1
dw1 = np. pi / 8
w1 = np.arange(dw1, np.pi / 2, dw1)
w2 = w_sum - w1

# k es fijo por ser el ancho de banda fijo
ks = np.tan(theta_p / 2) / np.tan((w2 - w1) / 2)
alphas = np.cos((w2 + w1) / 2) / np.cos((w2 - w1) / 2)

gamma1 = 2 * alphas * ks / (ks + 1)
gamma2 = (ks - 1) / (ks + 1)

# Relación entre las frecuencias
Nfrecs = 256
w = np.linspace(0, np.pi, Nfrecs)

thetas = -np.arctan2(2 * ks * (np.cos(w[:, np.newaxis]) - alphas) * np.sin(w[:, np.newaxis]),
                     -(ks * (np.cos(w[:, np.newaxis]) - alphas)) ** 2 + np.sin(w[:, np.newaxis]) ** 2)

#thetas = -np.pi - 2 * np.arctan2(-np.sin(w[:, np.newaxis]), ks * (np.cos(w[:, np.newaxis]) - alphas))

# Transformación de un filtro pasabajos
delta_p = 0.1
delta_s = 0.05
tp = theta_p
ts = 0.6 * np.pi

gpass = -20 * np.log10(1 - delta_p)
gstop = -20 * np.log10(delta_s)

# Filtro Butterworth
N, Wn = signal.buttord(tp, ts, gpass, gstop, analog=False, fs=2*np.pi)
b, a = signal.butter(N, Wn, btype='low', analog=False, output='ba', fs=2*np.pi)
theta, H = signal.freqz(b, a, worN=Nfrecs, whole=False, plot=None, fs=2*np.pi, include_nyquist=False)

# Filtros transformados
Hts = np.empty((Nfrecs, len(ks)))
for i, k in enumerate(ks):
    tnum = [-gamma2[i], gamma1[i], -1]
    tden = [1, -gamma1[i], gamma2[i]]
    bt, at = filter_transformation(b, a, tnum, tden)
    _, Ht = signal.freqz(bt, at, worN=Nfrecs, whole=False, plot=None, fs=2 * np.pi, include_nyquist=False)
    Hts[:, i] = np.abs(Ht)

# Gráfica
# Simetrización del filtro pasabajos
theta_sym = np.concatenate((-np.flip(theta[1:]), theta))
H = np.abs(H)
H_sym = np.concatenate((np.flip(H[1:]), H))


thmin = -np.pi
thmax = np.pi
wmin = 0
wmax = np.pi
gmin = 0
gmax = 1.1

fs = 11  # fontsize
grey = 0.7 * np.ones((3, ))
markersize = 5

thticks = np.linspace(thmin, thmax, 5)
thticks_labels = ['$-\pi$', '$-\dfrac{\pi}{2}$', '$0$', '$\dfrac{\pi}{2}$', '$\pi$']
wticks = np.linspace(wmin, wmax, 9)
wticks_labels = ['$0$', '$\dfrac{\pi}{8}$', '$\dfrac{\pi}{4}$', '$\dfrac{3\pi}{8}$', '$\dfrac{\pi}{2}$',
                 '$\dfrac{5\pi}{8}$', '$\dfrac{3\pi}{4}$', '$\dfrac{7\pi}{8}$', '$\pi$']
gticks = np.linspace(0, 1, int(1/0.2) + 1, endpoint=True)

leg = []
for i in np.arange(len(w1)):
    if np.abs(ks[i] - round(ks[i])) > 1e-10:
        s = r'$k\approx{:.3f}$'.format(ks[i])
    else:
        s = r'$k={:d}$'.format(int(np.round(ks[i])))
    s += r'$\;\Big(\omega_1=$'
    s += wticks_labels[i + 1]
    s += r'$\;\Big)$'
    leg.append(s)

fig = plt.figure(0, figsize=(8, 8), frameon=False)

ax1 = plt.subplot2grid((20, 20), (0, 0), rowspan=11, colspan=7)
plt.xlim(gmin, gmax)
plt.ylim(thmin, thmax)
plt.gca().invert_xaxis()
plt.plot(H_sym, theta_sym, 'k')
plt.plot([1 - delta_p, 1 - delta_p], [thmin, theta_p], '--', lw=1, color=grey, zorder=-1)
plt.plot([0, 1 - delta_p], [theta_p, theta_p], '--', lw=1, color=grey, zorder=-1)
plt.plot([0, 1 - delta_p], [-theta_p, -theta_p], '--', lw=1, color=grey, zorder=-1)
plt.plot([1 - delta_p, 1 - delta_p], [-theta_p, theta_p], marker='o', markersize=markersize, color='k', linestyle='')
plt.xticks(gticks, usetex=True)
plt.yticks(thticks, thticks_labels, usetex=True)
plt.xlabel(r'$|H_\textrm{lp}(e^{j\theta})|$', fontsize=fs)
plt.text(-0.175, 0, r'$\theta$', fontsize=fs, ha='center', va='center')
ax1.yaxis.tick_right()

ax2 = plt.subplot2grid((20, 20), (0, 9), rowspan=11, colspan=10)
plt.ylim(thmin, thmax)
plt.xlim(wmin, wmax)
plt.plot(w, thetas)
plt.gca().set_prop_cycle(None)
plt.plot([w1, w1], [theta_sym[0], -theta_p], '--', lw=1, zorder=-1)
plt.gca().set_prop_cycle(None)
plt.plot([w1], [-theta_p], marker='o', markersize=markersize)
plt.gca().set_prop_cycle(None)
plt.plot([w2, w2], [theta_sym[0], theta_p], '--', lw=1, zorder=-2)
plt.gca().set_prop_cycle(None)
plt.plot([w2], [theta_p], marker='o', markersize=markersize)
plt.plot([0, w2[0]], [theta_p, theta_p], '--', lw=1, color=grey, zorder=-3)
plt.plot([0, w1[-1]], [-theta_p, -theta_p], '--', lw=1, color=grey, zorder=-3)
plt.xticks(wticks, wticks_labels, usetex=True)
plt.yticks(thticks, thticks_labels, usetex=True)
plt.xlabel(r'$\omega$', fontsize=fs)

ax3 = plt.subplot2grid((20, 20), (13, 9), rowspan=6, colspan=10)
plt.xlim(wmin, wmax)
plt.ylim(gmin, gmax)
plt.plot(theta, Hts)
plt.gca().set_prop_cycle(None)
plt.plot([w1, w1], [0, 1 - delta_p], '--', lw=1)
plt.gca().set_prop_cycle(None)
plt.plot([w1], [1 - delta_p], marker='o', markersize=markersize)
plt.gca().set_prop_cycle(None)
plt.plot([w2, w2], [0, 1 - delta_p], '--', lw=1)
plt.gca().set_prop_cycle(None)
plt.plot([w2], [1 - delta_p], marker='o', markersize=markersize)
plt.plot([0, w2[0]], [1 - delta_p, 1 - delta_p], '--', lw=1, color=grey, zorder=-1)
plt.xticks(wticks, wticks_labels, usetex=True)
plt.yticks(gticks, usetex=True)
plt.xlabel(r'$\omega$', fontsize=fs)
plt.ylabel(r'$|H(e^{j\omega})|$', fontsize=fs)
plt.legend(leg, loc='center right', bbox_to_anchor=(-0.25, 0.5), frameon=False, framealpha=1)

label = r'$\alpha=0\;(\omega_2+\omega_1=\pi)$'
plt.text(-2.65, 0.2, label, fontsize=10, ha='left', va='top')
plt.text(-2.65, 0.05, r'$\theta_p=\dfrac{\pi}{2}$', fontsize=10, ha='left', va='top')

plt.savefig('filter_design_frequency_transformations_low_bandpass_k.pdf', bbox_inches='tight')

plt.show()
