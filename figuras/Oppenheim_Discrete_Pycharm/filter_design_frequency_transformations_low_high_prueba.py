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


nfrecs = 200
theta = np.linspace(0, np.pi, nfrecs)

# pasabajos a pasabajos
# frecuencia original
theta_p = np.pi / 2
# frecuencias transformadas
dw = np. pi / 8
w_p = np.arange(dw, np.pi, dw)

# alphas
# alphas = np.sin((theta_p - w_p) / 2) / np.sin((theta_p + w_p) / 2)
#ws = np.arctan2((1 - alphas ** 2) * np.sin(theta[:, np.newaxis]),
#                2 * alphas + (1 + alphas ** 2) * np.cos(theta[:, np.newaxis]))
alphas = -np.cos((theta_p - w_p) / 2) / np.cos((theta_p + w_p) / 2)
ws = -np.arctan2((1 - alphas ** 2) * np.sin(theta[:, np.newaxis]),
                -2 * alphas - (1 + alphas ** 2) * np.cos(theta[:, np.newaxis]))



## PRUEBA
#alphas = 1 / alphas
ws2 = -np.arctan2((alphas ** 2 - 1) * np.sin(theta[:, np.newaxis]),
                -2 * alphas - (1 + alphas ** 2) * np.cos(theta[:, np.newaxis]))

print(np.sum(np.abs(ws - ws2)))

alphas = -np.cos((theta_p + w_p) / 2) / np.cos((theta_p - w_p) / 2)
ws = np.arctan2((1 - alphas ** 2) * np.sin(theta[:, np.newaxis]),
                -2 * alphas - (1 + alphas ** 2) * np.cos(theta[:, np.newaxis]))
#ws = np.arctan((1 - alphas ** 2) * np.sin(theta[:, np.newaxis])/
#               (2 * alphas + (1 + alphas ** 2) * np.cos(theta[:, np.newaxis])))



### prueba para el cálculo de alpha
alphas_teo = -(np.exp(-1j * theta_p) + np.exp(-1j * w_p)) / (1 + np.exp(-1j * (theta_p + w_p)))
#alphas_teo = -(np.exp(-1j * theta_p) + np.exp(-1j * (np.pi - w_p))) / (1 + np.exp(-1j * (theta_p + (np.pi - w_p))))
print(np.real(alphas_teo))
print(alphas)

# Transformación de un filtro pasabajos
delta_p = 0.1
delta_s = 0.05
tp = theta_p
ts = 0.6 * np.pi

gpass = -20 * np.log10(1 - delta_p)
gstop = -20 * np.log10(delta_s)

# Filtro elíptico
Nfrec = 512
N, Wn = signal.buttord(tp, ts, gpass, gstop, analog=False, fs=2*np.pi)
b, a = signal.butter(N, Wn, btype='low', analog=False, output='ba', fs=2*np.pi)
w, H = signal.freqz(b, a, worN=Nfrec, whole=False, plot=None, fs=2*np.pi, include_nyquist=False)

print(N)

# Filtros transformados
Hts = np.empty((Nfrec, len(alphas)))
for i, alpha in enumerate(alphas):
    tnum = [-alpha, -1]
    tden = [1, alpha]
    bt, at = filter_transformation(b, a, tnum, tden)
    _, Ht = signal.freqz(bt, at, worN=Nfrec, whole=False, plot=None, fs=2 * np.pi, include_nyquist=False)
    Hts[:, i] = np.abs(Ht)


# Graficas
xmin = 0
xmax = np.pi
ymin = 0
ymax = np.pi

fs = 11  # fontsize
y_label_coords = -0.08
grey = 0.7 * np.ones((3, ))

xticks = np.linspace(xmin, xmax, 5)
xticks_labels = ['$0$', '$\dfrac{\pi}{4}$', '$\dfrac{\pi}{2}$', '$\dfrac{3\pi}{4}$', '$\pi$']
yticks = np.linspace(xmin, xmax, 9)
yticks_labels = ['$0$', '$\dfrac{\pi}{8}$', '$\dfrac{\pi}{4}$', '$\dfrac{3\pi}{8}$', '$\dfrac{\pi}{2}$',
                 '$\dfrac{5\pi}{8}$', '$\dfrac{3\pi}{4}$', '$\dfrac{7\pi}{8}$', '$\pi$']

leg = []
for i in np.arange(len(w_p)):
    if alphas[i] != np.floor(alphas[i]):
        s = r'$\alpha\approx{:.3f}$'.format(alphas[i])
    else:
        s = r'$\alpha={:d}$'.format(int(alphas[i]))
    s += r'$\;\Big(\omega_p=$'
    s += yticks_labels[i + 1]
    s += r'$\;\Big)$'
    leg.append(s)

fig = plt.figure(0, figsize=(7, 4), frameon=False)
#plt.ylim(ymin, ymax)
plt.xlim(xmin, xmax)
plt.gca().set_aspect('equal', adjustable='box')
plt.plot(theta, ws)
plt.gca().set_prop_cycle(None)
plt.plot([0, theta_p], [w_p, w_p], '--', lw=1)
plt.plot([theta_p, theta_p], [0, w_p[-1]], '--', lw=1, color=grey)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.yticks(yticks, yticks_labels, usetex=True)
plt.xlabel(r'$\theta$', fontsize=fs)
plt.ylabel(r'$\omega$', fontsize=fs)
plt.legend(leg, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, framealpha=1)

#plt.savefig('filter_design_frequency_transformations_low_high_mapping.pdf', bbox_inches='tight')

ymin = 0
ymax = 1.05

fig = plt.figure(1, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.plot(w, Hts)
plt.gca().set_prop_cycle(None)
plt.plot([w_p, w_p], [0, 1 - delta_p], '--', lw=1)
plt.plot([0, w_p[-1]], [1 - delta_p, 1 - delta_p], '--', lw=1, color=grey)
plt.xticks(yticks, yticks_labels, usetex=True)
plt.xlabel(r'$\omega$', fontsize=fs)
plt.ylabel(r'$|H(e^{j\omega})|$', fontsize=fs)
#plt.savefig('filter_design_frequency_transformations_low_high_magnitude_response.pdf', bbox_inches='tight')


## Prueba
idx = 0
tnum = [-alphas[idx], -1]
tden = [1, alphas[idx]]
bt, at = filter_transformation(b, a, tnum, tden)
z, p, k = signal.tf2zpk(bt, at)

xmax = 1.5
xmin = -xmax
ymax = 1.5
ymin = -1.2
ymax_ax = 2
ytm = 0.08
xtm = -0.24

fs2 = 10

# circulo unidad
theta = np.linspace(0, 2 * np.pi, 100)
xc = np.cos(theta)
yc = np.sin(theta)

zeros_marker_style = dict(marker='o', linestyle='', markersize=6, markerfacecolor='w', markeredgecolor='k',
                          markeredgewidth=1.5)
poles_marker_style = dict(marker='x', linestyle='', markersize=6, markeredgecolor='k', markeredgewidth=2)


fig = plt.figure(2, figsize=(5, 5), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.ylim(ymin, ymax_ax)
plt.xlim(xmin, xmax)
plt.gca().set_aspect('equal', adjustable='box')
# ejes
# axis arrows
plt.annotate("", xytext=(xmin, 0), xycoords='data', xy=(xmax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin), xycoords='data', xy=(0, ymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
# circulo unidad
ax.plot(xc, yc, 'k-', lw=1)
# polos y ceros
plt.plot(z.real, z.imag, **zeros_marker_style, zorder=10)
plt.plot(p.real, p.imag, **poles_marker_style)

# etiquetas de los ejes
#plt.text(xmax, xtm, r'$\textrm{Re}(z)$', fontsize=fontsize, ha='center', va='baseline')
#plt.text(ytm, ymax, r'$\textrm{Im}(z)$', fontsize=fontsize, ha='left', va='center')

#plt.text(0, 1.85, r'$\textrm{{El\'iptico (N={})}}$'.format(Ne),
#         fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

plt.show()

