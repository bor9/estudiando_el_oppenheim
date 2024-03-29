import matplotlib.pyplot as plt
import numpy as np

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Parámetros de los pasabajos
# frecuencia de corte
omega_c = 0.4 * np.pi
# retardo
alphas = np.array([5, 4.5, 4.3])
# rango de n
nmin = -5
nmax = 15
# muestras
n = np.arange(nmin, nmax)
# puntos teporales para las señales interpoladas
Ni = 500
ni = np.linspace(nmin, nmax, Ni)
# construcción de las señales
hh = np.zeros((len(alphas), len(n)))
hhi = np.zeros((len(alphas), len(ni)))
for i, alpha in enumerate(alphas):
    hh[i, :] = np.sin(omega_c * (n - alpha)) / (np.pi * (n - alpha))
np.nan_to_num(hh, copy=False, nan=omega_c / np.pi)
for i, alpha in enumerate(alphas):
    hhi[i, :] = np.sin(omega_c * (ni - alpha)) / (np.pi * (ni - alpha))
np.nan_to_num(hhi, copy=False, nan=omega_c / np.pi)

# limites del eje y
dy = 0.05
ymax = np.amax(hh) + dy
ymin = np.amin(hh) - dy

fs = 12
ms = 3

# Figuras
fig = plt.figure(0, figsize=(8.5, 3), frameon=False)

# retardo entero
ax = plt.subplot2grid((1, 6), (0, 0), rowspan=1, colspan=2)
plt.xlim(nmin, nmax)
plt.ylim(ymin, ymax)
i = 0
plt.plot(n, hh[i, :], ls='', marker='s', color='k', markersize=ms)
plt.plot(ni, hhi[i, :], 'k-', lw=1)
plt.xlabel(r"$n$", fontsize=fs)
plt.ylabel(r"$\textrm{Amplitud}$", fontsize=fs)
plt.text(14, 0.4, r'$\alpha={:.0f}$'.format(alphas[i]), fontsize=fs, ha='right', va='baseline')

# 2alpha entero
ax = plt.subplot2grid((1, 6), (0, 2), rowspan=1, colspan=2)
plt.xlim(nmin, nmax)
plt.ylim(ymin, ymax)
i = 1
plt.plot(n, hh[i, :], ls='', marker='s', color='k', markersize=ms)
plt.plot(ni, hhi[i, :], 'k-', lw=1)
plt.xlabel(r"$n$", fontsize=fs)
ax.yaxis.set_ticklabels([])
plt.text(14, 0.4, r'$\alpha={:.1f}$'.format(alphas[i]), fontsize=fs, ha='right', va='baseline')

# alpha arbitrario
ax = plt.subplot2grid((1, 6), (0, 4), rowspan=1, colspan=2)
plt.xlim(nmin, nmax)
plt.ylim(ymin, ymax)
i = 2
plt.plot(n, hh[i, :], ls='', marker='s', color='k', markersize=ms)
plt.plot(ni, hhi[i, :], 'k-', lw=1)
plt.xlabel(r"$n$", fontsize=fs)
ax.yaxis.set_ticklabels([])
plt.text(14, 0.4, r'$\alpha={:.1f}$'.format(alphas[i]), fontsize=fs, ha='right', va='baseline')

# save as pdf image
plt.savefig('transform_analysis_example_05_14_linear_phase.pdf', bbox_inches='tight')

plt.show()
