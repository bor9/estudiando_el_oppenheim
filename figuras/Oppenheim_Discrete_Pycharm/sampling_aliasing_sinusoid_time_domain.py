import matplotlib.pyplot as plt
import numpy as np

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Parámetros
# intervalo temporal
T = 1
# número de muestras de la señal continua
Nd = 1000
# número de muestras de la señal discreta
N = 23
t = np.linspace(0, T, Nd)
n = np.linspace(0, T, N)
fs = (N - 1) / T

# frecuencias de las sinusoides
f1 = 0
f2 = 2
f3 = 9
f4 = 21

# sinusoides continuas
x1 = np.cos(2 * np.pi * f1 * t)
x2 = np.cos(2 * np.pi * f2 * t)
x3 = np.cos(2 * np.pi * f3 * t)
x4 = np.cos(2 * np.pi * f4 * t)

# sinusoides discretas
xn1 = np.cos(2 * np.pi * f1 * n)
xn2 = np.cos(2 * np.pi * f2 * n)
xn3 = np.cos(2 * np.pi * f3 * n)
xn4 = np.cos(2 * np.pi * f4 * n)

# gráficas

ymin = -1.5
ymax = 1.5
xmin = 0
xmax = T

ms = 3
y_label_coords = -0.15
fontsize = 11

fig = plt.figure(0, figsize=(8, 5), frameon=False)
ax = plt.subplot2grid((2, 4), (0, 0), rowspan=1, colspan=2)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.plot(t, x1, 'k-', lw=1)
plt.plot(n, xn1, ls='', marker='s', color='k', markersize=ms)
ax.xaxis.set_ticklabels([])
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$\Omega_0=0$')
plt.ylabel(r"$\textrm{Amplitud}$", fontsize=fontsize)

ax = plt.subplot2grid((2, 4), (0, 2), rowspan=1, colspan=2)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.plot(t, x2, 'k-', lw=1)
plt.plot(n, xn2, ls='', marker='s', color='k', markersize=ms)
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$\Omega_0\approx{:.2f}\Omega_s$'.format(f2 / fs))
#plt.title(r'$\textrm{FIR tipo II: }M\textrm{ impar, simetr\'ia par}$')

ax = plt.subplot2grid((2, 4), (1, 0), rowspan=1, colspan=2)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.plot(t, x3, 'k-', lw=1)
plt.plot(n, xn3, ls='', marker='s', color='k', markersize=ms)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$\Omega_0\approx{:.2f}\Omega_s$'.format(f3 / fs))
plt.xlabel(r"$\textrm{Tiempo (s)}$", fontsize=fontsize)
plt.ylabel(r"$\textrm{Amplitud}$", fontsize=fontsize)

ax = plt.subplot2grid((2, 4), (1, 2), rowspan=1, colspan=2)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.plot(t, x4, 'k-', lw=1)
plt.plot(n, xn4, ls='', marker='s', color='k', markersize=ms)
ax.yaxis.set_ticklabels([])
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$\Omega_0\approx{:.2f}\Omega_s$'.format(f4 / fs))
plt.xlabel(r"$\textrm{Tiempo (s)}$", fontsize=fontsize)

plt.savefig('sampling_aliasing_sinusoid_time_domain.pdf', bbox_inches='tight')

plt.show()
