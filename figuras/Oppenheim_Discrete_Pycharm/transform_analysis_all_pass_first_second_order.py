import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from cycler import cycler

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# sistema de primer orden con cero en z=re^{j\theta}
# 1 - theta variable y r fijo
aa = [0.9, -0.9]
nw = 512

nt = len(aa)
# inicialización de matrices
magHs = np.zeros((nw, nt))
argHs = np.zeros((nw, nt))
grdHs = np.zeros((nw, nt))
for i in np.arange(nt):
    # vector del numerador y denominador de la función del sistema
    b = [-aa[i] , 1]
    a = [1, -aa[i]]
    # respuesta en frecuencia
    w, H = signal.freqz(b, a, nw, whole=True)
    # magnitud de la respuesta en frecuencia
    magHs[:, i] = 20 * np.log10(np.abs(H))
    # fase de la respuesta en frecuencia
    argHs[:, i] = np.angle(H)
    # retardo de grupo de la respuesta en frecuencia
    _, grdHs[:, i] = signal.group_delay((b, a), w)


########## Gráficas ##########

fs = 11  # fontsize
cy = cycler('color', ['black', 'red', 'blue', 'green'])
y_label_coords = -0.12

xmin = 0
xmax = 2 * np.pi
xticks = np.linspace(xmin, xmax, 5)
xticks_labels = ['$0$', '$\dfrac{\pi}{2}$', '$\pi$', '$\dfrac{3\pi}{2}$', '$2\pi$']


fig = plt.figure(0, figsize=(8, 6), frameon=False)
# Respuesta en magnitud
ax = plt.subplot2grid((6, 4), (0, 0), rowspan=2, colspan=2)
ax.set_prop_cycle(cy)
plt.plot(w, magHs, lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(-2, 2)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$\textrm{Primer orden}$', fontsize=fs)
plt.ylabel(r"$\textrm{dB}$", fontsize=fs)

leg = plt.legend([r'$z={}$'.format(aa[0]), r'$z={}$'.format(aa[1])], loc=2, frameon=False, fontsize=fs, framealpha=1)


# Respuesta en fase
ax = plt.subplot2grid((6, 4), (2, 0), rowspan=2, colspan=2)
ax.set_prop_cycle(cy)
plt.plot(w, argHs, lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.ylabel(r"$\textrm{Radianes}$", fontsize=fs)

# Retardo de grupo
ax = plt.subplot2grid((6, 4), (4, 0), rowspan=2, colspan=2)
ax.set_prop_cycle(cy)
plt.plot(w, grdHs, lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.xlabel(r"$\omega$", fontsize=fs)
plt.ylabel(r"$\textrm{Muestras}$", fontsize=fs)

##### Pasatodo de segundo orden #####
r = 0.9
theta = np.pi / 4

# vector del numerador y denominador de la función del sistema
b = [r ** 2, -2 * r * np.cos(theta), 1]
a = [1, -2 * r * np.cos(theta), r ** 2]
# respuesta en frecuencia
w, H2 = signal.freqz(b, a, nw, whole=True)
# magnitud de la respuesta en frecuencia
magH2 = 20 * np.log10(np.abs(H2))
# fase de la respuesta en frecuencia
argH2 = np.angle(H2)
# retardo de grupo de la respuesta en frecuencia
_, grdH2 = signal.group_delay((b, a), w)

# Respuesta en magnitud
ax = plt.subplot2grid((6, 4), (0, 2), rowspan=2, colspan=2)
ax.set_prop_cycle(cy)
plt.plot(w, magH2, lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(-2, 2)
ax.yaxis.set_ticklabels([])
plt.title(r'$\textrm{Segundo orden}$', fontsize=fs)

# Respuesta en fase
ax = plt.subplot2grid((6, 4), (2, 2), rowspan=2, colspan=2)
ax.set_prop_cycle(cy)
plt.plot(w, argH2, lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
ax.yaxis.set_ticklabels([])

# Retardo de grupo
ax = plt.subplot2grid((6, 4), (4, 2), rowspan=2, colspan=2)
ax.set_prop_cycle(cy)
plt.plot(w, grdH2, lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
ax.yaxis.set_ticklabels([])
plt.xlabel(r"$\omega$", fontsize=fs)

# save as pdf image
plt.savefig('transform_analysis_all_pass_first_second_order.pdf', bbox_inches='tight')

plt.show()
