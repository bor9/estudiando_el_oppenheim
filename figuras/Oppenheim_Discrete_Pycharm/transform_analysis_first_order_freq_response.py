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
thetas = [0, np.pi/2, np.pi]
r = 0.9
nw = 512

nt = len(thetas)
# vector del denominador de la función del sistema
a = 1
# inicialización de matrices
magHs = np.zeros((nw, nt))
argHs = np.zeros((nw, nt))
grdHs = np.zeros((nw, nt))
for i in np.arange(nt):
    # vector del numerador de la función del sistema
    b = [1, -r * np.exp(1j * thetas[i])]
    # respuesta en frecuencia
    w, H = signal.freqz(b, a, nw, whole=True)
    # magnitud de la respuesta en frecuencia
    magHs[:, i] = 20 * np.log10(np.abs(H))
    # fase de la respuesta en frecuencia
    argHs[:, i] = np.unwrap(np.angle(H))
    # retardo de grupo de la respuesta en frecuencia
    _, grdHs[:, i] = signal.group_delay((b, a), w)


########## Gráficas ##########

fs = 11  # fontsize
cy = cycler('color', ['black', 'red', 'blue', 'green'])
y_label_coords = -0.07

xmin = 0
xmax = 2 * np.pi
xticks = np.linspace(xmin, xmax, 5)
xticks_labels = ['$0$', '$\dfrac{\pi}{2}$', '$\pi$', '$\dfrac{3\pi}{2}$', '$2\pi$']


fig = plt.figure(0, figsize=(7, 7), frameon=False)
# Respuesta en magnitud
ax = plt.subplot2grid((6, 1), (0, 0), rowspan=2, colspan=1)
ax.set_prop_cycle(cy)
plt.plot(w, magHs, lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$\textrm{Magnitud, respuesta en fase y retardo de grupo}$', fontsize=fs)
plt.ylabel(r"$\textrm{dB}$", fontsize=fs)

# Respuesta en fase
ax = plt.subplot2grid((6, 1), (2, 0), rowspan=2, colspan=1)
ax.set_prop_cycle(cy)
plt.plot(w, argHs, lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.ylabel(r"$\textrm{Radianes}$", fontsize=fs)

# Retardo de grupo
ax = plt.subplot2grid((6, 1), (4, 0), rowspan=2, colspan=1)
ax.set_prop_cycle(cy)
plt.plot(w, grdHs, lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.xlabel(r"$\omega$", fontsize=fs)
plt.ylabel(r"$\textrm{Muestras}$", fontsize=fs)

leg = plt.legend([r'$\theta=0$', r'$\theta=\dfrac{\pi}{2}$', r'$\theta=\pi$'], loc=0, bbox_to_anchor=[0.9, 0.7],
                 frameon=False, fontsize=fs, framealpha=1)

# save as pdf image
plt.savefig('transform_analysis_first_order_freq_response_theta.pdf', bbox_inches='tight')

##
## Lo mismo pero variando r
##

rs = [0.5, 0.7, 0.9, 1]
theta = np.pi
nw = 512

nr = len(rs)
# vector del denominador de la función del sistema
a = 1
# inicialización de matrices
magHs = np.zeros((nw, nr))
argHs = np.zeros((nw, nr))
grdHs = np.zeros((nw, nr))
for i in np.arange(nr):
    # vector del numerador de la función del sistema
    b = [1, -rs[i] * np.exp(1j * theta)]
    # respuesta en frecuencia
    w, H = signal.freqz(b, a, nw, whole=True)
    # magnitud de la respuesta en frecuencia
    magHs[:, i] = 20 * np.log10(np.abs(H))
    # fase de la respuesta en frecuencia
    argHs[:, i] = np.unwrap(np.angle(H))
    # retardo de grupo de la respuesta en frecuencia
    _, grdHs[:, i] = signal.group_delay((b, a), w)


fig = plt.figure(1, figsize=(7, 7), frameon=False)
# Respuesta en magnitud
ax = plt.subplot2grid((6, 1), (0, 0), rowspan=2, colspan=1)
ax.set_prop_cycle(cy)
plt.plot(w, magHs, lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(-30, 10)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$\textrm{Magnitud, respuesta en fase y retardo de grupo}$', fontsize=fs)
plt.ylabel(r"$\textrm{dB}$", fontsize=fs)
leg = plt.legend([r'$r=0.5$', r'$r=0.7$', r'$r=0.9$', r'$r=1$'], loc=0, bbox_to_anchor=[0.9, 0.7],
                 frameon=False, fontsize=fs, framealpha=1)

# Respuesta en fase
ax = plt.subplot2grid((6, 1), (2, 0), rowspan=2, colspan=1)
ax.set_prop_cycle(cy)
plt.plot(w, argHs, lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.ylabel(r"$\textrm{Radianes}$", fontsize=fs)

# Retardo de grupo
ax = plt.subplot2grid((6, 1), (4, 0), rowspan=2, colspan=1)
ax.set_prop_cycle(cy)
plt.plot(w, grdHs[:, :-1], lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.xlabel(r"$\omega$", fontsize=fs)
plt.ylabel(r"$\textrm{Muestras}$", fontsize=fs)

# save as pdf image
plt.savefig('transform_analysis_first_order_freq_response_r.pdf', bbox_inches='tight')




plt.show()
