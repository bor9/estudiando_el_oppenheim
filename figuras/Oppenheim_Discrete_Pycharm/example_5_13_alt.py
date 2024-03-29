import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Implementación del ejemplo 5.13 del libro
# A. V. Oppenheim and R. W. Schafer, Discrete-Time Signal Processing. Prentice Hall, 3rd ed., 2009.

# parámetros
# ceros del sistema H_d(z): zd1, zd1^*, zd2, zd2^*
zd1 = 0.9 * np.exp(1j * 0.6 * np.pi)
zd2 = 1.25 * np.exp(1j * 0.8 * np.pi)

### Sistema Hd(z)

# Ceros
z = np.array([zd1, np.conjugate(zd1), zd2, np.conjugate(zd2)])
# Polos
p = np.zeros((4))
# Ganancia
k = 1

# Cálculo de los coeficientes del filtro a partir de los polos, ceros y la ganancia
b_d, a_d = signal.zpk2tf(z, p, k)

### Sistema de fase mínima H_min(z)

# Ceros
z = np.array([zd1, np.conjugate(zd1), 1 / np.conjugate(zd2), 1 / zd2])
# Ganancia
k = np.abs(zd2) ** 2
b_min, a_min = signal.zpk2tf(z, p, k)

### Sistema pasatodo H_ap(z)
z = np.array([zd2, np.conjugate(zd2)])
# Polos
p = np.array([1 / np.conjugate(zd2), 1 / zd2])
# Ganancia
k = 1 / (np.abs(zd2) ** 2)
b_ap, a_ap = signal.zpk2tf(z, p, k)

# número de frecuencias
nw = 512

# RESPUESTA EN FRECUENCIA DE H_d(z)
w, H_d = signal.freqz(b_d, a_d, nw, whole=True)
# magnitud de la respuesta en frecuencia
magH_d = 20 * np.log10(np.abs(H_d))
# fase de la respuesta en frecuencia
argH_d = np.unwrap(np.angle(H_d))
#argHd = np.angle(Hd)
# retardo de grupo de la respuesta en frecuencia
_, grdH_d = signal.group_delay((b_d, a_d), w)

# RESPUESTA EN FRECUENCIA DE H_min(z)
_, H_min = signal.freqz(b_min, a_min, nw, whole=True)
# magnitud de la respuesta en frecuencia
magH_min = 20 * np.log10(np.abs(H_min))
# fase de la respuesta en frecuencia
argH_min = np.unwrap(np.angle(H_min))
#argHd = np.angle(Hd)
# retardo de grupo de la respuesta en frecuencia
_, grdH_min = signal.group_delay((b_min, a_min), w)

# RESPUESTA EN FRECUENCIA DE H_ap(z)
_, H_ap = signal.freqz(b_ap, a_ap, nw, whole=True)
# magnitud de la respuesta en frecuencia
magH_ap = 20 * np.log10(np.abs(H_ap))
# fase de la respuesta en frecuencia
argH_ap = np.unwrap(np.angle(H_ap))
#argHd = np.angle(Hd)
# retardo de grupo de la respuesta en frecuencia
_, grdH_ap = signal.group_delay((b_ap, a_ap), w)

########## Gráficas ##########

fs = 11  # fontsize
y_label_coords = -0.08

xmin = 0
xmax = 2 * np.pi
xticks = np.linspace(xmin, xmax, 5)
xticks_labels = ['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$']
xticks_labels = ['$0$', '$\dfrac{\pi}{2}$', '$\pi$', '$\dfrac{3\pi}{2}$', '$2\pi$']

# Se establecen los valores máximos y mínimos del eje y para todas las gráficas
dy = 1
ymax_mag = np.amax(np.concatenate((magH_d, magH_min, magH_ap))) + dy
ymin_mag = np.amin(np.concatenate((magH_d, magH_min, magH_ap))) - dy
ymax_ph = np.amax(np.concatenate((argH_d, argH_min, argH_ap))) + dy
ymin_ph = np.amin(np.concatenate((argH_d, argH_min, argH_ap))) - dy
ymax_gr = np.amax(np.concatenate((grdH_d, grdH_min, grdH_ap))) + dy
ymin_gr = np.amin(np.concatenate((grdH_d, grdH_min, grdH_ap))) - dy

##### H_d(z) #####
fig = plt.figure(0, figsize=(7, 7), frameon=False)
# Respuesta en magnitud
ax = plt.subplot2grid((6, 1), (0, 0), rowspan=2, colspan=1)
plt.plot(w, magH_d, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_mag, ymax_mag)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$H_d(z):\,\textrm{magnitud, respuesta en fase desenrollada y retardo de grupo}$', fontsize=fs)
plt.ylabel(r"$\textrm{dB}$", fontsize=fs)

# Respuesta en fase
ax = plt.subplot2grid((6, 1), (2, 0), rowspan=2, colspan=1)
plt.plot(w, argH_d, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_ph, ymax_ph)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.ylabel(r"$\textrm{Radianes}$", fontsize=fs)

# Retardo de grupo
ax = plt.subplot2grid((6, 1), (4, 0), rowspan=2, colspan=1)
plt.plot(w, grdH_d, 'k-', lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_gr, ymax_gr)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.xlabel(r"$\omega$", fontsize=fs)
plt.ylabel(r"$\textrm{Muestras}$", fontsize=fs)

# save as pdf image
plt.savefig('example_5_13_H_d.pdf', bbox_inches='tight')

##### H_min(z) #####
fig = plt.figure(1, figsize=(7, 7), frameon=False)
# Respuesta en magnitud
ax = plt.subplot2grid((6, 1), (0, 0), rowspan=2, colspan=1)
plt.plot(w, magH_min, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_mag, ymax_mag)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'${H_d}_\textrm{min}(z):\,\textrm{magnitud, respuesta en fase desenrollada y retardo de grupo}$', fontsize=fs)
plt.ylabel(r"$\textrm{dB}$", fontsize=fs)

# Respuesta en fase
ax = plt.subplot2grid((6, 1), (2, 0), rowspan=2, colspan=1)
plt.plot(w, argH_min, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_ph, ymax_ph)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.ylabel(r"$\textrm{Radianes}$", fontsize=fs)

# Retardo de grupo
ax = plt.subplot2grid((6, 1), (4, 0), rowspan=2, colspan=1)
plt.plot(w, grdH_min, 'k-', lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_gr, ymax_gr)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.xlabel(r"$\omega$", fontsize=fs)
plt.ylabel(r"$\textrm{Muestras}$", fontsize=fs)

# save as pdf image
plt.savefig('example_5_13_H_min.pdf', bbox_inches='tight')

##### H_ap(z) #####
fig = plt.figure(2, figsize=(7, 7), frameon=False)
# Respuesta en magnitud
ax = plt.subplot2grid((6, 1), (0, 0), rowspan=2, colspan=1)
plt.plot(w, magH_ap, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_mag, ymax_mag)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$H_\textrm{ap}(z):\,\textrm{magnitud, respuesta en fase desenrrollada y retardo de grupo}$', fontsize=fs)
plt.ylabel(r"$\textrm{dB}$", fontsize=fs)

# Respuesta en fase
ax = plt.subplot2grid((6, 1), (2, 0), rowspan=2, colspan=1)
plt.plot(w, argH_ap, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_ph, ymax_ph)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.ylabel(r"$\textrm{Radianes}$", fontsize=fs)

# Retardo de grupo
ax = plt.subplot2grid((6, 1), (4, 0), rowspan=2, colspan=1)
plt.plot(w, grdH_ap, 'k-', lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_gr, ymax_gr)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.xlabel(r"$\omega$", fontsize=fs)
plt.ylabel(r"$\textrm{Muestras}$", fontsize=fs)

# save as pdf image
plt.savefig('example_5_13_H_ap.pdf', bbox_inches='tight')

plt.show()
