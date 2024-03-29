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
z_d = np.array([zd1, np.conjugate(zd1), zd2, np.conjugate(zd2)])
# Polos
p = np.zeros(4)
# Ganancia
k = 1

# Cálculo de los coeficientes del filtro a partir de los polos, ceros y la ganancia
b_d, a_d = signal.zpk2tf(z_d, p, k)

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
y_label_coords = -0.17

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
ax = plt.subplot2grid((6, 6), (0, 0), rowspan=2, colspan=2)
plt.plot(w, magH_d, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_mag, ymax_mag)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$H_d(e^{j\omega})$', fontsize=fs)
plt.ylabel(r"$\textrm{dB}$", fontsize=fs)

# Respuesta en fase
ax = plt.subplot2grid((6, 6), (2, 0), rowspan=2, colspan=2)
plt.plot(w, argH_d, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_ph, ymax_ph)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.ylabel(r"$\textrm{Radianes}$", fontsize=fs)

# Retardo de grupo
ax = plt.subplot2grid((6, 6), (4, 0), rowspan=2, colspan=2)
plt.plot(w, grdH_d, 'k-', lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_gr, ymax_gr)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.xlabel(r"$\omega$", fontsize=fs)
plt.ylabel(r"$\textrm{Muestras}$", fontsize=fs)


##### H_min(z) #####
# Respuesta en magnitud
ax = plt.subplot2grid((6, 6), (0, 2), rowspan=2, colspan=2)
plt.plot(w, magH_min, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_mag, ymax_mag)
ax.yaxis.set_ticklabels([])
plt.title(r'${H_d}_\textrm{min}(e^{j\omega})$', fontsize=fs)

# Respuesta en fase
ax = plt.subplot2grid((6, 6), (2, 2), rowspan=2, colspan=2)
plt.plot(w, argH_min, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_ph, ymax_ph)
ax.yaxis.set_ticklabels([])

# Retardo de grupo
ax = plt.subplot2grid((6, 6), (4, 2), rowspan=2, colspan=2)
plt.plot(w, grdH_min, 'k-', lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_gr, ymax_gr)
ax.yaxis.set_ticklabels([])
plt.xlabel(r"$\omega$", fontsize=fs)

##### H_ap(z) #####
# Respuesta en magnitud
ax = plt.subplot2grid((6, 6), (0, 4), rowspan=2, colspan=2)
plt.plot(w, magH_ap, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_mag, ymax_mag)
ax.yaxis.set_ticklabels([])
plt.title(r'$H_\textrm{ap}(e^{j\omega})$', fontsize=fs)

# Respuesta en fase
ax = plt.subplot2grid((6, 6), (2, 4), rowspan=2, colspan=2)
plt.plot(w, argH_ap, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_ph, ymax_ph)
ax.yaxis.set_ticklabels([])

# Retardo de grupo
ax = plt.subplot2grid((6, 6), (4, 4), rowspan=2, colspan=2)
plt.plot(w, grdH_ap, 'k-', lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_gr, ymax_gr)
ax.yaxis.set_ticklabels([])
plt.xlabel(r"$\omega$", fontsize=fs)

# save as pdf image
plt.savefig('example_5_13_H_freq_responses.pdf', bbox_inches='tight')


#### DIGRAMA DE POLOS Y CEROS

# circulo |z|=1
M = 200
thetas = np.linspace(0, 2 * np.pi, M)
# circulo unidad
xd = np.cos(thetas)
yd = np.sin(thetas)

max_ax = 1.5
# axis parameters
xmin_ax = -max_ax
xmax_ax = max_ax
ymin_ax = -max_ax
ymax_ax = max_ax

# x ticks labels margin
xtm = -0.18
ytm = -0.1
# font size
fontsize = 14
fontsize2 = 10

fig = plt.figure(1, figsize=(4, 4), frameon=False)
ax = fig.add_subplot(111)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.gca().set_aspect('equal', adjustable='box')


# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

# círculo unidad
plt.plot(xd, yd, 'k', lw=2)
# ceros
zeros_marker_style = dict(marker='o', linestyle='', markersize=7, markerfacecolor='w', markeredgecolor='k',
                          markeredgewidth=1.5)
plt.plot(z_d.real, z_d.imag, **zeros_marker_style, zorder=10)
# polos
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgecolor='k', markeredgewidth=2)
plt.plot(0, 0, **polos_marker_style)
# multiplicidad del polo
plt.text(-0.07, 0.07, '{:d}'.format(4), fontsize=fontsize2, ha='right', va='baseline')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')

plt.axis('off')

# save as pdf image
plt.savefig('example_5_13_zero_pole_plot.pdf', bbox_inches='tight')

plt.show()
