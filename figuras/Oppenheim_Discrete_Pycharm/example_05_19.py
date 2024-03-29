import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Implementación del ejemplo 5.19 del libro
# A. V. Oppenheim and R. W. Schafer, Discrete-Time Signal Processing. Prentice Hall, 3rd ed., 2009.

# parámetros
# ceros del sistema H_d(z): zd1, zd1^*, zd2, zd2^*
zd1 = 0.9 * np.exp(1j * 0.6 * np.pi)
zd2 = 1.25 * np.exp(1j * 0.8 * np.pi)

### Sistema de fase mínima H_min(z)

# Ceros
z_min = np.array([zd1, np.conjugate(zd1), 1 / np.conjugate(zd2), 1 / zd2])
# Polos
p_min = np.zeros(4)
# Ganancia
k_min = np.abs(zd2) ** 2
b_min, a_min = signal.zpk2tf(z_min, p_min, k_min)

### Sistema de fase máxima H_max(z)

# Ceros
z_max = np.array([1 / zd1, 1 / np.conjugate(zd1), zd2, np.conjugate(zd2)])
# Polos
p_max = np.zeros(4)
# Ganancia
k_max = np.abs(zd1) ** 2
b_max, a_max = signal.zpk2tf(z_max, p_max, k_max)

### Sistema global: H_min(z)H_max(z)

# Ceros
z = np.concatenate((z_min, z_max))
# Polos
p = np.concatenate((p_min, p_max))
# Ganancia
k = k_min * k_max
b, a = signal.zpk2tf(z, p, k)

# número de frecuencias
nw = 512

# RESPUESTA EN FRECUENCIA DE H_min(z)
w, H_min = signal.freqz(b_min, a_min, nw, whole=True)
# magnitud de la respuesta en frecuencia
magH_min = 20 * np.log10(np.abs(H_min))
# fase de la respuesta en frecuencia
argH_min = np.unwrap(np.angle(H_min))
#argHd = np.angle(Hd)
# retardo de grupo de la respuesta en frecuencia
_, grdH_min = signal.group_delay((b_min, a_min), w)

# RESPUESTA EN FRECUENCIA DE H_max(z)
_, H_max = signal.freqz(b_max, a_max, nw, whole=True)
# magnitud de la respuesta en frecuencia
magH_max = 20 * np.log10(np.abs(H_max))
# fase de la respuesta en frecuencia
argH_max = np.unwrap(np.angle(H_max))
#argHd = np.angle(Hd)
# retardo de grupo de la respuesta en frecuencia
_, grdH_max = signal.group_delay((b_max, a_max), w)

# RESPUESTA EN FRECUENCIA DE H(z)= H_min(z)H_max(z)
_, H = signal.freqz(b, a, nw, whole=True)
# magnitud de la respuesta en frecuencia
magH = 20 * np.log10(np.abs(H))
# fase de la respuesta en frecuencia
argH = np.unwrap(np.angle(H))
#argHd = np.angle(Hd)
# retardo de grupo de la respuesta en frecuencia
_, grdH = signal.group_delay((b, a), w)

########## Gráficas ##########

fs = 11  # fontsize
y_label_coords = -0.21

xmin = 0
xmax = 2 * np.pi
xticks = np.linspace(xmin, xmax, 5)
xticks_labels = ['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$']
xticks_labels = ['$0$', '$\dfrac{\pi}{2}$', '$\pi$', '$\dfrac{3\pi}{2}$', '$2\pi$']

# Se establecen los valores máximos y mínimos del eje y para todas las gráficas
dy = 1
ymax_mag = np.amax(np.concatenate((magH_min, magH_max, magH))) + dy
ymin_mag = np.amin(np.concatenate((magH_min, magH_max, magH))) - dy
ymax_ph = np.amax(np.concatenate((argH_min, argH_max, argH))) + dy
ymin_ph = np.amin(np.concatenate((argH_min, argH_max, argH))) - dy
ymax_gr = np.amax(np.concatenate((grdH_min, grdH_max, grdH))) + dy
ymin_gr = np.amin(np.concatenate((grdH_min, grdH_max, grdH))) - dy

##### H_min(z) #####
fig = plt.figure(0, figsize=(7, 7), frameon=False)

# Respuesta en magnitud
ax = plt.subplot2grid((6, 6), (0, 0), rowspan=2, colspan=2)
plt.plot(w, magH_min, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_mag, ymax_mag)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$H_\textrm{min}(e^{j\omega})$', fontsize=fs)
plt.ylabel(r"$\textrm{dB}$", fontsize=fs)

# Respuesta en fase
ax = plt.subplot2grid((6, 6), (2, 0), rowspan=2, colspan=2)
plt.plot(w, argH_min, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_ph, ymax_ph)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.ylabel(r"$\textrm{Radianes}$", fontsize=fs)

# Retardo de grupo
ax = plt.subplot2grid((6, 6), (4, 0), rowspan=2, colspan=2)
plt.plot(w, grdH_min, 'k-', lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_gr, ymax_gr)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.xlabel(r"$\omega$", fontsize=fs)
plt.ylabel(r"$\textrm{Muestras}$", fontsize=fs)


##### H_max(z) #####
# Respuesta en magnitud
ax = plt.subplot2grid((6, 6), (0, 2), rowspan=2, colspan=2)
plt.plot(w, magH_max, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_mag, ymax_mag)
ax.yaxis.set_ticklabels([])
plt.title(r'$H_\textrm{max}(e^{j\omega})$', fontsize=fs)

# Respuesta en fase
ax = plt.subplot2grid((6, 6), (2, 2), rowspan=2, colspan=2)
plt.plot(w, argH_max, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_ph, ymax_ph)
ax.yaxis.set_ticklabels([])

# Retardo de grupo
ax = plt.subplot2grid((6, 6), (4, 2), rowspan=2, colspan=2)
plt.plot(w, grdH_max, 'k-', lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_gr, ymax_gr)
ax.yaxis.set_ticklabels([])
plt.xlabel(r"$\omega$", fontsize=fs)

##### H(z) #####
# Respuesta en magnitud
ax = plt.subplot2grid((6, 6), (0, 4), rowspan=2, colspan=2)
plt.plot(w, magH, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_mag, ymax_mag)
ax.yaxis.set_ticklabels([])
plt.title(r'$H(e^{j\omega})$', fontsize=fs)

# Respuesta en fase
ax = plt.subplot2grid((6, 6), (2, 4), rowspan=2, colspan=2)
plt.plot(w, argH, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_ph, ymax_ph)
ax.yaxis.set_ticklabels([])

# Retardo de grupo
ax = plt.subplot2grid((6, 6), (4, 4), rowspan=2, colspan=2)
plt.plot(w, grdH, 'k-', lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin_gr, ymax_gr)
ax.yaxis.set_ticklabels([])
plt.xlabel(r"$\omega$", fontsize=fs)

# save as pdf image
plt.savefig('example_05_19_freq_responses.pdf', bbox_inches='tight')


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

# estilo de las marcas de los ceros y los polos
zeros_marker_style = dict(marker='o', linestyle='', markersize=7, markerfacecolor='w',
                          markeredgewidth=1.5)
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgecolor='k', markeredgewidth=2)

fig = plt.figure(1, figsize=(10, 4), frameon=False)

### H_min(z)
ax = plt.subplot2grid((4, 12), (0, 0), rowspan=4, colspan=4)

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
plt.plot(z_min.real, z_min.imag, **zeros_marker_style, markeredgecolor='b', zorder=10)
# polos
plt.plot(0, 0, **polos_marker_style)
for zz in z_min:
    plt.plot([0, zz.real], [0, zz.imag], 'k--', zorder=1, lw=1)

# multiplicidad del polo
plt.text(0.07, 0.07, '{:d}'.format(4), fontsize=fontsize2, ha='left', va='baseline')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')

# titulo
plt.text(xmax_ax, ymax_ax, r'$H_\textrm{min}(z)$', fontsize=fontsize, ha='right', va='baseline')

plt.axis('off')

### H_max(z)
ax = plt.subplot2grid((4, 12), (0, 4), rowspan=4, colspan=4)

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
plt.plot(z_max.real, z_max.imag, **zeros_marker_style, markeredgecolor='r', zorder=10)
# polos
plt.plot(0, 0, **polos_marker_style)
for zz in z_max:
    plt.plot([0, zz.real], [0, zz.imag], 'k--', zorder=1, lw=1)

# multiplicidad del polo
plt.text(0.07, 0.07, '{:d}'.format(4), fontsize=fontsize2, ha='left', va='baseline')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')

# titulo
plt.text(xmax_ax, ymax_ax, r'$H_\textrm{max}(z)$', fontsize=fontsize, ha='right', va='baseline')

plt.axis('off')

### H_(z)
ax = plt.subplot2grid((4, 12), (0, 8), rowspan=4, colspan=4)

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
plt.plot(z_min.real, z_min.imag, **zeros_marker_style, markeredgecolor='b', zorder=10)
plt.plot(z_max.real, z_max.imag, **zeros_marker_style, markeredgecolor='r', zorder=10)
# polos
plt.plot(0, 0, **polos_marker_style)
for zz in z_max:
    plt.plot([0, zz.real], [0, zz.imag], 'k--', zorder=1, lw=1)

# multiplicidad del polo
plt.text(0.07, 0.07, '{:d}'.format(8), fontsize=fontsize2, ha='left', va='baseline')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')

# titulo
plt.text(xmax_ax, ymax_ax, r'$H(z)$', fontsize=fontsize, ha='right', va='baseline')

plt.axis('off')

# save as pdf image
plt.savefig('example_05_19_zero_pole_plot.pdf', bbox_inches='tight')

plt.show()
