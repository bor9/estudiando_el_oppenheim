import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Construcción del filtro del ejemplo de la sección 5.1.2 del libro
#  A. V. Oppenheim and R. W. Schafer, Discrete-Time Signal Processing. Prentice Hall, 3rd ed., 2009.

# Polos
c_k = 0.95 * np.exp(1j * (0.15 * np.pi + 0.02 * np.pi * np.arange(1, 5)))

# Construcción del vector de polos y ceros
# Ceros dobles: 1/c_k, 1/c^*_k. Polos dobles: c_k, c^*_k.
# Polos
p = np.concatenate((c_k, c_k, np.conjugate(c_k), np.conjugate(c_k)))
# Ceros
z = np.concatenate((1 / c_k, 1 / c_k, 1 / np.conjugate(c_k), 1 / np.conjugate(c_k)))
# Ganancia
k = np.abs(np.prod(c_k)) ** 4 # el 4 es porque los polos son dobles.

# Cálculo de los coeficientes del filtro a partir de los polos, ceros y la ganancia
b, a = signal.zpk2tf(z, p, k)

# Construcción de la señal a filtrar
M = 71
nx = np.arange(M)
wn = np.hamming(M)
w1 = 0.2 * np.pi
w2 = 0.4 * np.pi
x1 = wn * np.cos(w1 * nx)
x2 = wn * np.cos(w2 * nx - np.pi / 2)
N = 300
n = np.arange(N)
x = np.zeros(N)
x[0: M] = x1
x[M: 2 * M] = x2

# Filtrado
y = signal.lfilter(b, a, x)

# Respuesta en frecuencia del filtro
nw = 1024
w, H = signal.freqz(b, a, nw)
# Magnitud de la respuesta en frecuencia
magH = np.abs(H)
# fase de la respuesta en frecuencia
argH = np.angle(H)
# fase desenvuelta
argH_unw = np.unwrap(np.angle(H))
# retardo de grupo de la respuesta en frecuencia
_, grdH = signal.group_delay((b, a), w)
# cálculo manual del retardo de grupo - la funcion signal.group_delay no funciona en este caso (??)
grdH_alt = -np.diff(argH_unw) * nw / np.pi

# Magnitud del espectro de la señal de entrada
_, magX = np.abs(signal.freqz(x, 1, nw))



########## Gráficas ##########

fs = 11  # fontsize

# Magnitud y fase de la respuesta en frecuencia
xmin = 0
xmax = np.pi
xticks = np.linspace(0, 1, num=6, endpoint=True)
xticks_labels = ["${:.1f}\pi$".format(xt) if xt != 0 else "${:.0f}$".format(xt) for xt in xticks]
xticks_labels[-1] = "$\pi$"
xticks = xticks * np.pi
yticks = [0, 1]

fig = plt.figure(0, figsize=(8, 5), frameon=False)
# Respuesta en magnitud
ax = plt.subplot2grid((4, 1), (0, 0), rowspan=2, colspan=1)
plt.plot(w, magH, lw=1.5, color='k')
plt.xticks(xticks, [], usetex=True)
plt.yticks(yticks, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(0, 1.2)
ax.yaxis.set_label_coords(-0.05, 0.5)
plt.title(r'$\textrm{Magnitud y fase de la respuesta en frecuencia del sistema}$', fontsize=fs)
plt.ylabel(r"$|H(e^{j\omega})|$", fontsize=fs)

# Respuesta en fase
yticks = [-np.pi, 0, np.pi]
yticks_labels = ['$-\pi$', '$0$', '$-\pi$']
ax = plt.subplot2grid((4, 1), (2, 0), rowspan=2, colspan=1)
plt.plot(w, argH, lw=1.5, color='k')
plt.xticks(xticks, xticks_labels, usetex=True)
plt.yticks(yticks, yticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(-np.pi, np.pi)
ax.yaxis.set_label_coords(-0.05, 0.5)
plt.xlabel(r"$\omega\;\rm(rad)$", fontsize=fs)
plt.ylabel(r"$\textrm{ARG}[H(e^{j\omega})]$", fontsize=fs)
plt.savefig('transform_analysis_example_5_1_2_group_delay_frequency_response.pdf', bbox_inches='tight')

# Fase desenvuelta y retardo de grupo del sistema y repuesta en frecuencia de la entrada
ylabels_ycoord = -0.08

yticks = np.linspace(-15, 0, num=4, endpoint=True)
yticks_labels = ["${:.0f}\pi$".format(xt) if xt != 0 else "${:.0f}$".format(xt) for xt in yticks]
yticks = yticks * np.pi
fig = plt.figure(1, figsize=(8, 6), frameon=False)
# Respuesta en fase desenvuelta
ax = plt.subplot2grid((6, 1), (0, 0), rowspan=2, colspan=1)
plt.plot(w, argH_unw, lw=1.5, color='k')
plt.xticks(xticks, [], usetex=True)
plt.yticks(yticks, yticks_labels, usetex=True)
plt.xlim(xmin, xmax)
ax.yaxis.set_label_coords(ylabels_ycoord, 0.5)
plt.title(r'$\textrm{Fase desenrollada y retardo de grupo del sistema y magnitud del espectro de la entrada}$',
          fontsize=fs)
plt.ylabel(r"$\textrm{arg}[H(e^{j\omega})]$", fontsize=fs)
plt.grid()

# Retardo de grupo
ax = plt.subplot2grid((6, 1), (2, 0), rowspan=2, colspan=1)
plt.plot(w[:-1] + np.pi / nw, grdH_alt, lw=1.5, color='k')
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
ax.yaxis.set_label_coords(ylabels_ycoord, 0.5)
plt.ylabel(r"$\textrm{grd}[H(e^{j\omega})]$", fontsize=fs)
plt.grid()

# Magnitud del espectro de la entrada
ax = plt.subplot2grid((6, 1), (4, 0), rowspan=2, colspan=1)
plt.plot(w, magX, lw=1.5, color='k')
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.xticks(xticks, xticks_labels, usetex=True)
ax.yaxis.set_label_coords(ylabels_ycoord, 0.5)
plt.ylabel(r"$|X(e^{j\omega})|$", fontsize=fs)
plt.xlabel(r"$\omega\;\rm(rad)$", fontsize=fs)
plt.grid()
plt.savefig('transform_analysis_example_5_1_2_group_delay.pdf', bbox_inches='tight')

# Entrada y la salida en el tiempo

# retardo de grupo en las frecuencias de la entrada
ms = 2.5
xmin = 0
xmax = N
ymax = 1.1
fig = plt.figure(2, figsize=(8, 5), frameon=False)
# Entrada
ax = plt.subplot2grid((4, 1), (0, 0), rowspan=2, colspan=1)
plt.plot(n, x, ls='-', marker='s', color='k', markersize=ms, lw=1)
plt.xlim(xmin, xmax)
plt.ylim(-ymax, ymax)
ax.set_xticklabels([])
ax.yaxis.set_label_coords(-0.05, 0.5)
plt.title(r'$\textrm{Entrada }x[n]\textrm{ y salida }y[n]$', fontsize=fs)
plt.ylabel(r"$x[n]$", fontsize=fs)

ax = plt.subplot2grid((4, 1), (2, 0), rowspan=2, colspan=1)
plt.plot(n, y, ls='-', marker='s', color='k', markersize=ms, lw=1)
plt.xlim(xmin, xmax)
plt.ylim(-ymax, ymax)
ax.yaxis.set_label_coords(-0.05, 0.5)
plt.ylabel(r"$y[n]$", fontsize=fs)
plt.xlabel(r"$n$", fontsize=fs)

plt.savefig('transform_analysis_example_5_1_2_group_delay_input_output.pdf', bbox_inches='tight')

# Diagrama de ceros y polos

# circulo |z|=1
thetas = np.linspace(0, 2 * np.pi, 200)
# circulo unidad
xd = np.cos(thetas)
yd = np.sin(thetas)
# polos
pp = np.concatenate((c_k, np.conjugate(c_k)))
# Ceros
zz = np.concatenate((1 / c_k, 1 / np.conjugate(c_k)))

# valores maximos y minimos de los ejes
max_ax = 1.2
xmin = -max_ax
xmax = max_ax
ymin = -max_ax
ymax = max_ax

# axis parameters
xmin_ax = xmin
xmax_ax = xmax
ymin_ax = ymin
ymax_ax = ymax

# length of the ticks for all subplot (6 pixels)
display_length = 6  # in pixels
# x ticks labels margin
xtm = -0.13
ytm = -0.07
# font size
fontsize = 14
fontsize2 = 10

fig = plt.figure(3, figsize=(5, 5), frameon=False)
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
ms = 6
zeros_marker_style = dict(marker='o', linestyle='', markersize=7, markerfacecolor='w', markeredgecolor='k',
                          markeredgewidth=1.5)
plt.plot(zz.real, zz.imag, **zeros_marker_style, zorder=10)
# polos
polos_marker_style = dict(marker='x', linestyle='', markersize=6, markeredgecolor='k', markeredgewidth=1.5)
plt.plot(pp.real, pp.imag, **polos_marker_style)
# multiplicidades
dm = 0.1
for c in c_k:
    # posición del texto en el plano complejo
    tp = (np.abs(c) - dm) * np.exp(1j * np.angle(c))
    plt.text(tp.real, tp.imag, '2', fontsize=fontsize2, ha='center', va='center')
    plt.text(tp.real, -tp.imag, '2', fontsize=fontsize2, ha='center', va='center')
    tp = (1 / np.abs(c) + 0.09) * np.exp(1j * np.angle(c))
    plt.text(tp.real, tp.imag, '2', fontsize=fontsize2, ha='center', va='center')
    plt.text(tp.real, -tp.imag, '2', fontsize=fontsize2, ha='center', va='center')

# puntos 1
ms = 8
plt.plot(1, 0, 'k.', ms=ms)

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')

# circle labels
plt.text(1.03, xtm, '$1$', fontsize=fontsize, ha='left', va='baseline')

plt.axis('off')

# save as pdf image
plt.savefig('transform_analysis_example_5_1_2_group_delay_poles_zeros_plot.pdf', bbox_inches='tight')

plt.show()
