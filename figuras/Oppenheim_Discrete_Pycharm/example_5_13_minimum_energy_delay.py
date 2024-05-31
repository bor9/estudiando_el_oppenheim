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
# a: ceros del sistema de fase mínima: zd1, zd1^*, zd2, zd2^*
za1 = 0.9 * np.exp(1j * 0.6 * np.pi)
za2 = 1.25 * np.exp(1j * 0.8 * np.pi)

### Ceros y ganancia de cada sistema

# ceros del sistema de fase mínima H_a(z)
zz = np.array([za1, np.conjugate(za1), 1 / np.conjugate(za2), 1 / za2])
# ganancia
kk = np.abs(za2) ** 2
# orden del sistema
s_ord = zz.shape[0]
# matriz con los ceros de los 4 sistemas
# se agregan los ceros del sistema con de fase máxima H_b(z)
zz = np.vstack([zz, np.array([1 / np.conjugate(za1), 1 / za1, za2, np.conjugate(za2)])])
kk = np.append(kk, np.abs(za1) ** 2)
# se agregan los ceros de H_c(z)
zz = np.vstack([zz, np.array([1 / np.conjugate(za1), 1 / za1, 1 / np.conjugate(za2), 1 / za2])])
kk = np.append(kk, np.abs(za1 * za2) ** 2)
# se agregan los ceros de H_d(z)
zz = np.vstack([zz, np.array([za1, np.conjugate(za1), za2, np.conjugate(za2)])])
kk = np.append(kk, 1)

# número de sistemas
n_sis = zz.shape[0]

# vector de polos
p = np.zeros(s_ord)

# cálculo de los coeficientes del filtro a partir de los polos, ceros y la ganancia
bb = np.zeros((n_sis, s_ord + 1))
for i in np.arange(n_sis):
    bb[i, :], _ = signal.zpk2tf(zz[i, :], p, kk[i])

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

fig = plt.figure(0, figsize=(8, 8), frameon=False)

# H_a(z)
ax = plt.subplot2grid((8, 8), (0, 0), rowspan=4, colspan=4)

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
i = 0
plt.plot(zz[i, :].real, zz[i, :].imag, **zeros_marker_style, zorder=10)
# polos
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgecolor='k', markeredgewidth=2)
plt.plot(0, 0, **polos_marker_style)
# multiplicidad del polo
plt.text(-0.07, 0.07, '{:d}'.format(4), fontsize=fontsize2, ha='right', va='baseline')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')
# etiqueta de la gráfica
plt.text(xmin_ax, ymax_ax, '$H_a(z)$', fontsize=fontsize, ha='left', va='top')
plt.axis('off')

# H_b(z)
ax = plt.subplot2grid((8, 8), (0, 4), rowspan=4, colspan=4)

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
i = 1
plt.plot(zz[i, :].real, zz[i, :].imag, **zeros_marker_style, zorder=10)
# polos
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgecolor='k', markeredgewidth=2)
plt.plot(0, 0, **polos_marker_style)
# multiplicidad del polo
plt.text(-0.07, 0.07, '{:d}'.format(4), fontsize=fontsize2, ha='right', va='baseline')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')
# etiqueta de la gráfica
plt.text(xmin_ax, ymax_ax, '$H_b(z)$', fontsize=fontsize, ha='left', va='top')
plt.axis('off')


# H_c(z)
ax = plt.subplot2grid((8, 8), (4, 0), rowspan=4, colspan=4)

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
i = 2
plt.plot(zz[i, :].real, zz[i, :].imag, **zeros_marker_style, zorder=10)
# polos
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgecolor='k', markeredgewidth=2)
plt.plot(0, 0, **polos_marker_style)
# multiplicidad del polo
plt.text(-0.07, 0.07, '{:d}'.format(4), fontsize=fontsize2, ha='right', va='baseline')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')
# etiqueta de la gráfica
plt.text(xmin_ax, ymax_ax, '$H_c(z)$', fontsize=fontsize, ha='left', va='top')
plt.axis('off')

# H_d(z)
ax = plt.subplot2grid((8, 8), (4, 4), rowspan=4, colspan=4)

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
i = 3
plt.plot(zz[i, :].real, zz[i, :].imag, **zeros_marker_style, zorder=10)
# polos
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgecolor='k', markeredgewidth=2)
plt.plot(0, 0, **polos_marker_style)
# multiplicidad del polo
plt.text(-0.07, 0.07, '{:d}'.format(4), fontsize=fontsize2, ha='right', va='baseline')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')
# etiqueta de la gráfica
plt.text(xmin_ax, ymax_ax, '$H_d(z)$', fontsize=fontsize, ha='left', va='top')
plt.axis('off')

# save as pdf image
plt.savefig('example_5_13_minimum_energy_zero_pole_plot.pdf', bbox_inches='tight')

#### RESPUESTAS AL IMPULSO

# cantidad del relleno de ceros al cominezo y al final
n_pad = 2
hh = np.concatenate((np.zeros((n_sis, n_pad)), bb, np.zeros((n_sis, n_pad))), axis=1)
n = np.arange(hh.shape[1]) - n_pad

# rango de los ejes
ymax = np.amax(hh)
ymin = np.amin(hh)

delta_n = 1
nmin_ax = n[0] - delta_n
nmax_ax = n[-1] + delta_n
delta_y = 0.5
ymax_ax = ymax + delta_y
ymin_ax = ymin - delta_y

baseline = -0.5
fontsize1 = 10
fontsize2 = 13
y_sep = 0.4
delta_l = 0.2


fig = plt.figure(1, figsize=(10, 5), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
i = 0
(markers, stemlines, bl) = plt.stem(n, hh[i, :], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)
plt.text(nmax_ax, baseline, '$n$', fontsize=fontsize1, ha='center', va='baseline')
for k in n:
    plt.text(k, baseline, '${}$'.format(k), fontsize=fontsize1, ha='center', va='baseline')
for j in np.arange(s_ord + 1):
    plt.text(j + delta_l, hh[i, j + n_pad], '${:,.2f}$'.format(hh[i, j + n_pad]), fontsize=fontsize1, ha='left', va='center')
plt.text(nmin_ax, ymax_ax, '$h_a[n]$', fontsize=fontsize2, ha='left', va='top')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (0, 2), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
i = 1
(markers, stemlines, bl) = plt.stem(n, hh[i, :], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)
plt.text(nmax_ax, baseline, '$n$', fontsize=fontsize1, ha='center', va='baseline')
for k in n:
    plt.text(k, baseline, '${}$'.format(k), fontsize=fontsize1, ha='center', va='baseline')
for j in np.arange(s_ord + 1):
    plt.text(j + delta_l, hh[i, j + n_pad], '${:,.2f}$'.format(hh[i, j + n_pad]), fontsize=fontsize1, ha='left', va='center')
plt.text(nmin_ax, ymax_ax, '$h_b[n]$', fontsize=fontsize2, ha='left', va='top')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 0), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
i = 2
(markers, stemlines, bl) = plt.stem(n, hh[i, :], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)
plt.text(nmax_ax, baseline, '$n$', fontsize=fontsize1, ha='center', va='baseline')
for k in n:
    plt.text(k, baseline, '${}$'.format(k), fontsize=fontsize1, ha='center', va='baseline')
for j in np.arange(s_ord + 1):
    plt.text(j + delta_l, hh[i, j + n_pad], '${:,.2f}$'.format(hh[i, j + n_pad]), fontsize=fontsize1, ha='left', va='center')
plt.text(nmin_ax, ymax_ax, '$h_c[n]$', fontsize=fontsize2, ha='left', va='top')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 2), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
i = 3
(markers, stemlines, bl) = plt.stem(n, hh[i, :], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)
plt.text(nmax_ax, baseline, '$n$', fontsize=fontsize1, ha='center', va='baseline')
for k in n:
    plt.text(k, baseline, '${}$'.format(k), fontsize=fontsize1, ha='center', va='baseline')
for j in np.arange(s_ord + 1):
    plt.text(j + delta_l, hh[i, j + n_pad], '${:,.2f}$'.format(hh[i, j + n_pad]), fontsize=fontsize1, ha='left', va='center')
plt.text(nmin_ax, ymax_ax, '$h_d[n]$', fontsize=fontsize2, ha='left', va='top')
plt.axis('off')

# save as pdf image
plt.savefig('example_5_13_minimum_energy_impulse_responses.pdf', bbox_inches='tight')

#### ENERGIA PARCIAL

E = np.cumsum(np.square(hh), axis=1)

nmin = 0
nmax = n[-2]
fs = 11

fig = plt.figure(2, figsize=(6, 4), frameon=False)
ax = fig.add_subplot(111)
plt.plot(n[n_pad:-1], E[:, n_pad:-1].T, ls='-', marker='s', markersize=4)
plt.xlim([nmin, nmax])
plt.ylim(bottom=0)
leg = plt.legend([r'$E_a[n]\textrm{ (fase m\'inima)}$', r'$E_b[n]\textrm{ (fase m\'axima)}$',
                  r'$E_c[n]$', r'$E_d[n]$'], loc=4, frameon=False, fontsize=fs, framealpha=1)
plt.xlabel(r'$n$', fontsize=fs)
plt.ylabel(r'$\textrm{Energ\'ia parcial}$', fontsize=fs)

plt.savefig('example_5_13_minimum_energy_partial_energy.pdf', bbox_inches='tight')

plt.show()
