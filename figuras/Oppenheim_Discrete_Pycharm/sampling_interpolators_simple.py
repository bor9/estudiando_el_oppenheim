import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy import signal
import math

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# auxiliar function for plot ticks of equal length in x and y axis despite its scales.
def convert_display_to_data_coordinates(transData, length=10):
    # create a transform which will take from display to data coordinates
    inv = transData.inverted()
    # transform from display coordinates to data coordinates in x axis
    data_coords = inv.transform([(0, 0), (length, 0)])
    # get the length of the segment in data units
    yticks_len = data_coords[1, 0] - data_coords[0, 0]
    # transform from display coordinates to data coordinates in y axis
    data_coords = inv.transform([(0, 0), (0, length)])
    # get the length of the segment in data units
    xticks_len = data_coords[1, 1] - data_coords[0, 1]
    return xticks_len, yticks_len


def cubic_interpolation_filter(L, a):
    h_cub = np.zeros((2 * L,))
    nh = np.arange(L + 1)
    h_cub[nh] = (a + 2) * (np.absolute(nh / L) ** 3) - (a + 3) * (np.absolute(nh / L) ** 2) + 1
    nh = np.arange(L + 1, 2 * L)
    h_cub[nh] = a * (np.absolute(nh / L) ** 3) - 5 * a * (np.absolute(nh / L) ** 2) + 8 * a * np.absolute(
        nh / L) - 4 * a
    return np.concatenate([np.flip(h_cub[1:]), h_cub])


# se crea una señal suave definiendo algunos puntos xn en los instantes n
# y luego empleando interpolación cúbica.
xn = [3.5, 3, 2.6, 3.4, 3.2, 3, 4.4, 3.5, 2, 3, 4, 1.5, 2, 2.2, 3, 4.5, 3.8, 2.7, 2.5]
N = len(xn)
n = np.arange(N)
f = interpolate.interp1d(n, xn, kind='cubic')
Ni = 4000
ni = np.linspace(0, N-1, Ni)
# xni: señal, ni: indice de muestras 0: Ni - 1
xni = f(ni)
# se reduce la amplitud de las muestras
xni -= 1
ni = np.arange(Ni)

#
# parametros
#
# muestra correspondiente al tiempo 0 en la señal continua
s0 = 2200
# periodo de muestreo (muestras)
T = 200

# s son los indices de las muestras en nT en la señal continua xni
# n son los indices de la señal muestreada (entero)
n = np.arange(-math.floor(s0 / T), math.floor((Ni - s0) / T))
s = s0 + n * T

# xn es la señal discreta: xn = xc[nT]
xn = xni[s]

# xc es la señal continua, t es el tiempo de la señal continua
# alineado con n
xc = xni[s[0]: s[-1] + 1]
t = np.linspace(n[0], n[-1], s[-1] - s[0] + 1)

#
# fin parametros
#

#
# Grafica de prueba
#
grey = [0.7, 0.7, 0.7]
fig = plt.figure(0, figsize=(9, 3), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
# señal continua
plt.plot(t, xc, '-', lw=1, color=grey)
# señal discreta
(markers, stemlines, bl) = plt.stem(n, xn, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

#
# sobremuestreo
#

# factor de sobremuestreo
L = 5

xe = np.zeros((L * (xn.shape[0] + 1) - 1, ))
xe[np.arange(1, xn.shape[0] + 1) * L - 1] = xn

# filtros interpoladores
# lineal
h_lin = signal.windows.triang(2 * L - 1)
# cúbico
a = -0.5
h_cub = cubic_interpolation_filter(L, a) # largo = 4L - 1

# respuestas en frecuencia
N = 512
w, H_lin = signal.freqz(h_lin, a=1, worN=N, whole=False)
_, H_cub = signal.freqz(h_cub, a=1, worN=N, whole=False)

# filtrado
x_lin = signal.lfilter(h_lin, 1, xe, axis=-1, zi=None)
x_cub = signal.lfilter(h_cub, 1, xe, axis=-1, zi=None)

n = np.arange(xe.shape[0])

fig = plt.figure(1, figsize=(9, 3), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
# señal discreta
(markers, stemlines, bl) = plt.stem(n - (L - 1), x_lin, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
(markers, stemlines, bl) = plt.stem(n - (2 * L - 1), x_cub, linefmt='r', markerfmt='sr', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
(markers, stemlines, bl) = plt.stem(n, xe, linefmt='b', markerfmt='sb', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

#
# Gráficas
#
fontsize = 11
fontsize2 = 12
# x y ticks labels margin
xtm = -0.15

nmax = 2 * L
nmin = -nmax
nmax_ax = nmax + 2
nmin_ax = -nmax_ax

n = np.arange(nmin, nmax + 1)
ymin = -0.2
ymax = 1.2

h_lin_aux = np.zeros(n.shape)
h_lin_aux[nmax - (L-1): nmax + L] = h_lin
h_cub_aux = np.zeros(n.shape)
h_cub_aux[nmax - (2 * L-1): nmax + 2 * L] = h_cub

fig = plt.figure(2, figsize=(9, 5), frameon=False)
ax = plt.subplot2grid((2, 4), (0, 0), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, h_lin_aux, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-L, xtm, '$-L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
dx = 0.4
for i in np.arange(L):
    if i == 0:
        plt.text(i + dx, h_lin[L - 1 + i], '$1$', fontsize=fontsize, ha='left', va='bottom')
    else:
        plt.text(i + dx, h_lin[L - 1 + i], '${}/{}$'.format(L - i, L), fontsize=fontsize, ha='left', va='bottom')
plt.text(nmin_ax, 1, '$\\textrm{Lineal}$\n$L=5$', fontsize=fontsize2, ha='left', va='center')
plt.title(r'$\textrm{Respuesta al impulso}$', fontsize=fontsize2)
plt.axis('off')

ax = plt.subplot2grid((2, 4), (1, 0), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, h_cub_aux, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# etiquetas
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * L, xtm, '$2L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-2 * L, xtm, '$-2L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(dx, 1, '$1$', fontsize=fontsize, ha='left', va='bottom')
plt.text(nmin_ax, 1, '$\mathrm{C\\acute{u}bico}$\n$L=5,\,a=-0.5$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')


# respuesta en frecuencia
# x y ticks labels margin
xtm = -0.25
xtm2 = -0.39
ytm = -0.1
display_length = 6

xmin = 0
xmax = np.pi
dx = 0.4
xmax_ax = xmax + 0.3
xmin_ax = xmin - dx
ymax = L
ymin = 0
dy = 0.5
ymax_ax = ymax + dy
ymin_ax = ymin - dy

ax = plt.subplot2grid((2, 4), (0, 2), rowspan=2, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
# ideal
plt.plot([0, np.pi / L], [L, L], ls='-', color=grey, lw=1.5, label='$\mathrm{Ideal}$')
plt.plot([np.pi / L, np.pi / L], [0, L], ls='-', color=grey, lw=1.5)
plt.plot(w, np.absolute(H_lin), 'k-', lw=2, label='$\mathrm{Lineal}$')
plt.plot(w, np.absolute(H_cub), 'r-', lw=2, label='$\mathrm{C\\acute{u}bico}$')

# ticks y  etiquetas
plt.text(ytm, xtm, '$0$', fontsize=fontsize, ha='right', va='baseline')
for i in np.arange(1, L):
    plt.plot([i * np.pi / L, i * np.pi / L], [0, ytl], 'k-', lw=1)
    si = '' if i == 1 else '{}'.format(i)
    plt.text(i * np.pi / L, xtm2, '$\dfrac{{{}\pi}}{{L}}$'.format(si), fontsize=fontsize, ha='center', va='baseline')
plt.plot([np.pi, np.pi], [0, ytl], 'k-', lw=1)
plt.text(np.pi, xtm, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
#plt.text(ytm, L, '${}$'.format(L), fontsize=fontsize, ha='right', va='center')
plt.text(ytm, L, '$L$', fontsize=fontsize, ha='right', va='center')

plt.title(r'$\textrm{Respuesta en frecuencia}$', fontsize=fontsize2)
leg = plt.legend(loc=1, frameon=False, fontsize=fontsize2, framealpha=1)

plt.axis('off')
plt.savefig('sampling_interpolators_responses.pdf', bbox_inches='tight')

#
# Gráficas con interpolación
#
# muestra central
N1 = 14
n0 = N1 * L - 1
# cantidad de muestras a ambos lados
N2 = 5
nl = N2 * L
n = np.arange(0, 2 * nl + 1)
# señal extrapolada y señales filtradas
xe = xe[n0 - nl: n0 + nl + 1]
x_lin = x_lin[n0 - nl + L - 1: n0 + nl + L]
x_cub = x_cub[n0 - nl + 2 * L - 1: n0 + nl + 2 * L]
xe_aux = 1 * xe # esto se necesita para crear una copia profunda
xe_aux[xe_aux == 0] = 'nan'

# Parámetros de las gráficas
# x y ticks labels margin
xtm = -0.6

nmax = n[-1]
nmin = n[0]
dn = 2
nmax_ax = nmax + 3
nmin_ax = nmin - dn
ymin = -0.5
ymax = np.max(xe) + 0.2

ticks = np.arange(0, n[-1] + 1, L)
ticks_labels = ['${}L$'.format(t // L) for t in ticks]
ticks_labels[0] = 0
ticks_labels[1] = '$L$'

#### Señal continua
N1 = 14
N2 = 5
Ni = N1 - N2 - 1
Nf = N1 + N2 - 1
xc = xni[s[Ni]: s[Nf] + 1]
tt = np.linspace(n[0], n[-1], xc.shape[0])

fig = plt.figure(3, figsize=(9, 5), frameon=False)
ax = plt.subplot2grid((3, 4), (0, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, xe, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
plt.plot(tt, xc, ls='-', color=grey, lw=1, zorder=-10)
for i, t in enumerate(ticks):
    plt.text(t, xtm, ticks_labels[i], fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, ymax, '$x_e[n]$\n$L=5$', fontsize=fontsize2, ha='left', va='top')
plt.axis('off')

ax = plt.subplot2grid((3, 4), (1, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x_lin, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
plt.plot(tt, xc, ls='-', color=grey, lw=1, zorder=-10)
(markers, stemlines, bl) = plt.stem(n, xe_aux, linefmt='r', markerfmt='sr', use_line_collection=True)
plt.setp(markers, markersize=3.5, zorder=10)
plt.setp(bl, visible=False)
for i, t in enumerate(ticks):
    plt.text(t, xtm, ticks_labels[i], fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, ymax, r'$x_\textrm{lin}[n]$', fontsize=fontsize2, ha='left', va='top')
plt.axis('off')

ax = plt.subplot2grid((3, 4), (2, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x_cub, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
plt.plot(tt, xc, ls='-', color=grey, lw=1, zorder=-10)
(markers, stemlines, bl) = plt.stem(n, xe_aux, linefmt='r', markerfmt='sr', use_line_collection=True)
plt.setp(markers, markersize=3.5, zorder=10)
plt.setp(bl, visible=False)
for i, t in enumerate(ticks):
    plt.text(t, xtm, ticks_labels[i], fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, ymax, r'$x_\textrm{cub}[n]$', fontsize=fontsize2, ha='left', va='top')
plt.axis('off')

plt.savefig('sampling_interpolators_interpolation.pdf', bbox_inches='tight')


plt.show()
