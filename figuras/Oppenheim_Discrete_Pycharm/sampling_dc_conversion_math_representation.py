import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import math
from matplotlib.lines import Line2D

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


# se crea una señal suave definiendo algunos puntos xn en los instantes n
# y luego empleando interpolación cúbica.
xn = [3.5, 3, 2.6, 3.4, 3.2, 3, 4.4, 3.5, 2, 3, 4, 1.5, 2, 2.2, 3, 4.5, 3.8, 2.7, 2.5]
N = len(xn)
n = np.arange(N)
f = interpolate.interp1d(n, xn, kind='cubic')
Ni = 4000
ni = np.linspace(0, N-1, Ni)
xni = f(ni)

ni = np.arange(Ni)

# parametros
# muestra correspondiente al tiempo 0 en la señal interpolada
s0 = 2100
# periodo de muestreo:
T = 90

# muestras en nT
n1 = np.arange(-math.floor(s0 / T), math.floor((Ni - s0) / T) + 1)
s1 = s0 + n1 * T

ds = 500
xmin = s0 - ds
xmax = s0 + ds
dx = T / 2
xmin_ax = xmin - dx
xmax_ax = xmax + dx

xni -= 0.5
ymax_ax = np.amax(xni) + 0.2
ymin_ax = -0.8

nn = ni[xmin: xmax]
xx = xni[xmin: xmax]

baseline = -0.55
fontsize = 10
fontsize2 = 12
grey = [0.6, 0.6, 0.6]
display_length = 6

fig = plt.figure(0, figsize=(10, 7), frameon=False)
#
# Señal discreta
#
ax = plt.subplot2grid((6, 4), (0, 0), rowspan=2, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.plot(nn, xx, '-', lw=1, color=grey)
ss = np.where((s1 > xmin) & (s1 < xmax))[0]
(markers, stemlines, bl) = plt.stem(s1[ss], xni[s1[ss]], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)

nmax = ds // T
for i in np.arange(-nmax, nmax + 1):
    plt.text(s0 + i * T, baseline, '${}$'.format(i), fontsize=fontsize, ha='center', va='baseline')

plt.text(xmax_ax, baseline, '$n$', fontsize=fontsize, ha='center', va='baseline')

ax.legend(['$x_c(t)$', '$x[n]=x_c(nT)$'], loc=1, frameon=False, fontsize=fontsize2, framealpha=1)

plt.text(xmin_ax, ymax_ax, r'$\textrm{Secuencia discreta}$', fontsize=fontsize2, ha='left', va='top')

plt.axis('off')

#
# Tren de impulsos
#
ax = plt.subplot2grid((6, 4), (2, 0), rowspan=2, colspan=4)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

plt.plot(nn, xx, '-', lw=1, color=grey)
for si in ss:
    plt.annotate(s='', xytext=(s1[si], 0), xy=(s1[si], xni[s1[si]]),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', facecolor='black',
                                 shrinkA=0, shrinkB=0))
for i in np.arange(-nmax, nmax + 1):
    if i == 0:
        plt.text(s0, baseline, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.abs(i) == 1:
        sg = '-' if i == -1 else ''
        plt.text(s0 + i * T, baseline, '${}T$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(s0 + i * T, baseline, '${}T$'.format(i), fontsize=fontsize, ha='center', va='baseline')

plt.text(s0, baseline, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, baseline, '$t$', fontsize=fontsize, ha='center', va='baseline')

custom_lines = [Line2D([0], [0], color=grey, lw=1),
                Line2D([0], [0], color='k', lw=1)]
ax.legend(custom_lines, ['$x_c(t)$', '$x_s(t)$'], loc=1, frameon=False, fontsize=fontsize, framealpha=1)

plt.text(xmin_ax, ymax_ax, r'$\textrm{Conversi\'on de secuencia a tren de impulsos}$', fontsize=fontsize2,
         ha='left', va='top')

plt.axis('off')

#
# Reconstrucción
#
ax = plt.subplot2grid((6, 4), (4, 0), rowspan=2, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
ss = np.where((s1 > xmin) & (s1 < xmax))[0]
plt.plot(s1[ss], xni[s1[ss]], 'ks', markersize=4, label='$x[n]$')
plt.plot(nn, xx, 'k-', lw=1, label='$x_r(t)$')

baseline2 = -0.8
for i in np.arange(-nmax, nmax + 1):
    sinc = xni[s1[ss]][i + nmax] * np.sin(np.pi * (nn - s0 - i * T) / T) / (np.pi * (nn - s0 - i * T) / T)
    plt.plot(nn, sinc, color=grey, lw=1, zorder=-1)
    plt.plot([s0 + i * T, s0 + i * T], [0, xtl], 'k-', lw=1)
    if i == 0:
        plt.text(s0, baseline2, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.abs(i) == 1:
        sg = '-' if i == -1 else ''
        plt.text(s0 + i * T, baseline2, '${}T$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(s0 + i * T, baseline2, '${}T$'.format(i), fontsize=fontsize, ha='center', va='baseline')

i = 1
plt.text(s1[ss][i + nmax], xni[s1[ss]][i + nmax]+0.2, '$x[1]$', fontsize=fontsize, ha='center', va='bottom')
sinc = xni[s1[ss]][i + nmax] * np.sin(np.pi * (nn - s0 - i * T) / T) / (np.pi * (nn - s0 - i * T) / T)
plt.plot(nn, sinc, color='r', lw=1, zorder=-1)
nn = 1.35 * T
sinc = xni[s1[ss]][i + nmax] * np.sin(np.pi * (nn - i * T) / T) / (np.pi * (nn - i * T) / T)
plt.annotate('$x[1]\dfrac{\sin[\pi(t-T)/T]}{\pi(t-T)/T}$', xytext=(s0 + 2.2 * T, xni[s1[ss]][i + nmax]),
             xycoords='data', xy=(s0 + nn, sinc),
             textcoords='data', fontsize=fontsize, va="center", ha="left",
             arrowprops=dict(arrowstyle="-|>, head_width=0.2, head_length=0.6", facecolor='black', relpos=(0, 0.5),
                             patchA=None, patchB=None, shrinkA=4, shrinkB=0))

plt.text(xmax_ax, baseline, '$t$', fontsize=fontsize, ha='center', va='baseline')

ax.legend(loc=1, frameon=False, fontsize=fontsize2, framealpha=1)

plt.text(xmin_ax, ymax_ax, r'$\textrm{Tren de impulsos filtrado pasabajos}$', fontsize=fontsize2,
         ha='left', va='top')

plt.axis('off')

plt.savefig('sampling_dc_conversion_math_representation.pdf', bbox_inches='tight')

plt.show()
