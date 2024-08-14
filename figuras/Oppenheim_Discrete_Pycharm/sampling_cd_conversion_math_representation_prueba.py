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
# periodo de muestreo: T1=T y T2=2T
T = 90
M = 3

# muestras en nT1 y nT2
n1 = np.arange(-math.floor(s0 / T), math.floor((Ni - s0) / T) + 1)
s1 = s0 + n1 * T
n2 = np.arange(-math.floor(s0 / (M * T)), math.floor((Ni - s0) / (M * T)) + 1)
s2 = s0 + n2 * M * T

ds = 1000
xmin = s0 - ds
xmax = s0 + ds
dx = T
xmin_ax = xmin - dx
xmax_ax = xmax + dx

xni -= 0.5
ymax_ax = np.amax(xni) + 0.4
ymin_ax = -0.5

nn = ni[xmin: xmax]
xx = xni[xmin: xmax]

baseline = -0.38
fontsize = 10
fontsize2 = 12
grey = [0.6, 0.6, 0.6]

fig = plt.figure(1, figsize=(10, 6), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=2, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

plt.plot(nn, xx, '-', lw=1, color=grey)
ss = np.where((s1 > xmin) & (s1 < xmax))[0]
for si in ss:
    plt.annotate(s='', xytext=(s1[si], 0), xy=(s1[si], xni[s1[si]]),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', facecolor='black',
                                 shrinkA=0, shrinkB=0))
for i in np.arange(2, 11, 2):
    plt.text(s0 - i * T, baseline, '${}T$'.format(-i), fontsize=fontsize, ha='center', va='baseline')
    plt.text(s0 + i * T, baseline, '${}T$'.format(i), fontsize=fontsize, ha='center', va='baseline')

plt.text(s0, baseline, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, baseline, '$t$', fontsize=fontsize, ha='center', va='baseline')

plt.text(1800, 3.8, '$T=T_1$', fontsize=fontsize2, ha='center', va='center')

custom_lines = [Line2D([0], [0], color=grey, lw=1),
                Line2D([0], [0], color='k', lw=1)]
ax.legend(custom_lines, ['$x_c(t)$', '$x_s(t)$'], loc=1, frameon=False, fontsize=fontsize2, framealpha=1)

plt.axis('off')

#
# Señal discreta con periodo de muestreo T
#
ax = plt.subplot2grid((4, 4), (2, 0), rowspan=2, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

plt.plot(nn, xx, '-', lw=1, color=grey)
(markers, stemlines, bl) = plt.stem(s1[ss], xni[s1[ss]], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

for i in np.arange(-10, 11, 2):
    plt.text(s0 + i * T, baseline, '${}$'.format(i), fontsize=fontsize, ha='center', va='baseline')

plt.text(xmax_ax, baseline, '$n$', fontsize=fontsize, ha='center', va='baseline')

ax.legend(['$x_c(t)$', '$x[n]$'], loc=1, frameon=False, fontsize=fontsize2, framealpha=1)

plt.axis('off')

#
# Con período de muestreo 2T
#
ax = plt.subplot2grid((4, 4), (0, 2), rowspan=2, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

plt.plot(nn, xx, '-', lw=1, color=grey)

ss = np.where((s2 > xmin) & (s2 < xmax))[0]
for si in ss:
    plt.annotate(s='', xytext=(s2[si], 0), xy=(s2[si], xni[s2[si]]),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', facecolor='black',
                                 shrinkA=0, shrinkB=0))
for i in np.arange(0, 6):
    if i == 0:
        plt.text(s0, baseline, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif i == 1:
        plt.text(s0 - i * 2 * T, baseline, '$T$'.format(-i), fontsize=fontsize, ha='center', va='baseline')
        plt.text(s0 + i * 2 * T, baseline, '$-T$'.format(i), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(s0 - i * 2 * T, baseline, '${}T$'.format(-i), fontsize=fontsize, ha='center', va='baseline')
        plt.text(s0 + i * 2 * T, baseline, '${}T$'.format(i), fontsize=fontsize, ha='center', va='baseline')


plt.text(1800, 3.8, '$T=2T_1$', fontsize=fontsize2, ha='center', va='center')
plt.text(xmax_ax, baseline, '$t$', fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

#
# Señal discreta con periodo de muestreo 2T
#
ds = 3000
xmin = s0 - ds
xmax = s0 + ds
dx = 2 * T
xmin_ax = xmin - dx
xmax_ax = xmax + dx

nn2 = ni[xmin: xmax]
xx2 = xni[xmin: xmax]


ss = np.where((s2 > xmin) & (s2 < xmax))[0]

ax = plt.subplot2grid((4, 4), (2, 2), rowspan=2, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

plt.plot(nn, xx, '-', lw=1, color=grey)
plt.plot(nn2, xx2, ':', lw=1, color=grey)
(markers, stemlines, bl) = plt.stem(s2[ss], xni[s2[ss]], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

for i in np.arange(-10, 11, 2):
    plt.text(s0 + i * 2 * T, baseline, '${}$'.format(i), fontsize=fontsize, ha='center', va='baseline')

plt.text(xmax_ax, baseline, '$n$', fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

plt.savefig('sampling_cd_conversion_math_representation.pdf', bbox_inches='tight')

plt.show()
