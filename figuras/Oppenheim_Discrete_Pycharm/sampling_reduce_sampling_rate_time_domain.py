import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import math

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
# xni: señal, ni: indice de muestras 0: Ni - 1
xni = f(ni)
# se reduce la amplitud de las muestras
xni -= 0.5
ni = np.arange(Ni)

# parametros
# muestra correspondiente al tiempo 0 en la señal interpolada
s0 = 2200
# periodo de muestreo: T1=T y T2=MT
T = 50
# factor de reducción de la tasa de muestreo
M = 3

# s1 y s2 son los indices de las muestras en nT y nMT en la señal xni
n1 = np.arange(-math.floor(s0 / T), math.floor((Ni - s0) / T))
s1 = s0 + n1 * T
n2 = np.arange(-math.floor(s0 / (M * T)), math.floor((Ni - s0) / (M * T)))
s2 = s0 + n2 * M * T


# rango de muestras en torno a la muestra s0 de la señal continua
ds = 640
xmin = s0 - ds
xmax = s0 + ds
# señal x[n]
# indice absoluto de las muestras en el rango de interés [xmin, xmax]
ss1 = np.where((s1 > xmin) & (s1 < xmax))[0]
# indice absoluto y valor de las muestras
na1 = s1[ss1]
x1 = xni[s1[ss1]]
# indice relativo de las muestras
n1 = n1[ss1]
# señal x[nM]
# indice absoluto de las muestras en el rango de interés [xmin, xmax]
ss2 = np.where((s2 > xmin) & (s2 < xmax))[0]
# indice absoluto y valor de las muestras
na2 = s2[ss2]
x2 = xni[s2[ss2]]
# indice relativo de las muestras
n2 = n2[ss2]

# rango de los ejes de las gráficas
dx = T
xmin_ax = xmin - dx
xmax_ax = xmax + dx

ymax_ax = np.amax(xni) + 0.4
ymin_ax = -0.5

# señal en tiempo continuo en el rango de interés
nn = ni[xmin: xmax]
xx = xni[xmin: xmax]

baseline = -0.38
fontsize = 10
fontsize2 = 12
grey = [0.6, 0.6, 0.6]

fig = plt.figure(1, figsize=(9, 6), frameon=False)
#
# Señal discreta con periodo de muestreo T
#
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=2, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.plot(nn, xx, '-', lw=1, color=grey)
# señal x[n]
(markers, stemlines, bl) = plt.stem(na1, x1, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
for i, ni in enumerate(n1):
    plt.text(na1[i], baseline, '${}$'.format(ni), fontsize=fontsize, ha='center', va='baseline')
# señal x[nM]
(markers, stemlines, bl) = plt.stem(na2, x2, linefmt='r', markerfmt='sr', use_line_collection=True)
plt.setp(markers, markersize=3.5, zorder=10)
plt.setp(bl, visible=False)

plt.text(xmax_ax, baseline, '$n$', fontsize=fontsize, ha='center', va='baseline')

plt.text(xmin_ax, ymax_ax, '$x[n]$', fontsize=fontsize2, ha='left', va='top')

plt.axis('off')

#
# Señal discreta con periodo de muestreo MT
#

# valores de los límites de los ejes: deben ser M veces mas grande que el anterior
# para mantener las proporciones de las gráficas
ds = M * ds
xmin = s0 - ds
xmax = s0 + ds
dx = M * T
xmin_ax = xmin - dx
xmax_ax = xmax + dx

ax = plt.subplot2grid((4, 4), (2, 0), rowspan=2, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.plot(nn, xx, '-', lw=1, color=grey)
(markers, stemlines, bl) = plt.stem(na2, x2, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# indice relativo de las muestras
for i, ni in enumerate(n2):
    plt.text(na2[i], baseline, '${}$'.format(ni), fontsize=fontsize, ha='center', va='baseline')

plt.text(xmax_ax, baseline, '$n$', fontsize=fontsize, ha='center', va='baseline')

plt.text(xmin_ax, ymax_ax, '$x_d[n]=x[nM]$\n$M=3$', fontsize=fontsize2, ha='left', va='top')

plt.axis('off')

plt.savefig('sampling_reduce_sampling_rate_time_domain.pdf', bbox_inches='tight')

plt.show()
