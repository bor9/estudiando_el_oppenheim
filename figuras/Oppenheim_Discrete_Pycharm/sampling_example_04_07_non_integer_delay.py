import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
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
s0 = 2180
# periodo de muestreo:
T = 90
Delta = 0.5

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
ymin_ax = -0.1

nn = ni[xmin: xmax]
xc = xni[xmin: xmax]
DT = int(T * Delta)
yc = xni[xmin - DT: xmax - DT]


baseline = -0.35
fontsize = 10
fontsize2 = 12
grey = [0.6, 0.6, 0.6]
display_length = 6
ms = 5

fig = plt.figure(0, figsize=(8, 3.5), frameon=False)
#
# Señal discreta
#
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.plot(nn, xc, '-', lw=2, color=grey, label='$x_c(t)$')
ss = np.where((s1 > xmin) & (s1 < xmax))[0]
(markers, stemlines, bl) = plt.stem(s1[ss], xni[s1[ss]], linefmt='k', markerfmt='sk', use_line_collection=True,
                                    label='$x[n]=x_c(nT)$')
plt.setp(markers, markersize=ms, zorder=10)
plt.setp(bl, visible=False)

plt.plot(nn, yc, '-', lw=2, color='salmon', label='$y_c(t)=x_c(t-T\Delta)$')
ss = np.where((s1 > xmin) & (s1 < xmax))[0]
(markers, stemlines, bl) = plt.stem(s1[ss], xni[s1[ss] - DT], linefmt='r', markerfmt='sr', use_line_collection=True,
                                    label='$y[n]=y_c(nT)$')
plt.setp(markers, markersize=ms, zorder=10)
plt.setp(bl, visible=False)

nmax = ds // T
for i in np.arange(-nmax, nmax + 1):
    plt.text(s0 + i * T, baseline, '${}$'.format(i), fontsize=fontsize2, ha='center', va='baseline')

plt.text(xmax_ax, baseline, '$n$', fontsize=fontsize2, ha='center', va='baseline')

handles, labels = plt.gca().get_legend_handles_labels()
order = [2, 0, 1, 3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], frameon=False, fontsize=fontsize2,
           framealpha=1, loc=1)

plt.axis('off')

plt.savefig('sampling_example_04_07_non_integer_delay.pdf', bbox_inches='tight')

plt.show()
