import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import polynomial as P
from scipy import special
from scipy import signal
from matplotlib.patches import Polygon

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


# Especificaciones del filtro
Amax = 0.1
Amin = 40
wp = 20
ws = 26

w = np.linspace(0, 60, 500)

### Filtro elíptico
# ellipord(wp, ws, gpass, gstop, analog=False, fs=None)
N, Wn = signal.ellipord(wp, ws, Amax, Amin, analog=True)
# ellip(N, rp, rs, Wn, btype='low', analog=False, output='ba', fs=None)
b, a = signal.ellip(N, Amax, Amin, Wn, btype='low', analog=True, output='ba')
_, H = signal.freqs(b, a, worN=w)

z, p, k = signal.ellip(N, Amax, Amin, Wn, btype='low', analog=True, output='zpk')

print('p = {}'.format(p))
print('z = {}'.format(z))


###################################

delta_p = 1 - 10 ** (-Amax / 20)
delta_s = 10 ** (-Amin / 20)

### Parámetros de la gráfica
fontsize = 12
# x y ticks labels margin
xtm = -0.07
ytm = -0.04
display_length = 6

xmin_ax = 0
xmax_ax = w[-1]
ymin_ax = 0
ymax_ax = 1.1

grey = [0.9, 0.9, 0.9]

fig = plt.figure(0, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(H), label=r'$\textrm{{El\'iptico (N={})}}$'.format(N))
# mascara
plt.plot([0, ws], [1, 1], 'k-')
plt.plot([0, wp], [1 - delta_p, 1 - delta_p], 'k-')
plt.plot([ws, w[-1]], [delta_s, delta_s], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p], 'k-')
plt.plot([ws, ws], [delta_s, 1], 'k-')
# región pintada
vert = np.vstack(([0, wp, wp, 0], [0, 0, 1 - delta_p, 1 - delta_p]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
mask_max = ymax_ax
vert = np.vstack(([0, ws, ws, 0], [1, 1, mask_max, mask_max]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax_ax, xmax_ax, ws], [delta_s, delta_s, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)
plt.legend(bbox_to_anchor=(0.95, 0.9), loc='upper right', frameon=False, fontsize=11, framealpha=1)
#plt.axis('off')

#plt.savefig('filter_design_specifications_example.pdf', bbox_inches='tight')

plt.show()
