import matplotlib.pyplot as plt
import numpy as np
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
delta_p = 0.1
delta_s = 0.05
wp = 0.25 * np.pi
ws = 0.4 * np.pi
# Conversión a decibeles para uso en las función ellipord
gpass = -20 * np.log10(1 - delta_p)
gstop = -20 * np.log10(delta_s)

### Filtro elíptico
# ellipord(wp, ws, gpass, gstop, analog=False, fs=None)
N, Wn = signal.ellipord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
# ellip(N, rp, rs, Wn, btype='low', analog=False, output='ba', fs=None)
b, a = signal.ellip(N, gpass, gstop, Wn, btype='low', analog=False, output='ba', fs=2*np.pi)
w, He = signal.freqz(b, a, worN=512, whole=False, plot=None, fs=2*np.pi, include_nyquist=False)


### Parámetros de la gráfica
fontsize = 12
fontsize2 = 11
# x y ticks labels margin
xtm = -0.07
ytm = -0.04
display_length = 6

xmin = 0
xmax = np.pi
dx = 0.15
xmax_ax = xmax + dx
xmin_ax = xmin - dx
ymin_ax = -0.1
ymax_ax = 1.2

grey = [0.9, 0.9, 0.9]

fig = plt.figure(0, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
# magnitud de la respuesta en frecuencia
plt.plot(w, np.absolute(He), 'r-', lw=1.5)
# mascara
plt.plot([0, ws], [1, 1], 'k-')
plt.plot([0, wp], [1 - delta_p, 1 - delta_p], 'k-')
plt.plot([ws, np.pi], [delta_s, delta_s], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p], 'k-')
plt.plot([ws, ws], [delta_s, 1], 'k-')
# región pintada
vert = np.vstack(([0, wp, wp, 0], [0, 0, 1 - delta_p, 1 - delta_p]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [1, 1, mask_max, mask_max]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [delta_s, delta_s, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)

# etiquetas
# eje x
plt.plot([wp, wp], [0, xtl], 'k-', lw=1)
plt.text(wp, xtm, '$\omega_p$', fontsize=fontsize, ha='center', va='baseline')
plt.plot([ws, ws], [0, xtl], 'k-', lw=1)
plt.text(ws, xtm, '$\omega_s$', fontsize=fontsize, ha='center', va='baseline')
plt.plot([np.pi, np.pi], [0, xtl], 'k-', lw=1)
plt.text(np.pi, xtm, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm, xtm, '$0$', fontsize=fontsize, ha='right', va='baseline')
# eje y
plt.plot([0, xtl], [1 - delta_p, 1 - delta_p], 'k-', lw=1)
plt.text(ytm, 1 - delta_p, '$1-\delta_{p_2}$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, xtl], [1, 1], 'k-', lw=1)
plt.text(ytm, 1, '$1+\delta_{p_1}$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, xtl], [delta_s, delta_s], 'k-', lw=1)
plt.text(ytm, delta_s, '$\delta_s$', fontsize=fontsize, ha='right', va='center')
# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.text(0.06, ymax_ax, '$|H(e^{j\omega})|$', fontsize=fontsize, ha='left', va='center')
#
plt.text(wp / 2, 0.5, '$\mathrm{Banda}$\n$\mathrm{pasante}$', fontsize=fontsize2, ha='center', va='center')
plt.text((xmax + ws) / 2, 0.5, '$\mathrm{Banda}$\n$\mathrm{suprimida}$', fontsize=fontsize2, ha='center', va='center')
plt.text((wp + ws) / 2, 1 - delta_p / 2, '$\mathrm{Banda\,de}$\n$\mathrm{transici\\acute{o}n}$',
         fontsize=fontsize2, ha='center', va='top')


plt.axis('off')

plt.savefig('filter_design_specifications.pdf', bbox_inches='tight')
plt.show()
