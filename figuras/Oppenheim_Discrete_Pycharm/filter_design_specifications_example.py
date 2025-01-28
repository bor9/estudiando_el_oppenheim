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
delta_p = 0.15
delta_s = 0.1
wp = np.pi / 4
ws = np.pi / 3

gpass = -20 * np.log10(1 - delta_p)
gstop = -20 * np.log10(delta_s)

print('gpass = {:.4f}'.format(gpass))
print('gstop = {:.4f}'.format(gstop))

### Filtro elíptico
# ellipord(wp, ws, gpass, gstop, analog=False, fs=None)
Ne, Wn = signal.ellipord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
# ellip(N, rp, rs, Wn, btype='low', analog=False, output='ba', fs=None)
b, a = signal.ellip(Ne, gpass, gstop, Wn, btype='low', analog=False, output='ba', fs=2*np.pi)
w, He = signal.freqz(b, a, worN=512, whole=False, plot=None, fs=2*np.pi, include_nyquist=False)

### Filtro Chebyshev tipo I
Nc1, Wn = signal.cheb1ord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
b, a = signal.cheby1(Nc1, gpass, Wn, btype='low', analog=False, output='ba', fs=2*np.pi)
_, Hc1 = signal.freqz(b, a, worN=512, whole=False, plot=None, fs=2*np.pi, include_nyquist=False)

### Filtro Chebyshev tipo II
Nc2, Wn = signal.cheb2ord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
b, a = signal.cheby2(Nc2, gstop, Wn, btype='low', analog=False, output='ba', fs=2*np.pi)
_, Hc2 = signal.freqz(b, a, worN=512, whole=False, plot=None, fs=2*np.pi, include_nyquist=False)

### Filtro Butterworth
Nb, Wn = signal.buttord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
b, a = signal.butter(Nb, Wn, btype='low', analog=False, output='ba', fs=2*np.pi)
_, Hb = signal.freqz(b, a, worN=512, whole=False, plot=None, fs=2*np.pi, include_nyquist=False)
print(b)
print(a)


### Parámetros de la gráfica
fontsize = 12
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
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hb), label=r'$\textrm{{Butterworth}}\;(N={})$'.format(Nb))
plt.plot(w, np.absolute(Hc1), label=r'$\textrm{{Chebyshev tipo I}}\;(N={})$'.format(Nc1))
plt.plot(w, np.absolute(Hc2), label=r'$\textrm{{Chebyshev tipo II}}\;(N={})$'.format(Nc2))
plt.plot(w, np.absolute(He), label=r'$\textrm{{El\'iptico}}\;(N={})$'.format(Ne))
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
xtm2 = -0.11
plt.plot([wp, wp], [0, xtl], 'k-', lw=1)
plt.text(wp, xtm2, '$\dfrac{\pi}{4}$', fontsize=fontsize, ha='center', va='baseline')
plt.plot([ws, ws], [0, xtl], 'k-', lw=1)
plt.text(ws, xtm2, '$\dfrac{\pi}{3}$', fontsize=fontsize, ha='center', va='baseline')
plt.plot([np.pi, np.pi], [0, xtl], 'k-', lw=1)
plt.text(np.pi, xtm, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm, xtm, '$0$', fontsize=fontsize, ha='right', va='baseline')
# eje y
plt.plot([0, xtl], [1 - delta_p, 1 - delta_p], 'k-', lw=1)
plt.text(ytm, 1 - delta_p, '${:.2f}$'.format(1 - delta_p), fontsize=fontsize, ha='right', va='center')
plt.plot([0, xtl], [1, 1], 'k-', lw=1)
plt.text(ytm, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, xtl], [delta_s, delta_s], 'k-', lw=1)
plt.text(ytm, delta_s, '${:.1f}$'.format(delta_s), fontsize=fontsize, ha='right', va='center')
# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.text(0.06, ymax_ax, '$|H(e^{j\omega})|$', fontsize=fontsize, ha='left', va='center')

plt.legend(bbox_to_anchor=(0.95, 0.9), loc='upper right', frameon=False, fontsize=11, framealpha=1)
plt.axis('off')

plt.savefig('filter_design_specifications_example.pdf', bbox_inches='tight')
plt.show()





