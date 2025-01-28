import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from matplotlib.patches import Polygon
from matplotlib.patches import ConnectionPatch

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


### Parámetros de la gráfica
fontsize = 11
# x y ticks labels margin
xtm1 = -0.65
ytm1 = -0.02
xtm2 = -0.13
ytm2 = -0.06

display_length = 6

# Ejes de las gráficas
wmin = 0
wmax = np.pi
dx = 0.25
wmax_ax = wmax + dx
wmin_ax = wmin

Wmin = 0
Wmax = 3 * np.pi
dx = 0.15
Wmax_ax = Wmax + dx
Wmin_ax = Wmin - 0.8

ymin = 0
ymin_ax = -0.2
ymax = 1.1
ymax_ax = 1.2

# Conversión de los requerimientos al tiempo continuo
Td = 1
Wp = (2 / Td) * np.tan(wp / 2)
Ws = (2 / Td) * np.tan(ws / 2)

dW = 0.3
wmax2 = 2 * np.arctan((Wmax - dW) * Td / 2)
w = np.linspace(0, wmax2, 200)
W = (2 / Td) * np.tan(w / 2)

### Filtro elíptico analógico
N, Wn = signal.ellipord(Wp, Ws, gpass, gstop, analog=True, fs=None)
b, a = signal.ellip(N, gpass, gstop, Wn, btype='low', analog=True, output='ba', fs=None)
_, He = signal.freqs(b, a, worN=W, plot=None)

grey = [0.9, 0.9, 0.9]

fig = plt.figure(0, figsize=(8, 6), frameon=False)

ax1 = plt.subplot2grid((8, 8), (0, 0), rowspan=5, colspan=4)
plt.xlim(ymin_ax, ymax_ax)
plt.ylim(Wmin_ax, Wmax_ax)
xtl, ytl = convert_display_to_data_coordinates(ax1.transData, length=display_length)
plt.annotate("", xytext=(ymin, 0), xycoords='data', xy=(ymax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, Wmin), xycoords='data', xy=(0, Wmax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.gca().invert_xaxis()
plt.plot(np.absolute(He), W, 'r-', lw=1.5)
plt.plot([1, 1], [0, Ws], 'k-', lw=1)
plt.plot([1 - delta_p, 1 - delta_p], [0, Wp], 'k-', lw=1)
plt.plot([delta_s, delta_s], [Ws, Wmax - dW], 'k-', lw=1)
# región pintada
vert = np.vstack(([0, 0, 1 - delta_p, 1 - delta_p], [0, Wp, Wp, 0]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-100)
ax1.add_artist(p1)
mask_max = 1.1
vert = np.vstack(([1, 1, mask_max, mask_max], [0, Ws, Ws, 0]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-100)
ax1.add_artist(p2)
vert = np.vstack(([delta_s, delta_s, mask_max, mask_max], [Ws, Wmax - dW, Wmax - dW, Ws]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-100)
ax1.add_artist(p3)
# etiquetas
plt.text(1 - delta_p - 0.05, xtm1, '$1-\delta_p$', fontsize=fontsize, ha='center', va='baseline')
plt.text(1, xtm1, '$1$', fontsize=fontsize, ha='center', va='baseline')
plt.plot([delta_s, delta_s], [0, xtl], 'k-', lw=1)
plt.text(delta_s, xtm1, '$\delta_s$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ymax_ax, xtm1, '$|H_c(j\Omega)|$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm1, Wp - 0.15, '$\Omega_p$', fontsize=fontsize, ha='left', va='top')
plt.text(ytm1, Ws + 0.15, '$\Omega_s$', fontsize=fontsize, ha='left', va='bottom')
plt.text(-0.035, Wmax_ax, '$\Omega$', fontsize=fontsize, ha='left', va='top')
plt.text(ytm1, xtm1, '$0$', fontsize=fontsize, ha='left', va='baseline')
plt.text(0.9, 2 * Wmax_ax / 3, r'$\Omega_p=\dfrac{2}{T_d}\tan\left(\dfrac{\omega_p}{2}\right)$',
         fontsize=fontsize, ha='left', va='baseline')
plt.text(0.9, 2 * Wmax_ax / 3 * 0.78, r'$\Omega_s=\dfrac{2}{T_d}\tan\left(\dfrac{\omega_s}{2}\right)$',
         fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

ax2 = plt.subplot2grid((8, 8), (0, 4), rowspan=5, colspan=4)
plt.xlim(wmin_ax, wmax_ax)
plt.ylim(Wmin_ax, Wmax_ax)
xtl, ytl = convert_display_to_data_coordinates(ax2.transData, length=display_length)
plt.annotate("", xytext=(wmin, 0), xycoords='data', xy=(wmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, Wmin), xycoords='data', xy=(0, Wmax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(w, W, 'k-', lw=1.5)
# etiquetas
plt.text(wp - 0.03, xtm1, '$\omega_p$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ws + 0.03, xtm1, '$\omega_s$', fontsize=fontsize, ha='left', va='baseline')
plt.plot([np.pi, np.pi], [0, xtl], 'k-', lw=1)
plt.text(np.pi, xtm1, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.text(wmax_ax, xtm1, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm1, Wp - 0.15, '$\Omega_p$', fontsize=fontsize, ha='right', va='top')
plt.text(ytm1, Ws + 0.15, '$\Omega_s$', fontsize=fontsize, ha='right', va='bottom')
plt.text(-0.07, Wmax_ax, '$\Omega$', fontsize=fontsize, ha='right', va='top')
plt.text(ytm1, xtm1, '$0$', fontsize=fontsize, ha='right', va='baseline')
ww = 2 * np.pi / 3
WW = (2 / Td) * np.tan(ww / 2)
plt.text(ww, WW, r'$\Omega=\dfrac{2}{T_d}\tan\left(\dfrac{\omega}{2}\right)$',
         fontsize=fontsize, ha='left', va='top')
plt.axis('off')

ax3 = plt.subplot2grid((8, 8), (5, 4), rowspan=3, colspan=4)
plt.xlim(wmin_ax, wmax_ax)
plt.ylim(ymin, ymax_ax)
xtl, ytl = convert_display_to_data_coordinates(ax3.transData, length=display_length)
plt.annotate("", xytext=(wmin, 0), xycoords='data', xy=(wmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(w, np.absolute(He), 'r-', lw=1.5)
plt.plot([0, ws], [1, 1], 'k-', lw=1)
plt.plot([0, wp], [1 - delta_p, 1 - delta_p], 'k-', lw=1)
plt.plot([ws, np.pi], [delta_s, delta_s], 'k-', lw=1)
# región pintada
vert = np.vstack(([0, wp, wp, 0], [0, 0, 1 - delta_p, 1 - delta_p]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-100)
ax3.add_artist(p1)
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [1, 1, mask_max, mask_max]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-100)
ax3.add_artist(p2)
vert = np.vstack(([ws, wmax, wmax, ws], [delta_s, delta_s, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-100)
ax3.add_artist(p3)
# etiquetas
plt.text(wp - 0.03, xtm2, '$\omega_p$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ws + 0.03, xtm2, '$\omega_s$', fontsize=fontsize, ha='center', va='baseline')
plt.plot([np.pi, np.pi], [0, xtl], 'k-', lw=1)
plt.text(np.pi, xtm2, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm1, xtm2, '$0$', fontsize=fontsize, ha='right', va='baseline')
plt.text(wmax_ax, xtm2, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm2, 1 - delta_p - 0.02, '$1-\delta_p$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm2, delta_s, '$\delta_s$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, ytl], [delta_s, delta_s], 'k-', lw=1)
plt.text(ytm2, ymax_ax, '$|H(e^{j\omega})|$', fontsize=fontsize, ha='right', va='top')

plt.axis('off')

con = ConnectionPatch(xyB=(wp, Wp), xyA=(wp, 0), coordsA="data", coordsB="data",
                      axesB=ax2, axesA=ax3, color='k', lw=1, linestyle='dashed', zorder=100)
ax3.add_artist(con)
con = ConnectionPatch(xyB=(ws, Ws), xyA=(ws, 0), coordsA="data", coordsB="data",
                      axesB=ax2, axesA=ax3, color='k', lw=1, linestyle='dashed', zorder=100)
ax3.add_artist(con)
con = ConnectionPatch(xyA=(1 - delta_p, Wp), xyB=(wp, Wp), coordsA="data", coordsB="data",
                      axesA=ax1, axesB=ax2, color='k', lw=1, linestyle='dashed', zorder=100)
ax1.add_artist(con)
con = ConnectionPatch(xyA=(1, Ws), xyB=(ws, Ws), coordsA="data", coordsB="data",
                      axesA=ax1, axesB=ax2, color='k', lw=1, linestyle='dashed', zorder=100)
ax1.add_artist(con)

plt.savefig('filter_design_bilinear_specifications_mapping.pdf', bbox_inches='tight')

plt.show()
