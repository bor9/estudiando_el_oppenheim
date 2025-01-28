import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import numpy.polynomial.polynomial as poly
from matplotlib.patches import Polygon
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)


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
wp = 0.4 * np.pi
ws = 0.6 * np.pi
delta_p1 = 0.01
delta_p2 = 0.01
delta_s = 0.001

# normalización de las especificaciones para adaptarlas a las funciones de diseño de filtros IIR
delta_p2p = (delta_p1 + delta_p1) / (1 + delta_p1)
delta_sp = delta_s / (1 + delta_p1)

# Atenuación máxima en banda pasante y atenuación mínima en banda atenuada (dB)
gpass = -20 * np.log10(1 - delta_p2p)
gstop = -20 * np.log10(delta_sp)

print('delta_p2p = {:.6f}'.format(delta_p2p))
print('delta_sp = {:.6f}'.format(delta_sp))

### Filtro Butterworth
Nb, Wn = signal.buttord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
zb, pb, kb = signal.butter(Nb, Wn, btype='low', analog=False, output='zpk', fs=2*np.pi)
w, Hb = signal.freqz_zpk(zb, pb, kb, worN=512, whole=False, fs=2*np.pi)
Hb *= (1 + delta_p1)

### Filtro Chebyshev tipo I
Nc1, Wn = signal.cheb1ord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
zc1, pc1, kc1 = signal.cheby1(Nc1, gpass, Wn, btype='low', analog=False, output='zpk', fs=2*np.pi)
_, Hc1 = signal.freqz_zpk(zc1, pc1, kc1, worN=512, whole=False, fs=2*np.pi)
Hc1 *= (1 + delta_p1)

### Filtro Chebyshev tipo II
Nc2, Wn = signal.cheb2ord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
zc2, pc2, kc2 = signal.cheby2(Nc2, gstop, Wn, btype='low', analog=False, output='zpk', fs=2*np.pi)
_, Hc2 = signal.freqz_zpk(zc2, pc2, kc2, worN=512, whole=False, fs=2*np.pi)
Hc2 *= (1 + delta_p1)

### Filtro elíptico
Ne, Wn = signal.ellipord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
ze, pe, ke = signal.ellip(Ne, gpass, gstop, Wn, btype='low', analog=False, output='zpk', fs=2*np.pi)
_, He = signal.freqz_zpk(ze, pe, ke, worN=512, whole=False, fs=2*np.pi)
He *= (1 + delta_p1)

### Parámetros de la gráfica
fontsize = 11
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
plt.plot([0, ws], [1 + delta_p1, 1 + delta_p1], 'k-')
plt.plot([0, wp], [1 - delta_p2, 1 - delta_p2], 'k-')
plt.plot([ws, np.pi], [delta_s, delta_s], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p2], 'k-')
plt.plot([ws, ws], [delta_s, 1 + delta_p1], 'k-')
# región pintada
vert = np.vstack(([0, wp, wp, 0], [0, 0, 1 - delta_p2, 1 - delta_p2]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [1 + delta_p1, 1 + delta_p1, mask_max, mask_max]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [delta_s, delta_s, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)

# etiquetas
# eje x
xtm2 = -0.11
plt.plot([wp, wp], [0, xtl], 'k-', lw=1)
plt.text(wp, xtm, '${:.1f}\pi$'.format(wp / np.pi), fontsize=fontsize, ha='center', va='baseline')
plt.plot([ws, ws], [0, xtl], 'k-', lw=1)
plt.text(ws, xtm, '${:.1f}\pi$'.format(ws / np.pi), fontsize=fontsize, ha='center', va='baseline')
plt.plot([np.pi, np.pi], [0, xtl], 'k-', lw=1)
plt.text(np.pi, xtm, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm, xtm, '$0$', fontsize=fontsize, ha='right', va='baseline')
# eje y
plt.plot([0, xtl], [1 - delta_p2, 1 - delta_p2], 'k-', lw=1)
plt.text(ytm, 1 - delta_p2, '${:.2f}$'.format(1 - delta_p2), fontsize=fontsize, ha='right', va='top')
plt.plot([0, xtl], [1 + delta_p1, 1 + delta_p1], 'k-', lw=1)
plt.text(ytm, 1 + delta_p1, '${:.2f}$'.format(1 + delta_p1), fontsize=fontsize, ha='right', va='bottom')
plt.plot([0, xtl], [delta_s, delta_s], 'k-', lw=1)
plt.text(ytm, delta_s, '${:.3f}$'.format(delta_s), fontsize=fontsize, ha='right', va='bottom')
# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.text(0.06, ymax_ax, '$|H(e^{j\omega})|$', fontsize=fontsize, ha='left', va='center')

plt.legend(bbox_to_anchor=(0.035, 0.3), loc='lower left', frameon=False, fontsize=fontsize, framealpha=1)
plt.axis('off')
#plt.savefig('filter_design_example_07_05_specifications.pdf', bbox_inches='tight')

###############################################

# x y ticks labels margin
display_length = 6

xmin = 0
xmax = wp + 0.05
dx = 0
xmax_ax = xmax
xmin_ax = xmin
ymin_ax = 0.96
ymax_ax = 1.02

xticks = np.arange(0, xmax_ax / np.pi, 0.1)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[0] = '$0$'
xticks *= np.pi

fig = plt.figure(1, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hb), label=r'$\textrm{{Butterworth}}\;(N={})$'.format(Nb))
plt.plot(w, np.absolute(Hc1), label=r'$\textrm{{Chebyshev tipo I}}\;(N={})$'.format(Nc1))
plt.plot(w, np.absolute(Hc2), label=r'$\textrm{{Chebyshev tipo II}}\;(N={})$'.format(Nc2))
plt.plot(w, np.absolute(He), label=r'$\textrm{{El\'iptico}}\;(N={})$'.format(Ne))
# mascara
plt.plot([0, ws], [1 + delta_p1, 1 + delta_p1], 'k-')
plt.plot([0, wp], [1 - delta_p2, 1 - delta_p2], 'k-')
plt.plot([ws, np.pi], [delta_s, delta_s], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p2], 'k-')
plt.plot([ws, ws], [delta_s, 1 + delta_p1], 'k-')
# región pintada
vert = np.vstack(([0, wp, wp, 0], [0, 0, 1 - delta_p2, 1 - delta_p2]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [1 + delta_p1, 1 + delta_p1, mask_max, mask_max]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [delta_s, delta_s, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)

# etiquetas
plt.xticks(xticks, xticks_labels, usetex=True)
# etiquetas de los ejes
plt.xlabel('$\omega$', fontsize=fontsize)
plt.ylabel('$|H(e^{j\omega})|$', fontsize=fontsize)
plt.legend(bbox_to_anchor=(0.035, 0.06), loc='lower left', frameon=False, fontsize=fontsize, framealpha=1)

xmin = ws - 0.05
xmax = np.pi
xmax_ax = xmax
xmin_ax = xmin
ymin_ax = 0
ymax_ax = 0.003

xticks = np.arange(round(ws / np.pi, 1), 1.01, 0.1)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[-1] = '$\pi$'
xticks *= np.pi
yticks = np.linspace(ymin_ax, ymax_ax, int(ymax_ax/0.001 + 1))
yticks_labels = ['${:.3f}$'.format(yt) for yt in yticks]

ax = plt.subplot2grid((4, 4), (0, 2), rowspan=4, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hb), label=r'$\textrm{{Butterworth}}\;(N={})$'.format(Nb))
plt.plot(w, np.absolute(Hc1), label=r'$\textrm{{Chebyshev tipo I}}\;(N={})$'.format(Nc1))
plt.plot(w, np.absolute(Hc2), label=r'$\textrm{{Chebyshev tipo II}}\;(N={})$'.format(Nc2))
plt.plot(w, np.absolute(He), label=r'$\textrm{{El\'iptico}}\;(N={})$'.format(Ne))
# mascara
plt.plot([0, ws], [1 + delta_p1, 1 + delta_p1], 'k-')
plt.plot([0, wp], [1 - delta_p2, 1 - delta_p2], 'k-')
plt.plot([ws, np.pi], [delta_s, delta_s], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p2], 'k-')
plt.plot([ws, ws], [delta_s, 1 + delta_p1], 'k-')
# región pintada
vert = np.vstack(([0, wp, wp, 0], [0, 0, 1 - delta_p2, 1 - delta_p2]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [1 + delta_p1, 1 + delta_p1, mask_max, mask_max]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [delta_s, delta_s, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)

# etiquetas
plt.xticks(xticks, xticks_labels, usetex=True)
plt.yticks(yticks, yticks_labels, usetex=True)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
# etiquetas de los ejes
plt.xlabel('$\omega$', fontsize=fontsize)
plt.ylabel('$|H(e^{j\omega})|$', fontsize=fontsize)

plt.savefig('filter_design_example_07_05_specifications.pdf', bbox_inches='tight')


###############################################


#
# Diagrama de polos y ceros
#
xmax = 1.5
xmin = -xmax
ymax = 1.5
ymin = -1.2
ymax_ax = 2
ytm = 0.08
xtm = -0.24

fs2 = 10

# circulo unidad
theta = np.linspace(0, 2 * np.pi, 100)
xc = np.cos(theta)
yc = np.sin(theta)

zeros_marker_style = dict(marker='o', linestyle='', markersize=6, markerfacecolor='w', markeredgecolor='k',
                          markeredgewidth=1.5)
poles_marker_style = dict(marker='x', linestyle='', markersize=6, markeredgecolor='k', markeredgewidth=2)


fig = plt.figure(2, figsize=(8, 7), frameon=False)
# Butterworth
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=2, colspan=2)
plt.ylim(ymin, ymax_ax)
plt.xlim(xmin, xmax)
plt.gca().set_aspect('equal', adjustable='box')
# ejes
# axis arrows
plt.annotate("", xytext=(xmin, 0), xycoords='data', xy=(xmax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin), xycoords='data', xy=(0, ymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
# circulo unidad
ax.plot(xc, yc, 'k-', lw=1)
# polos y ceros
plt.plot(zb.real, zb.imag, **zeros_marker_style, zorder=10)
plt.plot(pb.real, pb.imag, **poles_marker_style)
# multiplicidad del cero
plt.text(zb[0].real-0.05, zb[0].imag+0.05, '{:d}'.format(Nb), fontsize=fs2, ha='right', va='bottom')

# etiquetas de los ejes
plt.text(xmax, xtm, r'$\textrm{Re}(z)$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm, ymax, r'$\textrm{Im}(z)$', fontsize=fontsize, ha='left', va='center')

plt.text(0, 1.85, r'$\textrm{{Butterworth}}\;(N={})$'.format(Nb),
         fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

# Cheby 1
ax = plt.subplot2grid((4, 4), (0, 2), rowspan=2, colspan=2)
plt.ylim(ymin, ymax_ax)
plt.xlim(xmin, xmax)
plt.gca().set_aspect('equal', adjustable='box')
# ejes
# axis arrows
plt.annotate("", xytext=(xmin, 0), xycoords='data', xy=(xmax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin), xycoords='data', xy=(0, ymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
# circulo unidad
ax.plot(xc, yc, 'k-', lw=1)
# polos y ceros
plt.plot(zc1.real, zc1.imag, **zeros_marker_style, zorder=10)
plt.plot(pc1.real, pc1.imag, **poles_marker_style)
# multiplicidad del cero
plt.text(zb[0].real-0.05, zb[0].imag+0.05, '{:d}'.format(Nc1), fontsize=fs2, ha='right', va='bottom')

# etiquetas de los ejes
plt.text(xmax, xtm, r'$\textrm{Re}(z)$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm, ymax, r'$\textrm{Im}(z)$', fontsize=fontsize, ha='left', va='center')

plt.text(0, 1.85, r'$\textrm{{Chebyshev tipo I}}\;(N={})$'.format(Nc1),
         fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

# Cheby 2
ax = plt.subplot2grid((4, 4), (2, 0), rowspan=2, colspan=2)
plt.ylim(ymin, ymax_ax)
plt.xlim(xmin, xmax)
plt.gca().set_aspect('equal', adjustable='box')
# ejes
# axis arrows
plt.annotate("", xytext=(xmin, 0), xycoords='data', xy=(xmax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin), xycoords='data', xy=(0, ymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
# circulo unidad
ax.plot(xc, yc, 'k-', lw=1)
# polos y ceros
plt.plot(zc2.real, zc2.imag, **zeros_marker_style, zorder=10)
plt.plot(pc2.real, pc2.imag, **poles_marker_style)

# etiquetas de los ejes
plt.text(xmax, xtm, r'$\textrm{Re}(z)$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm, ymax, r'$\textrm{Im}(z)$', fontsize=fontsize, ha='left', va='center')

plt.text(0, 1.85, r'$\textrm{{Chebyshev tipo II}}\;(N={})$'.format(Nc2),
         fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

# Elíptico
ax = plt.subplot2grid((4, 4), (2, 2), rowspan=2, colspan=2)
plt.ylim(ymin, ymax_ax)
plt.xlim(xmin, xmax)
plt.gca().set_aspect('equal', adjustable='box')
# ejes
# axis arrows
plt.annotate("", xytext=(xmin, 0), xycoords='data', xy=(xmax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin), xycoords='data', xy=(0, ymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
# circulo unidad
ax.plot(xc, yc, 'k-', lw=1)
# polos y ceros
plt.plot(ze.real, ze.imag, **zeros_marker_style, zorder=10)
plt.plot(pe.real, pe.imag, **poles_marker_style)

# etiquetas de los ejes
plt.text(xmax, xtm, r'$\textrm{Re}(z)$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm, ymax, r'$\textrm{Im}(z)$', fontsize=fontsize, ha='left', va='center')

plt.text(0, 1.85, r'$\textrm{{El\'iptico}}\;(N={})$'.format(Ne),
         fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

plt.savefig('filter_design_example_07_05_pole_zero.pdf', bbox_inches='tight')


# Magnitud de la función del sistema \H(z)\
xmin = -1.5
xmax = 1.5
ymin = -1.5
ymax = 1.5
xr = np.linspace(xmin, xmax, 200)
yr = np.linspace(ymin, ymax, 200)
X, Y = np.meshgrid(xr, yr)

b, a = signal.zpk2tf(zb, pb, kb)
Hb_mag = np.absolute(poly.polyval(X + 1j * Y, np.flipud(b)) / poly.polyval(X + 1j * Y, np.flipud(a)))
b, a = signal.zpk2tf(zc1, pc1, kc1)
Hc1_mag = np.absolute(poly.polyval(X + 1j * Y, np.flipud(b)) / poly.polyval(X + 1j * Y, np.flipud(a)))
b, a = signal.zpk2tf(zc2, pc2, kc2)
Hc2_mag = np.absolute(poly.polyval(X + 1j * Y, np.flipud(b)) / poly.polyval(X + 1j * Y, np.flipud(a)))
b, a = signal.zpk2tf(ze, pe, ke)
He_mag = np.absolute(poly.polyval(X + 1j * Y, np.flipud(b)) / poly.polyval(X + 1j * Y, np.flipud(a)))

zmax = 1.8
zmin = -1
xmax = 1.5
xmin = -xmax
ymax = 1.5
ymin = -ymax
xmax_ax = 2

levels_surf = np.linspace(0, zmax, 200)

xH = np.cos(w)
yH = np.sin(w)

xticks = [-1, 0, 1]
yticks = xticks
zticks = [0, 0.5, 1, 1.5]

fontsize = 12
fontsize2 = 9

# altura del titulos
htitle = 0.93

fig = plt.figure(3, figsize=(10, 9), frameon=False)
# Butterworth
ax = plt.subplot2grid((16, 16), (0, 0), rowspan=8, colspan=8, projection='3d')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(zmin, zmax)
plt.contour(X, Y, Hb_mag, levels=levels_surf, cmap=cm.winter, vmin=0, vmax=zmax, zorder=-100, linewidths = 1)
# circulo unidad
ax.plot(xc, yc, zs=zmin, zdir='z', color='k', lw=1)
ax.plot3D(xH, yH, np.absolute(Hb), color='r', zorder=500)
ax.plot3D(xH, -yH, np.absolute(Hb), color='r', zorder=500)
# polos y ceros
plt.plot(zb.real, zb.imag, **zeros_marker_style, zorder=10, zs=zmin, zdir='z')
plt.plot(pb.real, pb.imag, **poles_marker_style, zorder=10, zs=zmin, zdir='z')
ax.text(zb[0].real - 0.1, zb[0].imag - 0.2, zmin, '{:d}'.format(Nb), fontsize=fontsize2, ha='right', va='top')

# axis arrows
arw = Arrow3D([ymin, ymax], [0, 0], [zmin, zmin], arrowstyle="-|>, head_width=0.5, head_length=1",
              lw=1, mutation_scale=5, facecolor='black')
ax.add_artist(arw)
arw = Arrow3D([0, 0], [xmin, xmax], [zmin, zmin], arrowstyle="-|>, head_width=0.5, head_length=1",
              lw=1, mutation_scale=5, facecolor='black')
ax.add_artist(arw)

ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_zticks(zticks)

ax.view_init(elev=21, azim=-111)
#ax.view_init(elev=17, azim=-70)
ax.set_xlabel(r'$\textrm{Re}(z)$', fontsize=fontsize)
ax.set_ylabel(r'$\textrm{Im}(z)$', fontsize=fontsize)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel(r'$|H(z)|$', fontsize=fontsize, rotation=90, labelpad=2)

ax.text2D(0.5, htitle, r'$\textrm{{Butterworth}}\;(N={})$'.format(Nb), fontsize=fontsize,
          ha='center', va='baseline', transform=ax.transAxes)

# Cheby 1
ax = plt.subplot2grid((16, 16), (0, 8), rowspan=8, colspan=8, projection='3d')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(zmin, zmax)
plt.contour(X, Y, Hc1_mag, levels=levels_surf, cmap=cm.winter, vmin=0, vmax=zmax, zorder=-100, linewidths = 1)
# circulo unidad
ax.plot(xc, yc, zs=zmin, zdir='z', color='k', lw=1)
ax.plot3D(xH, yH, np.absolute(Hc1), color='r', zorder=500)
ax.plot3D(xH, -yH, np.absolute(Hc1), color='r', zorder=500)
# polos y ceros
plt.plot(zc1.real, zc1.imag, **zeros_marker_style, zorder=10, zs=zmin, zdir='z')
plt.plot(pc1.real, pc1.imag, **poles_marker_style, zorder=10, zs=zmin, zdir='z')
ax.text(zc1[0].real - 0.1, zc1[0].imag - 0.2, zmin, '{:d}'.format(Nc1), fontsize=fontsize2, ha='right', va='top')

# axis arrows
arw = Arrow3D([ymin, ymax], [0, 0], [zmin, zmin], arrowstyle="-|>, head_width=0.5, head_length=1",
              lw=1, mutation_scale=5, facecolor='black')
ax.add_artist(arw)
arw = Arrow3D([0, 0], [xmin, xmax], [zmin, zmin], arrowstyle="-|>, head_width=0.5, head_length=1",
              lw=1, mutation_scale=5, facecolor='black')
ax.add_artist(arw)

ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_zticks(zticks)

ax.view_init(elev=21, azim=-111)
#ax.view_init(elev=17, azim=-70)
ax.set_xlabel(r'$\textrm{Re}(z)$', fontsize=fontsize)
ax.set_ylabel(r'$\textrm{Im}(z)$', fontsize=fontsize)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel(r'$|H(z)|$', fontsize=fontsize, rotation=90, labelpad=2)

ax.text2D(0.5, htitle, r'$\textrm{{Chebyshev tipo I}}\;(N={})$'.format(Nc1), fontsize=fontsize,
          ha='center', va='baseline', transform=ax.transAxes)

# Cheby 2
ax = plt.subplot2grid((16, 16), (8, 0), rowspan=8, colspan=8, projection='3d')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(zmin, zmax)
plt.contour(X, Y, Hc2_mag, levels=levels_surf, cmap=cm.winter, vmin=0, vmax=zmax, zorder=-100, linewidths = 1)
# circulo unidad
ax.plot(xc, yc, zs=zmin, zdir='z', color='k', lw=1)
ax.plot3D(xH, yH, np.absolute(Hc2), color='r', zorder=500)
ax.plot3D(xH, -yH, np.absolute(Hc2), color='r', zorder=500)
# polos y ceros
plt.plot(zc2.real, zc2.imag, **zeros_marker_style, zorder=10, zs=zmin, zdir='z')
plt.plot(pc2.real, pc2.imag, **poles_marker_style, zorder=10, zs=zmin, zdir='z')

# axis arrows
arw = Arrow3D([ymin, ymax], [0, 0], [zmin, zmin], arrowstyle="-|>, head_width=0.5, head_length=1",
              lw=1, mutation_scale=5, facecolor='black')
ax.add_artist(arw)
arw = Arrow3D([0, 0], [xmin, xmax], [zmin, zmin], arrowstyle="-|>, head_width=0.5, head_length=1",
              lw=1, mutation_scale=5, facecolor='black')
ax.add_artist(arw)

ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_zticks(zticks)

ax.view_init(elev=21, azim=-111)
#ax.view_init(elev=17, azim=-70)
ax.set_xlabel(r'$\textrm{Re}(z)$', fontsize=fontsize)
ax.set_ylabel(r'$\textrm{Im}(z)$', fontsize=fontsize)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel(r'$|H(z)|$', fontsize=fontsize, rotation=90, labelpad=2)

ax.text2D(0.5, htitle, r'$\textrm{{Chebyshev tipo II}}\;(N={})$'.format(Nc2), fontsize=fontsize,
          ha='center', va='baseline', transform=ax.transAxes)

# Elíptico
ax = plt.subplot2grid((16, 16), (8, 8), rowspan=8, colspan=8, projection='3d')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(zmin, zmax)
plt.contour(X, Y, He_mag, levels=levels_surf, cmap=cm.winter, vmin=0, vmax=zmax, zorder=-100, linewidths = 1)
# circulo unidad
ax.plot(xc, yc, zs=zmin, zdir='z', color='k', lw=1)
ax.plot3D(xH, yH, np.absolute(He), color='r', zorder=500)
ax.plot3D(xH, -yH, np.absolute(He), color='r', zorder=500)
# polos y ceros
plt.plot(ze.real, ze.imag, **zeros_marker_style, zorder=10, zs=zmin, zdir='z')
plt.plot(pe.real, pe.imag, **poles_marker_style, zorder=10, zs=zmin, zdir='z')

# axis arrows
arw = Arrow3D([ymin, ymax], [0, 0], [zmin, zmin], arrowstyle="-|>, head_width=0.5, head_length=1",
              lw=1, mutation_scale=5, facecolor='black')
ax.add_artist(arw)
arw = Arrow3D([0, 0], [xmin, xmax], [zmin, zmin], arrowstyle="-|>, head_width=0.5, head_length=1",
              lw=1, mutation_scale=5, facecolor='black')
ax.add_artist(arw)

ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_zticks(zticks)

ax.view_init(elev=21, azim=-111)
#ax.view_init(elev=17, azim=-70)
ax.set_xlabel(r'$\textrm{Re}(z)$', fontsize=fontsize)
ax.set_ylabel(r'$\textrm{Im}(z)$', fontsize=fontsize)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel(r'$|H(z)|$', fontsize=fontsize, rotation=90, labelpad=2)

ax.text2D(0.5, htitle, r'$\textrm{{El\'iptico}}\;(N={})$'.format(Ne), fontsize=fontsize,
          ha='center', va='baseline', transform=ax.transAxes)

#plt.savefig('filter_design_example_07_05_Hz_magnitude.pdf', bbox_inches='tight')
plt.show()





