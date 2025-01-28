import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import numpy.polynomial.polynomial as poly
import math
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
# Frecuencias críticas (rad)
wp = 0.22 * np.pi
ws = 0.29 * np.pi
# Atenuación máxima en banda pasante y atenuación mínima en banda atenuada (dB)
gpass = 1
gstop = 40

# conversión de la ganancia en dB a escala lineal
delta_p = 1 - 10 ** (-gpass / 20)
delta_s = 10 ** (-gstop / 20)

print('delta_p = {:.6f}'.format(delta_p))
print('delta_s = {:.6f}'.format(delta_s))

### Filtros IIR

nw = 1024
### Filtro Butterworth
Nb, Wn = signal.buttord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
zb, pb, kb = signal.butter(Nb, Wn, btype='low', analog=False, output='zpk', fs=2*np.pi)
w, Hb = signal.freqz_zpk(zb, pb, kb, worN=nw, whole=False, fs=2*np.pi)

### Filtro Chebyshev tipo I
Nc1, Wn = signal.cheb1ord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
zc1, pc1, kc1 = signal.cheby1(Nc1, gpass, Wn, btype='low', analog=False, output='zpk', fs=2*np.pi)
_, Hc1 = signal.freqz_zpk(zc1, pc1, kc1, worN=nw, whole=False, fs=2*np.pi)

### Filtro Chebyshev tipo II
Nc2, Wn = signal.cheb2ord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
zc2, pc2, kc2 = signal.cheby2(Nc2, gstop, Wn, btype='low', analog=False, output='zpk', fs=2*np.pi)
_, Hc2 = signal.freqz_zpk(zc2, pc2, kc2, worN=nw, whole=False, fs=2*np.pi)

### Filtro elíptico
Ne, Wn = signal.ellipord(wp, ws, gpass, gstop, analog=False, fs=2*np.pi)
ze, pe, ke = signal.ellip(Ne, gpass, gstop, Wn, btype='low', analog=False, output='zpk', fs=2*np.pi)
_, He = signal.freqz_zpk(ze, pe, ke, worN=nw, whole=False, fs=2*np.pi)

### Filtros FIR

# conversión de parámetros
delta_pp = delta_p / (2 - delta_p)
delta_sp = 2 * delta_s / (2 - delta_p)

### Ventana de Kaiser
delta = np.min([delta_pp, delta_sp])
A = -20 * np.log10(delta)
# frecuencia de corte
wc = (wp + ws) / 2
# ancho de la banda de transición
Dw = ws - wp

print('Parámetros del filtro de Kaiser')
print('A = {:.6f}'.format(A))
print('Dw = {:.6f}'.format(Dw))

# cálculo de beta
if A > 50:
    beta = 0.1102 * (A - 8.8)
elif A < 21:
    beta = 0
else:
    beta = 0.5842 * (A - 21) ** 0.4 + 0.07886 * (A - 21)

# cálculo de M
Mk = math.ceil((A - 7.95) / (2.285 * Dw))

print('beta = {:.6f}'.format(beta))
print('M = {}'.format(Mk))

# Construcción de la respuesta al impulso.
wk = signal.windows.kaiser(Mk + 1, beta)
n = np.arange(Mk + 1)
alpha = Mk / 2
hk = np.sin(wc * (n - alpha)) / (np.pi * (n - alpha)) * wk
hk /= (1 + delta_pp)
if Mk % 2 == 0:
    hk[Mk // 2] = wc / np.pi

# Cálculo de la respuesta en frecuencia
_, Hk = signal.freqz(hk, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
# Cálculo de los ceros
zk = poly.polyroots(hk)

### Parks-McClellan

# predicción del valor de M
Mpm = math.ceil((-10 * np.log10(delta_pp * delta_sp) - 13) / (2.324 * Dw))
# M calculado de esta forma no llega a cumplr las especificaciones. Se cambia.
Mpm = 44
print('Parámetros del filtro de Parks-McClellan')
print('M = {}'.format(Mpm))

# Construcción de la respuesta al impulso.
wb = [0, wp, ws, np.pi]
K = delta_pp / delta_sp

hpm = signal.remez(Mpm + 1, wb, [1, 0], [1 / K, 1], fs=2 * np.pi)
hpm /= (1 + delta_pp)

# Cálculo de la respuesta en frecuencia
_, Hpm = signal.freqz(hpm, 1, worN=w, whole=False, plot=None, fs=2*np.pi)
# Cálculo de los ceros
zpm = poly.polyroots(hpm)

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
plt.plot(w, np.absolute(Hk), label=r'$\textrm{{Kaiser}}\;(M={})$'.format(Mk))
plt.plot(w, np.absolute(Hpm), label=r'$\textrm{{Parks-McClellan}}\;(M={})$'.format(Mpm))
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
plt.text(wp, xtm, '${:.1f}\pi$'.format(wp / np.pi), fontsize=fontsize, ha='center', va='baseline')
plt.plot([ws, ws], [0, xtl], 'k-', lw=1)
plt.text(ws, xtm, '${:.1f}\pi$'.format(ws / np.pi), fontsize=fontsize, ha='center', va='baseline')
plt.plot([np.pi, np.pi], [0, xtl], 'k-', lw=1)
plt.text(np.pi, xtm, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm, xtm, '$0$', fontsize=fontsize, ha='right', va='baseline')
# eje y
plt.plot([0, xtl], [1 - delta_p, 1 - delta_p], 'k-', lw=1)
plt.text(ytm, 1 - delta_p, '${:.5f}$'.format(1 - delta_p), fontsize=fontsize, ha='right', va='top')
plt.plot([0, xtl], [1, 1], 'k-', lw=1)
plt.text(ytm, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, xtl], [delta_s, delta_s], 'k-', lw=1)
plt.text(ytm, delta_s, '${:.5f}$'.format(delta_s), fontsize=fontsize, ha='right', va='center')
# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.text(0.06, ymax_ax, '$|H(e^{j\omega})|$', fontsize=fontsize, ha='left', va='center')

plt.legend(bbox_to_anchor=(0.55, 0.5), loc='center left', frameon=False, fontsize=fontsize, framealpha=1)
plt.axis('off')
# plt.savefig('filter_design_upsampling_filter_sec_07_10.pdf', bbox_inches='tight')

#
# Diagrama de polos y ceros
#
xmax = 2.5
xmin = -1.2
ymax = 1.5
ymin = -1.2
ymax_ax = 2
ytm = 0.12
xtm = -0.3

fs2 = 10

# circulo unidad
theta = np.linspace(0, 2 * np.pi, 100)
xc = np.cos(theta)
yc = np.sin(theta)

zeros_marker_style = dict(marker='o', linestyle='', markersize=6, markerfacecolor='w', markeredgecolor='k',
                          markeredgewidth=1.5)
poles_marker_style = dict(marker='x', linestyle='', markersize=6, markeredgecolor='k', markeredgewidth=2)


fig = plt.figure(1, figsize=(8, 9), frameon=False)
# Butterworth
ax = plt.subplot2grid((6, 4), (0, 0), rowspan=2, colspan=2)
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

xtit = (xmin + xmax) / 2
plt.text(xtit, 1.85, r'$\textrm{{Butterworth}}\;(N={})$'.format(Nb),
         fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

# Cheby 1
ax = plt.subplot2grid((6, 4), (0, 2), rowspan=2, colspan=2)
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

plt.text(xtit, 1.85, r'$\textrm{{Chebyshev tipo I}}\;(N={})$'.format(Nc1),
         fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

# Cheby 2
ax = plt.subplot2grid((6, 4), (2, 0), rowspan=2, colspan=2)
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

plt.text(xtit, 1.85, r'$\textrm{{Chebyshev tipo II}}\;(N={})$'.format(Nc2),
         fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

# Elíptico
ax = plt.subplot2grid((6, 4), (2, 2), rowspan=2, colspan=2)
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

plt.text(xtit, 1.85, r'$\textrm{{El\'iptico}}\;(N={})$'.format(Ne),
         fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

# Kaiser
ax = plt.subplot2grid((6, 4), (4, 0), rowspan=2, colspan=2)
plt.ylim(ymin, ymax_ax)
plt.xlim(xmin, xmax)
plt.gca().set_aspect('equal', adjustable='box')
# ejes
# axis arrows
plt.plot([xmin, 1.5], [0, 0], 'k-', lw=1.5)
plt.annotate("", xytext=(1.6, 0), xycoords='data', xy=(xmax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.plot([1.5, 1.5], [-0.1, 0.1], 'k-', lw=1)
plt.plot([1.6, 1.6], [-0.1, 0.1], 'k-', lw=1)
plt.annotate("", xytext=(0, ymin), xycoords='data', xy=(0, ymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
# circulo unidad
ax.plot(xc, yc, 'k-', lw=1)
# polos y ceros
plt.plot(zk.real, zk.imag, **zeros_marker_style, zorder=10)
# cero grande
zx = 1.9
plt.plot(zx, 0, **zeros_marker_style, zorder=10)
plt.text(zx, xtm, r'${:.2f}$'.format(np.max(zk.real)), fontsize=9, ha='center', va='baseline')

# etiquetas de los ejes
plt.text(xmax, xtm, r'$\textrm{Re}(z)$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm, ymax, r'$\textrm{Im}(z)$', fontsize=fontsize, ha='left', va='center')

plt.text(xtit, 1.85, r'$\textrm{{Kaiser}}\;(M={})$'.format(Mk),
         fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

# Parks-McClellam
ax = plt.subplot2grid((6, 4), (4, 2), rowspan=2, colspan=2)
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
plt.plot(zpm.real, zpm.imag, **zeros_marker_style, zorder=10)

# etiquetas de los ejes
plt.text(xmax, xtm, r'$\textrm{Re}(z)$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm, ymax, r'$\textrm{Im}(z)$', fontsize=fontsize, ha='left', va='center')

plt.text(xtit, 1.85, r'$\textrm{{Parks-McClellan}}\;(M={})$'.format(Mpm),
         fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

plt.savefig('filter_design_upsampling_filter_sec_07_10_zero_pole.pdf', bbox_inches='tight')

###############################################

# x y ticks labels margin
display_length = 6

xmin = 0
dw = 0.15
xmax = wp + dw
xmax_ax = xmax
xmin_ax = xmin
ymin_ax = (1 - 1.5 * delta_p)
ymax_ax = 1 + 0.5 * delta_p

xticks = np.arange(0, xmax_ax / np.pi, 0.1)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[0] = '$0$'
xticks *= np.pi

fig = plt.figure(2, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hb), label=r'$\textrm{{Butterworth}}\;(N={})$'.format(Nb))
plt.plot(w, np.absolute(Hc1), label=r'$\textrm{{Chebyshev tipo I}}\;(N={})$'.format(Nc1))
plt.plot(w, np.absolute(Hc2), label=r'$\textrm{{Chebyshev tipo II}}\;(N={})$'.format(Nc2))
plt.plot(w, np.absolute(He), label=r'$\textrm{{El\'iptico}}\;(N={})$'.format(Ne))
plt.plot(w, np.absolute(Hk), label=r'$\textrm{{Kaiser}}\;(M={})$'.format(Mk))
plt.plot(w, np.absolute(Hpm), label=r'$\textrm{{Parks-McClellan}}\;(M={})$'.format(Mpm))
# mascara
plt.plot([0, ws], [1, 1], 'k-')
plt.plot([0, wp], [1 - delta_p, 1 - delta_p], 'k-')
plt.plot([ws, np.pi], [delta_s, delta_s], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p], 'k-')
plt.plot([ws, ws], [delta_s, 1 + delta_p], 'k-')
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
plt.xticks(xticks, xticks_labels, usetex=True)
# etiquetas de los ejes
plt.xlabel('$\omega$', fontsize=fontsize)
#plt.ylabel('$|H(e^{j\omega})|$', fontsize=fontsize)

xmin = ws - dw
xmax = np.pi
xmax_ax = xmax
xmin_ax = xmin
ymin_ax = 0
ymax_ax = 2.5 * delta_s

xticks = np.arange(round(ws / np.pi, 1), 1.01, 0.1)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[-1] = '$\pi$'
xticks *= np.pi

ax = plt.subplot2grid((4, 4), (0, 2), rowspan=4, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hb), label=r'$\textrm{{Butterworth}}\;(N={})$'.format(Nb))
plt.plot(w, np.absolute(Hc1), label=r'$\textrm{{Chebyshev tipo I}}\;(N={})$'.format(Nc1))
plt.plot(w, np.absolute(Hc2), label=r'$\textrm{{Chebyshev tipo II}}\;(N={})$'.format(Nc2))
plt.plot(w, np.absolute(He), label=r'$\textrm{{El\'iptico}}\;(N={})$'.format(Ne))
plt.plot(w, np.absolute(Hk), label=r'$\textrm{{Kaiser}}\;(M={})$'.format(Mk))
plt.plot(w, np.absolute(Hpm), label=r'$\textrm{{Parks-McClellan}}\;(M={})$'.format(Mpm))
# mascara
plt.plot([0, ws], [1, 1], 'k-')
plt.plot([0, wp], [1 - delta_p, 1 - delta_p], 'k-')
plt.plot([ws, np.pi], [delta_s, delta_s], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p], 'k-')
plt.plot([ws, ws], [delta_s, 1 + delta_p], 'k-')
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
plt.xticks(xticks, xticks_labels, usetex=True)
#plt.yticks(yticks, yticks_labels, usetex=True)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
# etiquetas de los ejes
plt.xlabel('$\omega$', fontsize=fontsize)
plt.legend(bbox_to_anchor=(0.1, 0.45), loc='lower left', frameon=False, fontsize=fontsize, framealpha=1)

plt.savefig('filter_design_upsampling_filter_sec_07_10_zoom.pdf', bbox_inches='tight')


plt.show()





