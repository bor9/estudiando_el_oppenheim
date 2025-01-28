import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from matplotlib.patches import Polygon
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


# Especificaciones del filtro pasa-altos
ws = 0.35 * np.pi
wp = 0.5 * np.pi
delta_s = 0.02
delta_p1 = 0.02
delta_p2 = 0.02

# Parámetros del filtro FIR
# ripple permitido
delta = np.min([delta_p1, delta_p2, delta_s])
A = -20 * np.log10(delta)
# frecuencia de corte
wc = (wp + ws) / 2
# ancho de la banda de transición
Dw = wp - ws

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
# se fuerza a ser un número par. en caso contrario el filtro
# tiene un cero en pi, por lo que no es un pasabajos.
# M1 es el numero entero inmediatamente superior. M2 es el numero entero par inmediatamente superior.
M = (A - 7.95) / (2.285 * Dw)
M1 = math.ceil(M)
M2 = math.ceil(M / 2) * 2

print('beta = {:.6f}'.format(beta))
print('M = {:.6f}'.format(M))
print('M (entero superior) = {}'.format(M1))
print('M (entero par superior) = {}'.format(M2))

# Diseño de filtro FIR por enventanado
# Construcción de la respuesta al impulso.
wk = signal.windows.kaiser(M2 + 1, beta)
n = np.arange(M2 + 1)
alpha = M2 / 2
hk = (np.sin(np.pi * (n - alpha)) / (np.pi * (n - alpha))-np.sin(wc * (n - alpha)) / (np.pi * (n - alpha))) * wk
if M2 % 2 == 0:
    hk[M2 // 2] = 1 - wc / np.pi

# Cálculo de la respuesta en frecuencia
nw = 2 * 1024
w = np.linspace(0, np.pi, nw)
_, Hk = signal.freqz(hk, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)

# Cálculo del error de aproximación
# Eliminación del componente de fase lineal
Ae = np.real(Hk * np.exp(1j * w * alpha))
Ae[w <= ws] = 0 - Ae[w <= ws]
Ae[w >= wp] = 1 - Ae[w >= wp]
Ae[np.logical_and(w > ws, w < wp)] = 'nan'

### Parámetros de la gráfica
fontsize = 10
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
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hk), label=r'$\textrm{{Kaiser}}\;(M={0:},\;\beta={1:.3f})$'.format(M2, beta), zorder=10)
# máscara
plt.plot([0, ws], [delta_s, delta_s], 'k-')
plt.plot([ws, np.pi], [1 + delta_p1, 1 + delta_p1], 'k-')
plt.plot([wp, np.pi], [1 - delta_p2, 1 - delta_p2], 'k-')
plt.plot([ws, ws], [delta_s, 1 + delta_p1], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p2], 'k-')
# región pintada
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [delta_s, delta_s, mask_max, mask_max]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
vert = np.vstack(([wp, xmax, xmax, wp], [0, 0, 1 - delta_p2, 1 - delta_p2]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [1 + delta_p1, 1 + delta_p1, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)

# etiquetas
# eje x
xtm2 = -0.11
plt.plot([wp, wp], [0, xtl], 'k-', lw=1)
plt.text(wp, xtm, '${:.1f}\pi$'.format(wp / np.pi), fontsize=fontsize, ha='center', va='baseline')
plt.plot([ws, ws], [0, xtl], 'k-', lw=1)
plt.text(ws, xtm, '${:.2f}\pi$'.format(ws / np.pi), fontsize=fontsize, ha='center', va='baseline')
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
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize2, ha='center', va='baseline')
plt.text(0.06, ymax_ax, '$|H(e^{j\omega})|$', fontsize=fontsize2, ha='left', va='center')

plt.legend(bbox_to_anchor=(0.5, 1), loc='center', frameon=False, fontsize=fontsize2, framealpha=1)
plt.axis('off')
#plt.savefig('filter_design_example_07_05_with_kaiser.pdf', bbox_inches='tight')

###############################################

# x y ticks labels margin
display_length = 6

xmin = 0
dw = 0.15
xmax = ws + dw
xmax_ax = xmax
xmin_ax = xmin
dd = 0.02
ymin_ax = 0
ymax_ax = delta + dd

xticks = np.arange(0, xmax_ax / np.pi, 0.1)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[0] = '$0$'
xticks *= np.pi

fig = plt.figure(1, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hk), zorder=10)
# máscara
plt.plot([0, ws], [delta_s, delta_s], 'k-')
plt.plot([ws, np.pi], [1 + delta_p1, 1 + delta_p1], 'k-')
plt.plot([wp, np.pi], [1 - delta_p2, 1 - delta_p2], 'k-')
plt.plot([ws, ws], [delta_s, 1 + delta_p1], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p2], 'k-')
# región pintada
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [delta_s, delta_s, mask_max, mask_max]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
vert = np.vstack(([wp, xmax, xmax, wp], [0, 0, 1 - delta_p2, 1 - delta_p2]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [1 + delta_p1, 1 + delta_p1, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)

# etiquetas
plt.xticks(xticks, xticks_labels, usetex=True)
# etiquetas de los ejes
plt.xlabel('$\omega$', fontsize=fontsize2)
plt.ylabel('$|H(e^{j\omega})|$', fontsize=fontsize2)

xmin = wp - dw
xmax = np.pi
xmax_ax = xmax
xmin_ax = xmin
ymin_ax = 1 - delta_p2 - dd / 2
ymax_ax = 1 + delta_p1 + dd / 2

xticks = np.arange(round(wp / np.pi, 1), 1.01, 0.1)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[-1] = '$\pi$'
xticks *= np.pi

ax = plt.subplot2grid((4, 4), (0, 2), rowspan=4, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hk), zorder=10)
# máscara
plt.plot([0, ws], [delta_s, delta_s], 'k-')
plt.plot([ws, np.pi], [1 + delta_p1, 1 + delta_p1], 'k-')
plt.plot([wp, np.pi], [1 - delta_p2, 1 - delta_p2], 'k-')
plt.plot([ws, ws], [delta_s, 1 + delta_p1], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p2], 'k-')
# región pintada
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [delta_s, delta_s, mask_max, mask_max]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
vert = np.vstack(([wp, xmax, xmax, wp], [0, 0, 1 - delta_p2, 1 - delta_p2]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [1 + delta_p1, 1 + delta_p1, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)

# etiquetas
plt.xticks(xticks, xticks_labels, usetex=True)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
# etiquetas de los ejes
plt.xlabel('$\omega$', fontsize=fontsize2)
plt.ylabel('$|H(e^{j\omega})|$', fontsize=fontsize2)

plt.savefig('filter_design_windowing_kaiser_high_pass_zoom.pdf', bbox_inches='tight')

#
# gráfica del error de aproximación
#
xmax_ax = np.pi
xmin_ax = 0
f = 1.3
ymax_ax = delta_s * f
ymin_ax = -delta_s * f

xticks = np.linspace(0, 1, 11)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[0] = '$0$'
xticks_labels[-1] = '$\pi$'
xticks *= np.pi

fig = plt.figure(2, figsize=(8, 3), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.plot(w, Ae, 'k-', lw=1.5)
plt.plot([0, np.pi], [delta, delta], 'k--', lw=1)
plt.plot([0, np.pi], [-delta, -delta], 'k--', lw=1)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel('$\omega$', fontsize=fontsize2)
plt.ylabel('$E_A(\omega)$', fontsize=fontsize2)

plt.savefig('filter_design_windowing_kaiser_high_pass_aprox_error.pdf', bbox_inches='tight')

# grafica de la ventana en el tiempo
# valores maximos y minimos de los ejes
dM = 3
xmax_ax = M2 + dM
xmin_ax = 0 - dM
ymax_ax = 0.6
ymin_ax = -0.2


fig = plt.figure(3, figsize=(5, 3), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)

plt.xlim(xmin_ax, xmax_ax)

(markers, stemlines, bl) = plt.stem(n, hk, linefmt='k', markerfmt='s', use_line_collection=True)
plt.setp(markers, markersize=3.5, markeredgecolor='k', markerfacecolor='k')
plt.setp(bl, visible=False)
plt.plot([xmin_ax, xmax_ax], [0, 0], 'k-', lw=1, zorder=-1)

plt.xlabel('$n$')

plt.savefig('filter_design_windowing_kaiser_high_pass_impulse_response.pdf', bbox_inches='tight')

#
# Respuesta en frecuencia con M impar.
#
wk = signal.windows.kaiser(M1 + 1, beta)
n = np.arange(M1 + 1)
alpha = M1 / 2
hk = (np.sin(np.pi * (n - alpha)) / (np.pi * (n - alpha))-np.sin(wc * (n - alpha)) / (np.pi * (n - alpha))) * wk
if M1 % 2 == 0:
    hk[M1 // 2] = 1 - wc / np.pi

# Cálculo de la respuesta en frecuencia
nw = 2 * 1024
w = np.linspace(0, np.pi, nw)
_, Hk = signal.freqz(hk, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)

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

fig = plt.figure(4, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hk), label=r'$\textrm{{Kaiser}}\;(M={0:},\;\beta={1:.3f})$'.format(M1, beta), zorder=10)
# máscara
plt.plot([0, ws], [delta_s, delta_s], 'k-')
plt.plot([ws, np.pi], [1 + delta_p1, 1 + delta_p1], 'k-')
plt.plot([wp, np.pi], [1 - delta_p2, 1 - delta_p2], 'k-')
plt.plot([ws, ws], [delta_s, 1 + delta_p1], 'k-')
plt.plot([wp, wp], [0, 1 - delta_p2], 'k-')
# región pintada
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [delta_s, delta_s, mask_max, mask_max]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
vert = np.vstack(([wp, xmax, xmax, wp], [0, 0, 1 - delta_p2, 1 - delta_p2]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [1 + delta_p1, 1 + delta_p1, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)

# etiquetas
# eje x
xtm2 = -0.11
plt.plot([wp, wp], [0, xtl], 'k-', lw=1)
plt.text(wp, xtm, '${:.1f}\pi$'.format(wp / np.pi), fontsize=fontsize, ha='center', va='baseline')
plt.plot([ws, ws], [0, xtl], 'k-', lw=1)
plt.text(ws, xtm, '${:.2f}\pi$'.format(ws / np.pi), fontsize=fontsize, ha='center', va='baseline')
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
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize2, ha='center', va='baseline')
plt.text(0.06, ymax_ax, '$|H(e^{j\omega})|$', fontsize=fontsize2, ha='left', va='center')

plt.legend(bbox_to_anchor=(0.5, 1), loc='center', frameon=False, fontsize=fontsize2, framealpha=1)
plt.axis('off')
plt.savefig('filter_design_windowing_kaiser_high_pass_freq_response_M_odd.pdf', bbox_inches='tight')


plt.show()

