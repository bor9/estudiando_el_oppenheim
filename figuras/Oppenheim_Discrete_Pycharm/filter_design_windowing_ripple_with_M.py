import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

__author__ = 'ernesto'

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


# parámetros
# largo del eneventanado (n=-(M-1)/2: (M-1)/2)
# M debe ser impar
Ms = np.array([6, 20, 40, 100])
wc = np.pi / 2

# rango de frecuencia
wmax = np.pi
wmin = -np.pi
nw = 256
w = np.linspace(0, wmax, nw)

Hs = np.zeros((2 * nw - 1, Ms.shape[0]))

for i, M in enumerate(Ms):
    n = np.arange(M / 2 + 1)
    hn = np.sin(wc * n) / (np.pi * n)
    hn[0] = wc / np.pi
    hn = np.concatenate((np.flip(hn[1: ]), hn))
    _, H = signal.freqz(hn, 1, worN=w, whole=False, plot=None, fs=2 * np.pi)
    # se elimina el componente de fase lineal para que Hz sea real
    H = np.real(H * np.exp(1j * w * M / 2))
    H = np.concatenate((np.flip(H[1: ]), H))
    Hs[:, i] = H

w = np.concatenate((np.flip(-w[1:]), w))

# valores maximos y minimos de los ejes
xmax_ax = wmax + 1
xmin_ax = -xmax_ax
ymax_ax = 1.4
ymin_ax = -0.4

# length of the ticks for all subplot (6 pixels)
display_length = 7  # in pixels
# x y ticks labels margin
xtm = -0.18
xtm2 = -0.25
ytm = -0.3
# font size
fontsize = 12

xticks = [-np.pi, -wc, 0, wc, np.pi]
xticks_labels = ['$-\pi$', '$-\omega_c$', '$0$', '$\omega_c$', '$\pi$']

fig = plt.figure(0, figsize=(8, 5), frameon=False)

ax = plt.subplot2grid((4, 4), (0, 0), rowspan=2, colspan=2)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)

# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
# Filtro ideal
plt.plot([wmin, -wc], [0, 0], 'k-', lw=1.5, zorder=1)
plt.plot([-wc, wc], [1, 1], 'k-', lw=1.5, zorder=1, label=r'$H_d(e^{j\omega})$')
plt.plot([wc, wmax], [0, 0], 'k-', lw=1.5, zorder=1)
plt.plot([-wc, -wc], [0, 1], 'k-', lw=1.5, zorder=1)
plt.plot([wc, wc], [0, 1], 'k-', lw=1.5, zorder=1)
# Filtro enventanado
i = 0
plt.plot(w, Hs[:, i], 'r-', lw=1.5, zorder=10, label=r'$H(e^{j\omega})$')
# ticks y etiquetas
for k, xtick in enumerate(xticks):
    plt.plot([xtick, xtick], [0, xtl], 'k-', lw=1)
    plt.text(xtick, xtm, xticks_labels[k], ha='center', va='baseline', fontsize=fontsize)

plt.text(xmax_ax, xtm, '$\omega$', ha='center', va='baseline', fontsize=fontsize)
plt.text(xmin_ax, ymax_ax, '$M={}$'.format(Ms[i]), ha='left', va='top', fontsize=fontsize)

plt.legend(loc='center', bbox_to_anchor=(1.1, 0.6), frameon=False, framealpha=1, fontsize=fontsize)

plt.axis('off')

#############################################################
#############################################################
ax = plt.subplot2grid((4, 4), (0, 2), rowspan=2, colspan=2)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)

# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
# Filtro ideal
plt.plot([wmin, -wc], [0, 0], 'k-', lw=1.5, zorder=1)
plt.plot([-wc, wc], [1, 1], 'k-', lw=1.5, zorder=1, label=r'$H_d(e^{j\omega})$')
plt.plot([wc, wmax], [0, 0], 'k-', lw=1.5, zorder=1)
plt.plot([-wc, -wc], [0, 1], 'k-', lw=1.5, zorder=1)
plt.plot([wc, wc], [0, 1], 'k-', lw=1.5, zorder=1)
# Filtro enventanado
i = 1
plt.plot(w, Hs[:, i], 'r-', lw=1.5, zorder=10, label=r'$H(e^{j\omega})$')
# ticks y etiquetas
for k, xtick in enumerate(xticks):
    plt.plot([xtick, xtick], [0, xtl], 'k-', lw=1)
    plt.text(xtick, xtm, xticks_labels[k], ha='center', va='baseline', fontsize=fontsize)

plt.text(xmax_ax, xtm, '$\omega$', ha='center', va='baseline', fontsize=fontsize)
plt.text(xmin_ax, ymax_ax, '$M={}$'.format(Ms[i]), ha='left', va='top', fontsize=fontsize)


plt.axis('off')

#############################################################
#############################################################
ax = plt.subplot2grid((4, 4), (2, 0), rowspan=2, colspan=2)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)

# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
# Filtro ideal
plt.plot([wmin, -wc], [0, 0], 'k-', lw=1.5, zorder=1)
plt.plot([-wc, wc], [1, 1], 'k-', lw=1.5, zorder=1, label=r'$H_d(e^{j\omega})$')
plt.plot([wc, wmax], [0, 0], 'k-', lw=1.5, zorder=1)
plt.plot([-wc, -wc], [0, 1], 'k-', lw=1.5, zorder=1)
plt.plot([wc, wc], [0, 1], 'k-', lw=1.5, zorder=1)
# Filtro enventanado
i = 2
plt.plot(w, Hs[:, i], 'r-', lw=1.5, zorder=10, label=r'$H(e^{j\omega})$')
# ticks y etiquetas
for k, xtick in enumerate(xticks):
    plt.plot([xtick, xtick], [0, xtl], 'k-', lw=1)
    plt.text(xtick, xtm, xticks_labels[k], ha='center', va='baseline', fontsize=fontsize)

plt.text(xmax_ax, xtm, '$\omega$', ha='center', va='baseline', fontsize=fontsize)
plt.text(xmin_ax, ymax_ax, '$M={}$'.format(Ms[i]), ha='left', va='top', fontsize=fontsize)


plt.axis('off')

#############################################################
#############################################################
ax = plt.subplot2grid((4, 4), (2, 2), rowspan=2, colspan=2)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)

# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
# Filtro ideal
plt.plot([wmin, -wc], [0, 0], 'k-', lw=1.5, zorder=1)
plt.plot([-wc, wc], [1, 1], 'k-', lw=1.5, zorder=1, label=r'$H_d(e^{j\omega})$')
plt.plot([wc, wmax], [0, 0], 'k-', lw=1.5, zorder=1)
plt.plot([-wc, -wc], [0, 1], 'k-', lw=1.5, zorder=1)
plt.plot([wc, wc], [0, 1], 'k-', lw=1.5, zorder=1)
# Filtro enventanado
i = 3
plt.plot(w, Hs[:, i], 'r-', lw=1.5, zorder=10, label=r'$H(e^{j\omega})$')
# ticks y etiquetas
for k, xtick in enumerate(xticks):
    plt.plot([xtick, xtick], [0, xtl], 'k-', lw=1)
    plt.text(xtick, xtm, xticks_labels[k], ha='center', va='baseline', fontsize=fontsize)

plt.text(xmax_ax, xtm, '$\omega$', ha='center', va='baseline', fontsize=fontsize)
plt.text(xmin_ax, ymax_ax, '$M={}$'.format(Ms[i]), ha='left', va='top', fontsize=fontsize)
plt.axis('off')

# save as pdf image
plt.savefig('filter_design_windowing_ripple_with_M.pdf', bbox_inches='tight')

plt.show()
