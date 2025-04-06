import matplotlib.pyplot as plt
import numpy as np
from scipy import fft

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}', r'\usepackage{amssymb}']


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

# largo de x[n]
L = 6
N = 2 * L
# relleno de zeros antes y despues
Npad = 3

n = np.arange(-Npad, N + Npad)

# pulsos x1[n] y x2[n] de largo L
x1 = np.zeros(N)
x1[0: L] = 1
x2 = np.zeros(N)
x2[0: L] = 1

# convolucion circular empleando la dft
x3 = fft.ifft(fft.fft(x1) * fft.fft(x2))
# se elimina la parte imaginaria
x3 = np.real(x3)
# señales intermedias
x2_modN0 = np.flip(x2)
x2_modN0 = np.roll(x2_modN0, 1)
x2_modN2 = np.roll(x2_modN0, 2)


# relleno de ceros para gráficas
x1_pad = np.zeros(N + 2 * Npad)
x1_pad[Npad : Npad + N] = x1
x2_pad = np.zeros(N + 2 * Npad)
x2_pad[Npad : Npad + N] = x2
x3_pad = np.zeros(N + 2 * Npad)
x3_pad[Npad : Npad + N] = x3
x2_modN0_pad = np.zeros(N + 2 * Npad)
x2_modN0_pad[Npad : Npad + N] = x2_modN0
x2_modN1_pad = np.zeros(N + 2 * Npad)
x2_modN1_pad[Npad : Npad + N] = x2_modN2


# Gráfica
fontsize = 12
fontsize2 = 9
xtm = -0.35
xtit = 1.4

ymin = -0.2
ymax = 1.6
nmin_ax = n[0] - 1
nmax_ax = n[-1] + 1


fig = plt.figure(0, figsize=(5, 7), frameon=False)

ax0 = plt.subplot2grid((5, 4), (0, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x2_pad, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-0.3, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x_1[n]$', fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

ax = plt.subplot2grid((5, 4), (1, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x1_pad, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-0.3, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x_2[n]$', fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

ax = plt.subplot2grid((5, 4), (2, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x2_modN0_pad, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-0.3, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x_2[((0-n))_N],\;0\leq n\leq N-1$', fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

plt.subplot2grid((5, 4), (3, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x2_modN1_pad, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-0.3, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x_2[((2-n))_N],\;0\leq n\leq N-1$', fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')


plt.subplot2grid((5, 4), (4, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x3_pad / L, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L-1, 1.06, '$L$', fontsize=fontsize, ha='center', va='bottom')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x_3[n]=x_1[n]\circledast x_2[n]$', fontsize=fontsize, ha='left', va='baseline')

plt.axis('off')
# save as pdf image
plt.savefig('dft_properties_circular_convolution_two_impulses.pdf', bbox_inches='tight')


plt.show()
