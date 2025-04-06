import matplotlib.pyplot as plt
import numpy as np
from scipy import fft
from scipy import signal

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
P = L
N = L
# relleno de zeros antes y despues
Npre = 7
Npos = 12
# largo total de las secuencias
M = N + Npre + Npos

n = np.arange(-Npre, N + Npos)

# pulsos x1[n] y x2[n] de largo L
x1 = np.ones(L)
x2 = np.ones(P)

# convolución lineal de x1[n] y x2[n]:  largo Nx3 = L + P - 1
x3 = signal.convolve(x1, x2, mode='full')
Nx3 = L + P - 1

# convolución circular de x1[n] y x2[n] de N = L puntos
x4 = fft.ifft(fft.fft(x1) * fft.fft(x2))
x4 = np.real(x4)
# convolución circular de x1[n] y x2[n] de N = 2L puntos
x5 = fft.ifft(fft.fft(x1, 2 * L) * fft.fft(x2, 2 * L))
x5 = np.real(x5)


# relleno de ceros para gráficas
x1_pad = np.zeros(M)
x1_pad[Npre: Npre + L] = x1
x2_pad = np.zeros(M)
x2_pad[Npre: Npre + P] = x2
x3_pad = np.zeros(M)
x3_pad[Npre: Npre + Nx3] = x3
# x3[n-N]
x3_1_pad = np.zeros(M)
x3_1_pad[Npre + N:  Npre + N + Nx3] = x3
# x3[n+N]
x3_2_pad = np.zeros(M)
x3_2_pad[Npre - N:  Npre - N + Nx3] = x3

x4_pad = np.zeros(M)
x4_pad[Npre: Npre + L] = x4

x5_pad = np.zeros(M)
x5_pad[Npre: Npre + 2 * L] = x5


# Gráfica
fontsize = 12
fontsize2 = 9
xtm = -0.35

ymin = -0.2
ymax = 1.6
nmin_ax = n[0] - 1
nmax_ax = n[-1] + 1

ytitle = 1.7

fig = plt.figure(0, figsize=(9, 8), frameon=False)

ax0 = plt.subplot2grid((6, 4), (0, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x1_pad, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L, xtm, '$L=P$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-0.3, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
tit = r'$' \
      r'\begin{{array}}{{l}}' \
      r'x_1[n]=x_2[n]\\' \
      r'L=P={}' \
      r'\end{{array}}' \
      r'$'.format(L)
plt.text(nmin_ax, ytitle, tit, fontsize=fontsize, ha='left', va='top')
plt.axis('off')

ax = plt.subplot2grid((6, 4), (1, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x3_pad / L, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * L - 1, xtm, '$2L-1$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L-1-0.3, 1, '$L$', fontsize=fontsize, ha='right', va='center')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, ytitle, r'$x_3[n] = x_1[n]*x_2[n]$', fontsize=fontsize, ha='left', va='top')
plt.axis('off')

ax = plt.subplot2grid((6, 4), (2, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x3_1_pad / L, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L+N-1-0.3, 1, '$L$', fontsize=fontsize, ha='right', va='center')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
tit = r'$' \
      r'\begin{{array}}{{l}}' \
      r'x_3[n-N]\\' \
      r'N=L={}' \
      r'\end{{array}}' \
      r'$'.format(L)
plt.text(nmin_ax, ytitle, tit, fontsize=fontsize, ha='left', va='top')
plt.axis('off')

plt.subplot2grid((6, 4), (3, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x3_2_pad / L, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-N, xtm, '$-N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L-N-1-0.3, 1, '$L$', fontsize=fontsize, ha='right', va='center')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
tit = r'$' \
      r'\begin{{array}}{{l}}' \
      r'x_3[n+N]\\' \
      r'N=L={}' \
      r'\end{{array}}' \
      r'$'.format(L)
plt.text(nmin_ax, ytitle, tit, fontsize=fontsize, ha='left', va='top')
plt.axis('off')


plt.subplot2grid((6, 4), (4, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x4_pad / L, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-0.3, 1, '$L$', fontsize=fontsize, ha='right', va='center')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
tit = r'$' \
      r'\begin{{array}}{{l}}' \
      r'x_1[n]\circledast x_2[n]\\' \
      r'N=L={}' \
      r'\end{{array}}' \
      r'$'.format(L)
plt.text(nmin_ax, ytitle, tit, fontsize=fontsize, ha='left', va='top')
plt.axis('off')

ax = plt.subplot2grid((6, 4), (5, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x5_pad / L, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * L, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L-1-0.3, 1, '$L$', fontsize=fontsize, ha='right', va='center')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
tit = r'$' \
      r'\begin{{array}}{{l}}' \
      r'x_1[n]\circledast x_2[n]\\' \
      r'N=2L={}' \
      r'\end{{array}}' \
      r'$'.format(2*L)
plt.text(nmin_ax, ytitle, tit, fontsize=fontsize, ha='left', va='top')
plt.axis('off')

# save as pdf image
plt.savefig('dft_circular_convolution_as_linear_with_aliasing_example_1.pdf', bbox_inches='tight')

plt.show()
