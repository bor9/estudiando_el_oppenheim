import matplotlib.pyplot as plt
import numpy as np

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


# largo de las secuencias
L = 14
P = 5
# relleno de zeros antes y despues
Npad = 9

n = np.arange(-Npad, L + Npad)

# pulsos x1[n] y x2[n] de largo L
x1 = np.ones(L)
x2 = np.arange(1, P + 1) / P

# relleno de ceros para gráficas
x1_pad = np.zeros(L + 2 * Npad)
x1_pad[Npad: Npad + L] = x1
# x2 en distintas posiciones
k = -P
x2_1 = np.zeros(L + 2 * Npad)
x2_1[Npad + k: Npad + k + P] = x2
k1 = 5
x2_2 = np.zeros(L + 2 * Npad)
x2_2[Npad + k1: Npad + k1 + P] = x2
k = L
x2_3 = np.zeros(L + 2 * Npad)
x2_3[Npad + k: Npad + k + P] = x2


# Gráfica
fontsize = 12
fontsize2 = 9
xtm = -0.27
xtit = 1.1

ymin = -0.2
ymax = 1.2
nmin_ax = n[0] - 1
nmax_ax = n[-1] + 1


fig = plt.figure(0, figsize=(9, 6), frameon=False)

ax0 = plt.subplot2grid((4, 4), (0, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x1_pad, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$m$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x_1[m]$', fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (1, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x2_1, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-1, xtm, '$-1$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-P, xtm, '$-P$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$m$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x_2[-1-m]$', fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x2_2, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.annotate('$n-P+1$', xytext=(k1, 1.6 * xtm), xycoords='data', xy=(k1, 0), textcoords='data',
             fontsize=fontsize, va="baseline", ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=0))
plt.text(k1 + P - 1, xtm, '$n$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$m$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x_2[n-m]$', fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

plt.subplot2grid((4, 4), (3, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x2_3, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.annotate('$L+P-1$', xytext=(L+P-1, 1.6 * xtm), xycoords='data', xy=(L+P-1, 0), textcoords='data',
             fontsize=fontsize, va="baseline", ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=0))
plt.text(L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$m$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x_2[L+P-1-m]$', fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

# save as pdf image
plt.savefig('dft_linear_convolution_support.pdf', bbox_inches='tight')


plt.show()
