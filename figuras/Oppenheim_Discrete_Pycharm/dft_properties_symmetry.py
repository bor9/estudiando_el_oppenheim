import matplotlib.pyplot as plt
import numpy as np
from scipy import fft
from matplotlib.patches import ConnectionPatch
from itertools import cycle

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

# largo total de las señales
Nx = 27
# largo de x[n]
N = 6

# se agregan muestras por los bordes. luego se eliminan.
Npad = N
Naux = Nx + Npad # tiene que ser un número impar

nmax = Naux // 2
n = np.arange(-nmax, nmax + 1)

# construcción de tilde x
# se construye un tren de impulsos periódico de período N
delta1 = np.zeros(Naux)
idx = np.where(n % N == 0)[0]
delta1[idx] = 1

# secuencia x[n]
x = N - np.arange(N)

# secuencia periódica
xt = np.convolve(delta1, x, mode='same')
xt_rev = np.flip(xt)
xt = np.roll(xt, -1)
xt_rev = np.roll(xt_rev, 1)

n = n + N // 2

x_pad = np.zeros(Naux)
x_pad[nmax - N // 2: nmax - N // 2 + N] = x
x_rev = xt_rev[nmax - N // 2: nmax - N // 2 + N]
x_pad_rev = np.zeros(Naux)
x_pad_rev[nmax - N // 2 - N + 1: nmax - N // 2 + 1] = np.flip(x)


n = n[Npad:-Npad]
x_pad = x_pad[Npad:-Npad]
xt = xt[Npad:-Npad]
xt_rev = xt_rev[Npad:-Npad]
x_pad_rev = x_pad_rev[Npad:-Npad]

# indice de n = 0
n0 = np.where(n == 0)[0][0]

xte = (xt + xt_rev) / 2
xto = (xt - xt_rev) / 2
xe = (x_pad + x_pad_rev) / 2
xo = (x_pad - x_pad_rev) / 2

fontsize = 12
fontsize2 = 9
xtm = -1.4

ymin = -1
ymax = N + 1
nmin_ax = n[0] - 1
nmax_ax = n[-1] + 1

xtit = ymax

fig = plt.figure(0, figsize=(9, 6), frameon=False)

ax0 = plt.subplot2grid((4, 4), (0, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x_pad, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x[n],\;N={}$'.format(N), fontsize=fontsize, ha='left', va='baseline')
for i in np.arange(N):
    plt.text(i, x_pad[n0 + i] + 0.4, '{:.0f}'.format(x_pad[n0 + i]), fontsize=fontsize2, ha='center', va='bottom')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (1, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x_pad_rev, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x[-n]$', fontsize=fontsize, ha='left', va='baseline')
for i in np.arange(N):
    plt.text(-i, x_pad_rev[n0 - i] + 0.4, '{:.0f}'.format(x_pad_rev[n0 - i]),
             fontsize=fontsize2, ha='center', va='bottom')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, xt, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$\tilde{x}[n]=x[((n))_N]$', fontsize=fontsize, ha='left', va='baseline')
for i in np.arange(N):
    plt.text(i, x[i] + 0.4, '{:.0f}'.format(x[i]), fontsize=fontsize2, ha='center', va='bottom')
plt.axis('off')


ax3 = plt.subplot2grid((4, 4), (3, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, xt_rev, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$\tilde{x}[-n]=x[((-n))_N]$', fontsize=fontsize, ha='left', va='baseline')
for i in np.arange(N):
    plt.text(i, xt_rev[n0 + i] + 0.4, '{:.0f}'.format(xt_rev[n0 + i]), fontsize=fontsize2, ha='center', va='bottom')
plt.axis('off')

con = ConnectionPatch(xyB=(-0.5, ymin), xyA=(-0.5, ymax), coordsA="data", coordsB="data",
                      axesB=ax3, axesA=ax0, color='k', lw=1, linestyle='dashed', zorder=100)
ax0.add_artist(con)
con = ConnectionPatch(xyB=(N-0.5, ymin), xyA=(N-0.5, ymax), coordsA="data", coordsB="data",
                      axesB=ax3, axesA=ax0, color='k', lw=1, linestyle='dashed', zorder=100)
ax0.add_artist(con)

# save as pdf image
plt.savefig('dft_properties_symmetry_signals.pdf', bbox_inches='tight')


######################################################
######################################################
######################################################

ymin2 = -(N + 1) / 2
ymax2 = (N + 1) / 2
xtit2 = ymax2


fig = plt.figure(1, figsize=(9, 6), frameon=False)

ax0 = plt.subplot2grid((4, 4), (0, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, xte, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$\tilde{x}_e[n]=(\tilde{x}[n]+\tilde{x}[-n])/2$',
         fontsize=fontsize, ha='left', va='baseline')
for i in np.arange(N):
    plt.text(i, xte[n0 + i] + 0.4, '{:.0f}'.format(xte[n0 + i]), fontsize=fontsize2, ha='center', va='bottom')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (1, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin2, ymax2)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, xto, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit2, r'$\tilde{x}_o[n]=(\tilde{x}[n]-\tilde{x}[-n])/2$',
         fontsize=fontsize, ha='left', va='baseline')
for i in np.arange(N):
    if xto[n0 + i] >= 0:
        plt.text(i, xto[n0 + i] + 0.4, '{:.0f}'.format(xto[n0 + i]), fontsize=fontsize2, ha='center', va='bottom')
    else:
        plt.text(i, xto[n0 + i] - 0.5, '{:.0f}'.format(xto[n0 + i]), fontsize=fontsize2, ha='center', va='top')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, xe, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x_e[n]=(x[n]+x[-n])/2$', fontsize=fontsize, ha='left', va='baseline')
for i in np.arange(N):
    if xe[n0 + i].is_integer():
        plt.text(i, xe[n0 + i] + 0.4, '{:.0f}'.format(xe[n0 + i]), fontsize=fontsize2, ha='center', va='bottom')
    else:
        plt.text(i, xe[n0 + i] + 0.4, '{:.1f}'.format(xe[n0 + i]), fontsize=fontsize2, ha='center', va='bottom')
plt.axis('off')

ax3 = plt.subplot2grid((4, 4), (3, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin2, ymax2)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, xo, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit2, r'$x_o[n]=(x[n]-x[-n])/2$', fontsize=fontsize, ha='left', va='baseline')
for i in np.arange(N):
    if xo[n0 + i].is_integer():
        plt.text(i, xo[n0 + i] + 0.4, '{:.0f}'.format(xo[n0 + i]), fontsize=fontsize2, ha='center', va='bottom')
    else:
        plt.text(i, xo[n0 + i] + 0.4, '{:.1f}'.format(xo[n0 + i]), fontsize=fontsize2, ha='center', va='bottom')
plt.axis('off')

con = ConnectionPatch(xyB=(-0.5, ymin2), xyA=(-0.5, ymax), coordsA="data", coordsB="data",
                      axesB=ax3, axesA=ax0, color='k', lw=1, linestyle='dashed', zorder=100)
ax0.add_artist(con)
con = ConnectionPatch(xyB=(N-0.5, ymin2), xyA=(N-0.5, ymax), coordsA="data", coordsB="data",
                      axesB=ax3, axesA=ax0, color='k', lw=1, linestyle='dashed', zorder=100)
ax0.add_artist(con)

# save as pdf image
plt.savefig('dft_properties_symmetry_decomposition.pdf', bbox_inches='tight')

plt.show()
