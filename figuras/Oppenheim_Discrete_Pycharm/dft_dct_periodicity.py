import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack
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
N = 4
Npre = 1
Nmax = 16
# largo total
M = Nmax + Npre + 1

n = np.arange(-Npre, Nmax + 1)
x = N - np.arange(N)

# se arma un período de cada señal
x1 = np.zeros(2 * N - 2)
x1[0: N] = x
x1[N: 2 * N - 2] = np.flip(x[1: N - 1])
x2 = np.zeros(2 * N)
x2[0: N] = x
x2[N:] = np.flip(x)
x3 = np.zeros(4 * N)
x3[0: N] = x
x3[N + 1: 2 * N] = -np.flip(x[1:])
x3[2 * N: 3 * N] = -x
x3[3 * N + 1: 4 * N] = np.flip(x[1:])
x4 = np.zeros(4 * N)
x4[0: N] = x
x4[N: 2 * N] = -np.flip(x)
x4[2 * N: 3 * N] = -x
x4[3 * N: 4 * N] = np.flip(x)


X_pad = np.zeros((M, 4))
for i, ni in enumerate(n):
    X_pad[i, 0] = x1[ni % (2 * N - 2)]
    X_pad[i, 1] = x2[ni % (2 * N)]
    X_pad[i, 2] = x3[ni % (4 * N)]
    X_pad[i, 3] = x4[ni % (4 * N)]


# Gráfica
fontsize = 12
fontsize2 = 10
xtm = -0.8
arrow_text = -2.8
arrow_point = -0.8

nmin_ax = n[0] - 1
nmax_ax = n[-1] + 1.5
ymax = np.max(x) + 1
ymin = -np.max(x) - 0.5

ytitle = 1.4

fig = plt.figure(0, figsize=(9, 5), frameon=False)

i = 0
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, X_pad[:, i], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
(markers, stemlines, bl) = plt.stem(np.arange(N), x, linefmt='r', markerfmt='sr', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(x1.shape[0], xtm, '${}$'.format(x1.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * x1.shape[0], xtm, '${}$'.format(2 * x1.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(x1):
    plt.text(i, xi + 0.4, '${:.0f}$'.format(xi), fontsize=fontsize2, ha='center', va='baseline')

plt.annotate(r'$\textrm{WS}$', xytext=(0, arrow_text), xycoords='data', xy=(0, arrow_point), textcoords='data',
             fontsize=fontsize, va='baseline', ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
plt.annotate(r'$\textrm{WS}$', xytext=(N - 1, arrow_text), xycoords='data', xy=(N - 1, arrow_point), textcoords='data',
             fontsize=fontsize, va='baseline', ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
plt.text(nmin_ax, ymax, r'$\tilde{x}_1[n]$', fontsize=fontsize, ha='center', va='top')

plt.axis('off')

i = 1
ax = plt.subplot2grid((4, 4), (0, 2), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, X_pad[:, i], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
(markers, stemlines, bl) = plt.stem(np.arange(N), x, linefmt='r', markerfmt='sr', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(x2.shape[0], xtm, '${}$'.format(x2.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * x2.shape[0], xtm, '${}$'.format(2 * x2.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(x2):
    plt.text(i, xi + 0.4, '${:.0f}$'.format(xi), fontsize=fontsize2, ha='center', va='baseline')

plt.annotate(r'$\textrm{HS}$', xytext=(-0.5, arrow_text), xycoords='data', xy=(-0.5, arrow_point), textcoords='data',
             fontsize=fontsize, va='baseline', ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
plt.annotate(r'$\textrm{HS}$', xytext=(N - 0.5, arrow_text), xycoords='data', xy=(N - 0.5, arrow_point),
             textcoords='data', fontsize=fontsize, va='baseline', ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
plt.text(nmin_ax, ymax, r'$\tilde{x}_2[n]$', fontsize=fontsize, ha='right', va='top')

plt.axis('off')

i = 2
ax = plt.subplot2grid((4, 4), (2, 0), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, X_pad[:, i], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
(markers, stemlines, bl) = plt.stem(np.arange(N), x, linefmt='r', markerfmt='sr', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(x3.shape[0], xtm, '${}$'.format(x3.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(x3[:2 * N + 1]):
    if xi >= 0:
        bl = xi + 0.4
        va = 'baseline'
    else:
        bl = xi - 0.4
        va = 'top'
    plt.text(i, bl, '${:.0f}$'.format(xi), fontsize=fontsize2, ha='center', va=va)

plt.annotate(r'$\textrm{WS}$', xytext=(0, arrow_text), xycoords='data', xy=(0, arrow_point), textcoords='data',
             fontsize=fontsize, va='baseline', ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
plt.annotate(r'$\textrm{WA}$', xytext=(N, arrow_text), xycoords='data', xy=(N, arrow_point), textcoords='data',
             fontsize=fontsize, va='baseline', ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
plt.text(nmin_ax, ymax, r'$\tilde{x}_3[n]$', fontsize=fontsize, ha='center', va='top')

plt.axis('off')

i = 3
ax = plt.subplot2grid((4, 4), (2, 2), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, X_pad[:, i], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
(markers, stemlines, bl) = plt.stem(np.arange(N), x, linefmt='r', markerfmt='sr', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(x4.shape[0], xtm, '${}$'.format(x4.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(x4[:2 * N]):
    if xi >= 0:
        bl = xi + 0.4
        va = 'baseline'
    else:
        bl = xi - 0.4
        va = 'top'
    if xi != -1:
        plt.text(i, bl, '${:.0f}$'.format(xi), fontsize=fontsize2, ha='center', va=va)

plt.annotate(r'$\textrm{HS}$', xytext=(-0.5, arrow_text), xycoords='data', xy=(-0.5, arrow_point), textcoords='data',
             fontsize=fontsize, va='baseline', ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
plt.annotate(r'$\textrm{HA}$', xytext=(N - 0.5, arrow_text), xycoords='data', xy=(N - 0.5, arrow_point),
             textcoords='data', fontsize=fontsize, va='baseline', ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
plt.text(nmin_ax, ymax, r'$\tilde{x}_4[n]$', fontsize=fontsize, ha='right', va='top')

plt.axis('off')

# save as pdf image
plt.savefig('dft_dct_periodicity_time.pdf', bbox_inches='tight')

#
# Construcción de las extensiones tipo 1 y tipo 2
#

xe1 = np.zeros(2 * N - 2)
xe1[0: N] = x
# x_alpha
alpha = np.ones(2 * N - 2)
alpha[[0, N - 1]] = 0.5
xa = xe1 * alpha

xe2 = np.zeros(2 * N)
xe2[0: N] = x

x1_pad = np.full(shape=M, fill_value=np.nan)
x1_pad[Npre: Npre + 2 * N - 2] = xe1
xa_pad = np.full(shape=M, fill_value=np.nan)
xa_pad[Npre: Npre + 2 * N - 2] = xa

x2_pad = np.full(shape=M, fill_value=np.nan)
x2_pad[Npre: Npre + 2 * N] = xe2

# x_alpha[-n]
xan = np.flip(xa)
xan = np.roll(xan, 1)
x2n = np.flip(xe2)

xa_per = np.zeros(M)
xan_per = np.zeros(M)
x2_per = np.zeros(M)
x2n_per = np.zeros(M)

for i, ni in enumerate(n):
    xa_per[i] = xa[ni % (2 * N - 2)]
    xan_per[i] = xan[ni % (2 * N - 2)]
    x2_per[i] = xe2[ni % (2 * N)]
    x2n_per[i] = x2n[ni % (2 * N)]


ymax = np.max(x) + 1
ymin = -1


fig = plt.figure(1, figsize=(9, 7), frameon=False)

######### DCT-1 #########

ax = plt.subplot2grid((4, 4), (0, 0), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x1_pad, linefmt='k', markerfmt='sk', use_line_collection=True, label='$x[n]$')
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
(markers, stemlines, bl) = plt.stem(n, xa_pad, linefmt='r', markerfmt='sr', use_line_collection=True,
                                    label=r'$x_\alpha[n]$')
plt.setp(markers, markersize=3.5, zorder=10)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(x1.shape[0], xtm, '${}$'.format(x1.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(xe1):
    plt.text(i, xi + 0.4, '${:.0f}$'.format(xi), fontsize=fontsize2, ha='center', va='baseline')
plt.legend(frameon=False, fontsize=fontsize, framealpha=1)
plt.axis('off')

ax = plt.subplot2grid((4, 4), (1, 0), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, xa_per, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(x1.shape[0], xtm, '${}$'.format(x1.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * x1.shape[0], xtm, '${}$'.format(2 * x1.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(xa):
    if xi.is_integer():
        plt.text(i, xi + 0.4, '${:.0f}$'.format(xi), fontsize=fontsize2, ha='center', va='baseline')
    else:
        plt.text(i, xi + 0.4, '${:.1f}$'.format(xi), fontsize=fontsize2, ha='center', va='baseline')

plt.text(nmax_ax, ymax, r'$x_\alpha[((n))_{2N-2}]$', fontsize=fontsize, ha='right', va='baseline')

plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 0), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, xan_per, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(x1.shape[0], xtm, '${}$'.format(x1.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * x1.shape[0], xtm, '${}$'.format(2 * x1.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(xan):
    if xi.is_integer():
        plt.text(i, xi + 0.4, '${:.0f}$'.format(xi), fontsize=fontsize2, ha='center', va='baseline')
    else:
        plt.text(i, xi + 0.4, '${:.1f}$'.format(xi), fontsize=fontsize2, ha='center', va='baseline')

plt.text(nmax_ax, ymax, r'$x_\alpha[((-n))_{2N-2}]$', fontsize=fontsize, ha='right', va='baseline')

plt.axis('off')

ax = plt.subplot2grid((4, 4), (3, 0), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, xa_per + xan_per, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(x1.shape[0], xtm, '${}$'.format(x1.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * x1.shape[0], xtm, '${}$'.format(2 * x1.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(x1):
    plt.text(i, xi + 0.4, '${:.0f}$'.format(xi), fontsize=fontsize2, ha='center', va='baseline')

plt.text(nmax_ax, ymax, r'$\tilde{x}_1[n]$', fontsize=fontsize, ha='right', va='baseline')

plt.axis('off')

######### DCT-2 #########

ax = plt.subplot2grid((4, 4), (0, 2), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x2_pad, linefmt='k', markerfmt='sk', use_line_collection=True, label='$x[n]$')
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xe2.shape[0], xtm, '${}$'.format(xe2.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(xe2):
    plt.text(i, xi + 0.4, '${:.0f}$'.format(xi), fontsize=fontsize2, ha='center', va='baseline')
plt.text(nmax_ax, 4, r'$x[n]$', fontsize=fontsize, ha='right', va='top')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (1, 2), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x2_per, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xe2.shape[0], xtm, '${}$'.format(xe2.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * xe2.shape[0], xtm, '${}$'.format(2 * xe2.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(xe2):
    plt.text(i, xi + 0.4, '${:.0f}$'.format(xi), fontsize=fontsize2, ha='center', va='baseline')
plt.text(nmax_ax, ymax, r'$x[((n))_{2N}]$', fontsize=fontsize, ha='right', va='baseline')

plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 2), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x2n_per, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xe2.shape[0], xtm, '${}$'.format(xe2.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * xe2.shape[0], xtm, '${}$'.format(2 * xe2.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(x2n):
    plt.text(i, xi + 0.4, '${:.0f}$'.format(xi), fontsize=fontsize2, ha='center', va='baseline')

plt.text(nmax_ax, ymax, r'$x[((-n-1))_{2N}]$', fontsize=fontsize, ha='right', va='baseline')

plt.axis('off')

ax = plt.subplot2grid((4, 4), (3, 2), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x2_per + x2n_per, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xe2.shape[0], xtm, '${}$'.format(xe2.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * xe2.shape[0], xtm, '${}$'.format(2 * xe2.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(x2):
    plt.text(i, xi + 0.4, '${:.0f}$'.format(xi), fontsize=fontsize2, ha='center', va='baseline')

plt.text(nmax_ax, ymax, r'$\tilde{x}_2[n]$', fontsize=fontsize, ha='right', va='baseline')

plt.axis('off')

# save as pdf image
plt.savefig('dft_dct_periodicity_types_1_2_construction.pdf', bbox_inches='tight')

#
# Cálculo de las DFTs. Se calcula a mano para evaluarla fuera de intervalo 0 \leq k \leq N-1.
#

k = np.arange(-Npre, Nmax + 1)
# x = N - np.arange(N)
Xc1 = np.zeros(Nmax + 1 + Npre)
Xc2 = np.zeros(Nmax + 1 + Npre)
xa = xa[: N]
for i, ki in enumerate(k):
    Xc1[i] = 2 * np.sum(xa * np.cos((np.pi * ki * np.arange(N)) / (N - 1)))
    Xc2[i] = 2 * np.sum(x * np.cos((np.pi * ki * (2 * np.arange(N) + 1)) / (2 * N)))

# DCT calculada con las funciones de scipy
Xc1s = fftpack.dct(x, type=1, norm=None)
Xc2s = fftpack.dct(x, type=2, norm=None)
print(np.sum(np.abs(Xc1[Npre: N + Npre]-Xc1s)))
print(np.sum(np.abs(Xc2[Npre: N + Npre]-Xc2s)))
print(Xc1)

xtm = -3.5
nmin_ax = k[0] - 1
nmax_ax = k[-1] + 1.5
ymax = np.max(np.abs(np.concatenate((Xc1, Xc2)))) * 1.2
ymin = -ymax
eps = 1e-12

fig = plt.figure(2, figsize=(9, 3), frameon=False)

ax = plt.subplot2grid((1, 4), (0, 0), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data', zorder=1,
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(k, Xc1, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
(markers, stemlines, bl) = plt.stem(np.arange(N), Xc1[Npre: N + Npre], linefmt='r', markerfmt='sr',
                                    use_line_collection=True, zorder=100)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(x1.shape[0], xtm, '${}$'.format(x1.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * x1.shape[0], xtm, '${}$'.format(2 * x1.shape[0]), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$k$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(Xc1):
    if np.abs(Xc1[i] - round(Xc1[i])) < eps:
        if round(Xc1[i]) >= 0:
            plt.text(k[i], xi + 1.5, '${:d}$'.format(round(xi)), fontsize=fontsize2, ha='center', va='baseline')
        else:
            plt.text(k[i], xi - 1.5, '${:d}$'.format(round(xi)), fontsize=fontsize2, ha='center', va='top')

plt.text(nmin_ax, ymax, r'$X^{c1}[k]$', fontsize=fontsize, ha='center', va='top')

plt.axis('off')

ax = plt.subplot2grid((1, 4), (0, 2), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data', zorder=1,
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(k, Xc2, linefmt='k', markerfmt='sk', use_line_collection=True, zorder=100)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
(markers, stemlines, bl) = plt.stem(np.arange(N), Xc2[Npre: N + Npre], linefmt='r', markerfmt='sr',
                                    use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)

# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(2 * N, -xtm, '${}$'.format(2 * N), fontsize=fontsize, ha='center', va='top')
plt.text(4 * N, xtm, '${}$'.format(4 * N), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$k$', fontsize=fontsize, ha='right', va='baseline')
for i, xi in enumerate(Xc2):
    if np.abs(Xc2[i] - round(Xc2[i])) < eps:
        if round(Xc2[i]) != 0:
            if round(Xc2[i]) >= 0:
                plt.text(k[i], xi + 1.5, '${:d}$'.format(round(xi)), fontsize=fontsize2, ha='center', va='baseline')
            else:
                plt.text(k[i], xi - 1.5, '${:d}$'.format(round(xi)), fontsize=fontsize2, ha='center', va='top')
plt.text(nmin_ax, ymax, r'$X^{c2}[k]$', fontsize=fontsize, ha='right', va='top')

plt.axis('off')

# save as pdf image
plt.savefig('dft_dct_periodicity_types_1_2_dct.pdf', bbox_inches='tight')


plt.show()
