import matplotlib.pyplot as plt
import numpy as np
from scipy import fft
from scipy import signal
from matplotlib.patches import ConnectionPatch

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


# largo de x[n] y h[n]
L = 18
P = 8
N = 3 * L
# relleno de zeros antes y despues
Npre = 3
Npos = P - 1
# largo total de las secuencias
M = N + Npre + Npos

# cantidad de bloques
Nb = N // L
n = np.arange(-Npre, N + Npos)

# generación de las señales
# respuesta al impulso: h[n] = e^{-kn}
h = np.exp(-0.2 * np.arange(P))
# x[n]: ruido blanco filtrado
# respuesta al impulso np.random.seed(31)para filtrar el ruido blanco y ruido blanco gaussiano
hx = signal.triang(7)
#np.random.seed(31)
#np.random.seed(34)
#np.random.seed(40)
np.random.seed(63)
u = np.random.randn(N + Npos)
x = signal.lfilter(hx, 1, u)
x /= np.max(np.abs(x))

# relleno de ceros para gráficas
h_pad = np.zeros(M)
h_pad[0: Npre] = 'nan'
h_pad[Npre: Npre + P] = h
x_pad = np.zeros(M)
x_pad[Npre: Npre + N + Npos] = x

#
# Overlap-add
#

# división en bloques
xr_oa = np.zeros((L + P - 1, Nb))
yr_oa = np.zeros((L + P - 1, Nb))
for i in np.arange(Nb):
    xi = x[i * L: (i + 1) * L]
    xr_oa[0: L, i] = xi
    yr_oa[:, i] = np.real(fft.ifft(fft.fft(xi, L + P - 1) * fft.fft(h, L + P - 1)))


# Gráfica
fontsize = 12
fontsize2 = 9
xtm = -0.3
xtm_arr = 2 * xtm

ymin = -1.2
ymax = 1.2
nmin_ax = n[0] - 3
nmax_ax = n[-1] + 3

fig = plt.figure(0, figsize=(9, 4), frameon=False)
ax = plt.subplot2grid((2, 4), (0, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
ymin2 = -0.2
ymax2 = 2 * ymax + ymin2
plt.ylim(ymin2, ymax2)
ytitle = ymax
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, 1.3 * h_pad, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, ytitle, '$h[n]$', fontsize=fontsize, ha='left', va='baseline')
plt.annotate('$P-1$', xytext=(P-1, xtm_arr), xycoords='data', xy=(P-1, 0), textcoords='data',
             fontsize=fontsize, va='baseline', ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
plt.axis('off')

ytitle = ymax
ax = plt.subplot2grid((2, 4), (1, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(-ymax, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x_pad, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
if x[0] > 0:
    plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
else:
    plt.text(0, -0.4 * xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, ytitle, r'$x[n]$', fontsize=fontsize, ha='left', va='baseline')
for i in np.arange(1, Nb + 1):
    label = '${}L$'.format(i) if i > 1 else '$L$'
    if x[i * L] > 0:
        va = 'baseline'
        xtm_arr = 2 * xtm
        rp = 1
    else:
        va = 'top'
        xtm_arr = -2.2 * xtm
        rp = 0
    plt.annotate(label, xytext=(i * L, xtm_arr), xycoords='data', xy=(i * L, 0), textcoords='data',
                 fontsize=fontsize, va=va, ha="center",
                 arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                                 patchA=None, patchB=None, shrinkA=3, shrinkB=1))
plt.axis('off')
# save as pdf image
plt.savefig('dft_lti_block_convolution_h_x.pdf', bbox_inches='tight')

xtm = -0.43
xtm2 = -0.4 * xtm
xtm_arr = 2 * xtm
fig = plt.figure(1, figsize=(9, 9), frameon=False)
for i in np.arange(Nb):
    ax = plt.subplot2grid((2 * Nb, 4), (i, 0), rowspan=1, colspan=4)
    plt.xlim(nmin_ax, nmax_ax)
    plt.ylim(-ymax, ymax)
    plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
                 arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
    (markers, stemlines, bl) = plt.stem(np.arange(i * L, (i + 1) * L + P - 1), xr_oa[:, i], linefmt='k', markerfmt='sk', use_line_collection=True)
    plt.setp(markers, markersize=3.5)
    plt.setp(bl, visible=False)
    # labels
    if xr_oa[0, i] > 0:
        plt.text(i * L, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(i * L, xtm2, '$0$', fontsize=fontsize, ha='center', va='baseline')
    if xr_oa[L - 1, i] > 0:
        plt.text((i + 1) * L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text((i + 1) * L, xtm2, '$L$', fontsize=fontsize, ha='center', va='baseline')
    plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
    plt.text(nmin_ax, ytitle, r'$x_{}[n]$'.format(i), fontsize=fontsize, ha='left', va='baseline')
    plt.annotate('$L+P-2$', xytext=((i + 1) * L + P - 2, xtm_arr), xycoords='data', xy=((i + 1) * L + P - 2, 0),
                 textcoords='data', fontsize=fontsize, va='baseline', ha="center",
                 arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                                 patchA=None, patchB=None, shrinkA=3, shrinkB=1))
    plt.axis('off')

fh = np.amax(np.abs(yr_oa))
for i in np.arange(Nb):
    ax = plt.subplot2grid((2 * Nb, 4), (i + Nb, 0), rowspan=1, colspan=4)
    plt.xlim(nmin_ax, nmax_ax)
    plt.ylim(-ymax, ymax)
    plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
                 arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
    (markers, stemlines, bl) = plt.stem(np.arange(i * L, (i + 1) * L + P - 1), yr_oa[:, i] / fh,
                                        linefmt='k', markerfmt='sk', use_line_collection=True)
    plt.setp(markers, markersize=3.5)
    plt.setp(bl, visible=False)
    # labels
    if yr_oa[0, i] > 0:
        plt.text(i * L, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(i * L, xtm2, '$0$', fontsize=fontsize, ha='center', va='baseline')
    if yr_oa[L, i] > 0:
        plt.text((i + 1) * L, xtm, '$L$', fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text((i + 1) * L, xtm2, '$L$', fontsize=fontsize, ha='center', va='baseline')
    plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
    plt.text(nmin_ax, ytitle, r'$y_{}[n]$'.format(i), fontsize=fontsize, ha='left', va='baseline')
    if yr_oa[L + P - 2, i] > 0:
        va = 'baseline'
        xtm_arr = 2 * xtm
        rp = 1
    else:
        va = 'top'
        xtm_arr = -2.2 * xtm
        rp = 0
    plt.annotate('$L+P-2$', xytext=((i + 1) * L + P - 2, xtm_arr), xycoords='data', xy=((i + 1) * L + P - 2, 0),
                 textcoords='data', fontsize=fontsize, va=va, ha="center",
                 arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                                 patchA=None, patchB=None, shrinkA=3, shrinkB=1))
    plt.axis('off')


plt.savefig('dft_lti_block_convolution_overlap_add.pdf', bbox_inches='tight')

#
# Overlap-save
#

# división en bloques
xr_os = np.zeros((L, Nb))
yr_os = np.zeros((L, Nb))
# relleno de P-1 zeros al comienzo de x[n]
x_tmp = np.zeros(N + Npos + P - 1)
x_tmp[P - 1:] = x

for i in np.arange(Nb):
    xi = x_tmp[i * (L - P + 1): i * (L - P + 1) + L]
    xr_os[0: L, i] = xi
    yr_os[:, i] = np.real(fft.ifft(fft.fft(xi, L) * fft.fft(h, L)))

fig = plt.figure(2, figsize=(9, 9), frameon=False)

i = 0
ax0 = plt.subplot2grid((2 * Nb, 4), (i, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(-ymax, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(np.arange(i * (L - P + 1), i * (L - P + 1) + L), xr_os[:, i], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
if xr_os[0, i] >= 0:
    plt.text(i * L, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
else:
    plt.text(i * L, xtm2, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, ytitle, r'$x_{}[n]$'.format(i), fontsize=fontsize, ha='left', va='baseline')
if xr_os[P - 2, i] >= 0:
    va = 'baseline'
    xtm_arr = 2 * xtm
    rp = 1
else:
    va = 'top'
    xtm_arr = -2.2 * xtm
    rp = 0
plt.annotate('$P-2$', xytext=(P - 2, xtm_arr), xycoords='data', xy=(P - 2, 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
if xr_os[L - 1, i] >= 0:
    va = 'baseline'
    xtm_arr = 2 * xtm
    rp = 1
else:
    va = 'top'
    xtm_arr = -2.2 * xtm
    rp = 0
plt.annotate('$L-1$', xytext=(L - 1 - 2, xtm_arr), xycoords='data', xy=(L - 1, 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
plt.axis('off')
if xr_os[L - (P - 1), i] >= 0:
    va = 'baseline'
    xtm_arr = 1.6 * (2 * xtm)
    rp = 1
else:
    va = 'top'
    xtm_arr = 1.6 * (-2.2 * xtm)
    rp = 0
plt.annotate('$L-(P-1)$', xytext=(L - (P - 1), xtm_arr), xycoords='data', xy=(L - (P - 1), 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))


i = 1
ax1 = plt.subplot2grid((2 * Nb, 4), (i, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(-ymax, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(np.arange(i * (L - P + 1), i * (L - P + 1) + L), xr_os[:, i], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
idx0 = i * (L - P + 1)
if xr_os[0, i] >= 0:
    plt.text(idx0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
else:
    plt.text(idx0, xtm2, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, ytitle, r'$x_{}[n]$'.format(i), fontsize=fontsize, ha='left', va='baseline')
if xr_os[P - 2, i] >= 0:
    va = 'baseline'
    xtm_arr = 2 * xtm
    rp = 1
else:
    va = 'top'
    xtm_arr = -2.2 * xtm
    rp = 0
plt.annotate('$P-2$', xytext=(idx0 + P - 2 - 2, xtm_arr), xycoords='data', xy=(idx0 + P - 2, 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
if xr_os[L - 1, i] >= 0:
    va = 'baseline'
    xtm_arr = 2 * xtm
    rp = 1
else:
    va = 'top'
    xtm_arr = -2.2 * xtm
    rp = 0
plt.annotate('$L-1$', xytext=(idx0 + L - 1 - 2, xtm_arr), xycoords='data', xy=(idx0 + L - 1, 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))

plt.axis('off')

con = ConnectionPatch(xyB=(L - 0.5, ymin), xyA=(L - 0.5, ymax), coordsA="data", coordsB="data",
                      axesB=ax1, axesA=ax0, color='k', lw=1, linestyle='dashed', zorder=100)
ax0.add_artist(con)


i = 2
ax2 = plt.subplot2grid((2 * Nb, 4), (i, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(-ymax, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(np.arange(i * (L - P + 1), i * (L - P + 1) + L), xr_os[:, i], linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
idx0 = i * (L - P + 1)
if xr_os[0, i] >= 0:
    plt.text(idx0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
else:
    plt.text(idx0, xtm2, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, ytitle, r'$x_{}[n]$'.format(i), fontsize=fontsize, ha='left', va='baseline')
if xr_os[P - 2, i] >= 0:
    va = 'baseline'
    xtm_arr = 2 * xtm
    rp = 1
else:
    va = 'top'
    xtm_arr = -2.2 * xtm
    rp = 0
plt.annotate('$P-2$', xytext=(idx0 + P - 2 - 2, xtm_arr), xycoords='data', xy=(idx0 + P - 2, 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
if xr_os[L - 1, i] >= 0:
    va = 'baseline'
    xtm_arr = 2 * xtm
    rp = 1
else:
    va = 'top'
    xtm_arr = -2.2 * xtm
    rp = 0
plt.annotate('$L-1$', xytext=(idx0 + L - 1, xtm_arr), xycoords='data', xy=(idx0 + L - 1, 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))

plt.axis('off')

con = ConnectionPatch(xyB=(idx0 + P - 1 - 0.5, ymin), xyA=(idx0 + P - 1 - 0.5, ymax), coordsA="data", coordsB="data",
                      axesB=ax2, axesA=ax1, color='k', lw=1, linestyle='dashed', zorder=100)
ax1.add_artist(con)

#########################################

i = 0
ax3 = plt.subplot2grid((2 * Nb, 4), (i + Nb, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(-ymax, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
idx0 = i * (L - P + 1)
(markers, stemlines, bl) = plt.stem(np.arange(idx0, idx0 + L), yr_os[:, i] / fh, linefmt='k',
                                    markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
(markers, stemlines, bl) = plt.stem(np.arange(idx0, idx0 + P - 1), yr_os[0: P-1, i] / fh, linefmt='r',
                                    markerfmt='sr', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
if yr_os[0, i] >= 0:
    plt.text(i * L, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
else:
    plt.text(i * L, xtm2, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, ytitle, r'$y_{{{}p}}[n]$'.format(i), fontsize=fontsize, ha='left', va='baseline')
if yr_os[P - 2, i] >= 0:
    va = 'baseline'
    xtm_arr = 2 * xtm
    rp = 1
else:
    va = 'top'
    xtm_arr = -2.2 * xtm
    rp = 0
plt.annotate('$P-2$', xytext=(P - 2, xtm_arr), xycoords='data', xy=(P - 2, 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
if yr_os[L - 1, i] >= 0:
    va = 'baseline'
    xtm_arr = 2 * xtm
    rp = 1
else:
    va = 'top'
    xtm_arr = -2.2 * xtm
    rp = 0
plt.annotate('$L-1$', xytext=(L - 1 - 2, xtm_arr), xycoords='data', xy=(L - 1, 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
plt.axis('off')
if yr_os[L - (P - 1), i] >= 0:
    va = 'baseline'
    xtm_arr = 1.6 * (2 * xtm)
    rp = 1
else:
    va = 'top'
    xtm_arr = 1.6 * (-2.2 * xtm)
    rp = 0
plt.annotate('$L-(P-1)$', xytext=(L - (P - 1), xtm_arr), xycoords='data', xy=(L - (P - 1), 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))


i = 1
ax4 = plt.subplot2grid((2 * Nb, 4), (i + Nb, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(-ymax, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
idx0 = i * (L - P + 1)
(markers, stemlines, bl) = plt.stem(np.arange(idx0, idx0 + L), yr_os[:, i] / fh, linefmt='k', markerfmt='sk',
                                    use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
(markers, stemlines, bl) = plt.stem(np.arange(idx0, idx0 + P - 1), yr_os[0: P-1, i] / fh, linefmt='r',
                                    markerfmt='sr', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
if yr_os[0, i] >= 0:
    plt.text(idx0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
else:
    plt.text(idx0, xtm2, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, ytitle, r'$y_{{{}p}}[n]$'.format(i), fontsize=fontsize, ha='left', va='baseline')
if yr_os[P - 2, i] >= 0:
    va = 'baseline'
    xtm_arr = 2 * xtm
    rp = 1
else:
    va = 'top'
    xtm_arr = -2.2 * xtm
    rp = 0
plt.annotate('$P-2$', xytext=(idx0 + P - 2 - 2, xtm_arr), xycoords='data', xy=(idx0 + P - 2, 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
if yr_os[L - 1, i] >= 0:
    va = 'baseline'
    xtm_arr = 2 * xtm
    rp = 1
else:
    va = 'top'
    xtm_arr = -2.2 * xtm
    rp = 0
plt.annotate('$L-1$', xytext=(idx0 + L - 1 - 2, xtm_arr), xycoords='data', xy=(idx0 + L - 1, 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))

plt.axis('off')

con = ConnectionPatch(xyB=(L - 0.5, ymin), xyA=(L - 0.5, ymax), coordsA="data", coordsB="data",
                      axesB=ax4, axesA=ax3, color='k', lw=1, linestyle='dashed', zorder=100)
ax3.add_artist(con)


i = 2
ax5 = plt.subplot2grid((2 * Nb, 4), (i + Nb, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(-ymax, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
idx0 = i * (L - P + 1)
(markers, stemlines, bl) = plt.stem(np.arange(idx0, idx0 + L), yr_os[:, i] / fh, linefmt='k', markerfmt='sk',
                                    use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
(markers, stemlines, bl) = plt.stem(np.arange(idx0, idx0 + P - 1), yr_os[0: P-1, i] / fh, linefmt='r',
                                    markerfmt='sr', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
if yr_os[0, i] >= 0:
    plt.text(idx0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
else:
    plt.text(idx0, xtm2, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, ytitle, r'$y_{{{}p}}[n]$'.format(i), fontsize=fontsize, ha='left', va='baseline')
if yr_os[P - 2, i] >= 0:
    va = 'baseline'
    xtm_arr = 2 * xtm
    rp = 1
else:
    va = 'top'
    xtm_arr = -2.2 * xtm
    rp = 0
plt.annotate('$P-2$', xytext=(idx0 + P - 2 - 2, xtm_arr), xycoords='data', xy=(idx0 + P - 2, 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))
if yr_os[L - 1, i] >= 0:
    va = 'baseline'
    xtm_arr = 2 * xtm
    rp = 1
else:
    va = 'top'
    xtm_arr = -2.2 * xtm
    rp = 0
plt.annotate('$L-1$', xytext=(idx0 + L - 1, xtm_arr), xycoords='data', xy=(idx0 + L - 1, 0),
             textcoords='data', fontsize=fontsize, va=va, ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, rp),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=1))

plt.axis('off')

con = ConnectionPatch(xyB=(idx0 + P - 1 - 0.5, ymin), xyA=(idx0 + P - 1 - 0.5, ymax), coordsA="data", coordsB="data",
                      axesB=ax5, axesA=ax4, color='k', lw=1, linestyle='dashed', zorder=100)
ax4.add_artist(con)


plt.savefig('dft_lti_block_convolution_overlap_save.pdf', bbox_inches='tight')


plt.show()
