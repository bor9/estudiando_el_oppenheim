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
Nx = 26
# largo de x[n]
N = 5

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
x = (N - np.arange(N)) / N

# secuencia periódica
xt = np.convolve(delta1, x, mode='same')
xt_rev = np.flip(xt)
xt_rev = np.roll(xt_rev, 1)

n = n + N // 2

x_pad = np.zeros(Naux)
x_pad[nmax - N // 2: nmax - N // 2 + N] = x
x_rev_pad = np.zeros(Naux)
x_rev = xt_rev[nmax - N // 2: nmax - N // 2 + N]
x_rev_pad[nmax - N // 2: nmax - N // 2 + N] = x_rev

n = n[Npad:-Npad]
x_pad = x_pad[Npad:-Npad]
xt = xt[Npad:-Npad]
xt_rev = xt_rev[Npad:-Npad]
x_rev_pad = x_rev_pad[Npad:-Npad]


fontsize = 12
fontsize2 = 9
xtm = -0.26
xtit = 1.3

ymin = -0.2
ymax = 1.4
nmin_ax = n[0] - 1
nmax_ax = n[-1] + 1


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
    if x[i].is_integer():
        plt.text(i, x[i] + 0.08, '{:.0f}'.format(x[i]), fontsize=fontsize2, ha='center', va='bottom')
    else:
        plt.text(i, x[i] + 0.08, '{:.1f}'.format(x[i]), fontsize=fontsize2, ha='center', va='bottom')

plt.axis('off')

ax = plt.subplot2grid((4, 4), (1, 0), rowspan=1, colspan=4)
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
    if x[i].is_integer():
        plt.text(i, x[i] + 0.08, '{:.0f}'.format(x[i]), fontsize=fontsize2, ha='center', va='bottom')
    else:
        plt.text(i, x[i] + 0.08, '{:.1f}'.format(x[i]), fontsize=fontsize2, ha='center', va='bottom')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 0), rowspan=1, colspan=4)
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
plt.text(nmax_ax, xtm, '$k$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$\tilde{X}_1[k]=N\tilde{x}[-k]$', fontsize=fontsize, ha='left', va='baseline')
for i in np.arange(N):
    if (N * x_rev[i]).is_integer():
        plt.text(i, x_rev[i] + 0.08, '{:.0f}'.format(N * x_rev[i]), fontsize=fontsize2, ha='center', va='bottom')
    else:
        plt.text(i, x_rev[i] + 0.08, '{:.1f}'.format(N * x_rev[i]), fontsize=fontsize2, ha='center', va='bottom')
plt.axis('off')


ax3 = plt.subplot2grid((4, 4), (3, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x_rev_pad, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$k$', fontsize=fontsize, ha='right', va='baseline')
latex_tit = r'$X_1[k]=' \
            r'\left\{' \
            r'\begin{array}{ll}' \
            r'\tilde{X}_1[k], & 0\leq k\leq N-1,\\' \
            r'0, & \textrm{en otro caso.}' \
            r'\end{array}' \
            r'\right.$'
plt.text(nmin_ax, xtit, latex_tit, fontsize=fontsize, ha='left', va='baseline')
for i in np.arange(N):
    if (N * x_rev[i]).is_integer():
        plt.text(i, x_rev[i] + 0.08, '{:.0f}'.format(N * x_rev[i]), fontsize=fontsize2, ha='center', va='bottom')
    else:
        plt.text(i, x_rev[i] + 0.08, '{:.1f}'.format(N * x_rev[i]), fontsize=fontsize2, ha='center', va='bottom')
plt.axis('off')

con = ConnectionPatch(xyB=(-0.5, ymin), xyA=(-0.5, ymax), coordsA="data", coordsB="data",
                      axesB=ax3, axesA=ax0, color='k', lw=1, linestyle='dashed', zorder=100)
ax0.add_artist(con)
con = ConnectionPatch(xyB=(N-0.5, ymin), xyA=(N-0.5, ymax), coordsA="data", coordsB="data",
                      axesB=ax3, axesA=ax0, color='k', lw=1, linestyle='dashed', zorder=100)
ax0.add_artist(con)

# save as pdf image
plt.savefig('dft_properties_duality_example_reversion_modulo.pdf', bbox_inches='tight')


## señal x y DFT X

### DFT de x
X = fft.fft(x)
# DTFT de x
N_dtft = 512
Xw = fft.fft(x, N_dtft)
w = 2 * np.pi * np.arange(N_dtft) / N_dtft

X_dft = fft.fft(X)
print(np.real(X))
print(np.imag(X))
print(X_dft)
Xk = 1 / (1 - np.exp(-1j*2*np.pi*np.arange(N)/N))
Xk[0] = (N+1)/2
print(np.real(Xk))
print(np.imag(Xk))

nmin = -1
nmax = N
nmin_ax = nmin - 0.5
nmax_ax = nmax + 1

ymax_ax = ymax + 1

x_pad = np.zeros(N + 1 - nmin)
x_pad[-nmin: -nmin + N] = x
n = np.arange(nmin, N + 1)

xtm = -0.3

fig = plt.figure(1, figsize=(9, 4.5), frameon=False)
ax = plt.subplot2grid((4, 8), (1, 0), rowspan=2, colspan=3)
plt.xlim(nmin_ax - 0.5, nmax_ax + 0.5)
plt.ylim(ymin, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x_pad, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
xtit = 1.5
plt.text(nmin_ax, xtit, r'$x[n],\;N={}$'.format(N), fontsize=fontsize, ha='left', va='baseline')
for i in np.arange(N):
    if x[i].is_integer():
        plt.text(i, x[i] + 0.08, '{:.0f}'.format(x[i]), fontsize=fontsize2, ha='center', va='bottom')
    else:
        plt.text(i, x[i] + 0.08, '{:.1f}'.format(x[i]), fontsize=fontsize2, ha='center', va='bottom')
for i in np.arange(N + 1):
    plt.text(i, xtm, '${}$'.format(i), fontsize=fontsize, ha='center', va='baseline')
plt.axis('off')

wmin = 0.1
wmax = 2 * np.pi
dw = 0.4
wmin_ax = wmin - dw
wmax_ax = wmax + dw

Re_min = -0.4
Re_max = 3.6
delta_y = 0.5

display_length = 6

ReX = np.zeros(N + 1)
ReX[0: N] = np.real(X)
nw = 2 * np.pi * np.arange(N + 1) / N

xtm1 = -1
xtm2 = -1.6
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = cycle(prop_cycle.by_key()['color'])
c = next(colors)

ax = plt.subplot2grid((4, 8), (0, 3), rowspan=2, colspan=5)
plt.xlim(wmin_ax, wmax_ax)
plt.ylim(Re_min - delta_y, Re_max)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)

plt.annotate("", xytext=(wmin_ax, 0), xycoords='data', xy=(wmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(nw, ReX, linefmt='k', markerfmt='sk', use_line_collection=True,
                                    label=r'$\operatorname{Re}\{X[k]\}$')
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
plt.plot(w, np.real(Xw), color=c, lw=1, label=r'$\operatorname{Re}\{X(e^{j\omega})\}$', zorder=1)

for i in np.arange(N + 1):
    plt.text(nw[i], xtm1, '${}$'.format(i), fontsize=fontsize, ha='center', va='baseline')

plt.text(wmax_ax, xtm2, '$\omega$', color=c, fontsize=fontsize, ha='center', va='baseline')
plt.text(nw[0], xtm2, '$0$', color=c, fontsize=fontsize, ha='center', va='baseline')
plt.text(nw[N], xtm2, '$2\pi$', color=c, fontsize=fontsize, ha='center', va='baseline')

plt.text(wmax_ax, xtm1, '$k$', fontsize=fontsize, ha='center', va='baseline')
plt.legend(bbox_to_anchor=(0.12, 1), loc='upper left', frameon=False, fontsize=fontsize, framealpha=1)
Xre = np.real(X)
for i in np.arange(N):
    if Xre[i].is_integer():
        plt.text(nw[i], Xre[i] + 0.12, '{:.0f}'.format(Xre[i]), fontsize=fontsize2, ha='center', va='bottom')
    else:
        plt.text(nw[i], Xre[i] + 0.12, '{:.1f}'.format(Xre[i]), fontsize=fontsize2, ha='center', va='bottom')
plt.axis('off')


Im_min = -2
Im_max = 2

ImX = np.zeros(N + 1)
ImX[0: N] = np.imag(X)

ax = plt.subplot2grid((4, 8), (2, 3), rowspan=2, colspan=5)

plt.xlim(wmin_ax, wmax_ax)
plt.ylim(Im_min, Im_max + delta_y)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)

plt.annotate("", xytext=(wmin_ax, 0), xycoords='data', xy=(wmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(nw, ImX, linefmt='k', markerfmt='sk', use_line_collection=True,
                                    label=r'$\operatorname{Im}\{X[k]\}$')
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
plt.plot(w, np.imag(Xw), color=c, lw=1,label=r'$\operatorname{Im}\{X(e^{j\omega})\}$', zorder=1)
plt.legend(bbox_to_anchor=(0.12, 1), loc='upper left', frameon=False, fontsize=fontsize, framealpha=1)
Xim = np.imag(X)
for i in np.arange(N):
    if Xim[i] < 0:
        desp = -0.12
        posv = 'top'
        posh = 'left'
    else:
        desp = 0.12
        posv = 'bottom'
        posh = 'right'
    if Xim[i].is_integer():
        plt.text(nw[i], Xim[i] + desp, '{:d}'.format(int(Xim[i])), fontsize=fontsize2, ha=posh, va=posv)
    else:
        plt.text(nw[i], Xim[i] + desp, '{:.5f}'.format(Xim[i]), fontsize=fontsize2, ha=posh, va=posv)
plt.axis('off')

plt.savefig('dft_properties_duality_example_dft_pair.pdf', bbox_inches='tight')

plt.show()
