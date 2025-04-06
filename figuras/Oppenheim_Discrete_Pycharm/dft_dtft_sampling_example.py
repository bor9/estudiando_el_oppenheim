import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


# largo total de las señales
N = 45
L = 9

# se agregan muestras por los bordes. luego se eliminan.
Npad = 8
Naux = N + Npad

nmax = Naux // 2
n = np.arange(-nmax, nmax+1)

# periodo de los tren de impulsos
N1 = 12
N2 = 7
delta1 = np.zeros(Naux)
delta2 = np.zeros(Naux)

idx = np.where(n % N1 == 0)[0]
delta1[idx] = 1
idx = np.where(n % N2 == 0)[0]
delta2[idx] = 1

x = signal.windows.triang(L)

x1 = np.convolve(delta1, x, mode='same')
x2 = np.convolve(delta2, x, mode='same')

n = n + L//2

x_pad = np.zeros(Naux)
x_pad[nmax - L//2: nmax - L//2 + L] = x


n = n[Npad:-Npad]
x_pad = x_pad[Npad:-Npad]
x1 = x1[Npad:-Npad]
x2 = x2[Npad:-Npad]

fontsize = 12
xtm = -0.21
xtit = 1.4

ymin = -0.2
ymax = 1.4
nmin_ax = n[0] - 1
nmax_ax = n[-1] + 1


fig = plt.figure(0, figsize=(9, 5), frameon=False)

ax = plt.subplot2grid((3, 4), (0, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x_pad, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L - 1, xtm, '${}$'.format(L - 1), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$x[n]$', fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

ax = plt.subplot2grid((3, 4), (1, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x1, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L - 1, xtm, '${}$'.format(L - 1), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.annotate('$N={}$'.format(N1), xytext=(N1, 2 * xtm), xycoords='data', xy=(N1, 0), textcoords='data',
             fontsize=fontsize, va="baseline", ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=0))
plt.text(nmin_ax, xtit, r'$\displaystyle\tilde{{x}}[n]=\sum_{{r=-\infty}}^\infty x[n-{}r]$'.format(N1),
         fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

ax = plt.subplot2grid((3, 4), (2, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, x2, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(L - 1, xtm, '${}$'.format(L - 1), fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.annotate('$N={}$'.format(N2), xytext=(N2, 2 * xtm), xycoords='data', xy=(N2, 0), textcoords='data',
             fontsize=fontsize, va="baseline", ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 1),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=0))
plt.text(nmin_ax, xtit, r'$\displaystyle\tilde{{x}}[n]=\sum_{{r=-\infty}}^\infty x[n-{}r]$'.format(N2),
         fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

# save as pdf image
plt.savefig('dft_dtft_sampling_example.pdf', bbox_inches='tight')


plt.show()
