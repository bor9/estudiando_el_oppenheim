import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# largo total de las señales
Nx = 35
# largo de x[n]
N = 6
# desplazamiento temporal
m = -2

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
x = 0.8 ** (np.arange(N))

# secuencia periódica
xt1 = np.convolve(delta1, x, mode='same')
# el siguiente desplazamiento hay que hacerlo para alinear la convoluciòn.
xt1 = np.roll(xt1, -1)
xt1_shift = np.roll(xt1, m)

n = n + N // 2

x_pad = np.zeros(Naux)
x_pad[nmax - N//2: nmax - N//2 + N] = x
x1_pad = np.zeros(Naux)
x1_pad[nmax - N//2: nmax - N//2 + N] = xt1_shift[nmax - N//2: nmax - N//2 + N]


n = n[Npad:-Npad]
x_pad = x_pad[Npad:-Npad]
xt1 = xt1[Npad:-Npad]
xt1_shift = xt1_shift[Npad:-Npad]
x1_pad = x1_pad[Npad:-Npad]

fontsize = 12
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
plt.axis('off')

ax = plt.subplot2grid((4, 4), (1, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, xt1, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$\tilde{x}[n]=x[((n))_N]$', fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 0), rowspan=1, colspan=4)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, xt1_shift, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5)
plt.setp(bl, visible=False)
# labels
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(N, xtm, '$N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
plt.text(nmin_ax, xtit, r'$\tilde{{x}}_1[n]=\tilde{{x}}[n-m],\;m={}$'.format(m), fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')


ax3 = plt.subplot2grid((4, 4), (3, 0), rowspan=1, colspan=4)
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
plt.text(nmax_ax, xtm, '$n$', fontsize=fontsize, ha='right', va='baseline')
latex_tit = r'$x_1[n]=' \
            r'\left\{' \
            r'\begin{array}{ll}' \
            r'\tilde{x}_1[n], & 0\leq n\leq N-1,\\' \
            r'0, & \textrm{en otro caso.}' \
            r'\end{array}' \
            r'\right.$'

plt.text(nmin_ax, xtit, latex_tit, fontsize=fontsize, ha='left', va='baseline')
plt.axis('off')

con = ConnectionPatch(xyB=(-0.5, ymin), xyA=(-0.5, ymax), coordsA="data", coordsB="data",
                      axesB=ax3, axesA=ax0, color='k', lw=1, linestyle='dashed', zorder=100)
ax0.add_artist(con)
con = ConnectionPatch(xyB=(N-0.5, ymin), xyA=(N-0.5, ymax), coordsA="data", coordsB="data",
                      axesB=ax3, axesA=ax0, color='k', lw=1, linestyle='dashed', zorder=100)
ax0.add_artist(con)

# save as pdf image
plt.savefig('dft_properties_circular_shift.pdf', bbox_inches='tight')


plt.show()
