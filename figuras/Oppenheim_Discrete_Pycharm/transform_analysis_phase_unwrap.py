import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from matplotlib.collections import LineCollection

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Parámetros de filtro eliptico pasabajos
# Orden y frecuencia de corte
ord = 6
rp = 0.1
rs = 30
wc = 0.3 # frecuencia de corte, 0 < wc < 1

# Construccion del filtro FIR mediante enventanado
# Filtro elíptico
b, a = signal.ellip(ord, rp, rs, wc)
# Obtención de la respuesta al impulso
delta = np.zeros(100)
delta[1] = 1
h = signal.lfilter(b, a, delta)
# Ventana
M = 30
win = np.hanning(2 * M + 1)
# Enventanado de la respuesta al impulso
h = h[:M + 1] * win[M:]

# Respuesta en frecuencia y retardo de grupo del filtro
ntheta = 1024
# Respuesta en frecuencia
w, H = signal.freqz(h, 1, ntheta)
# Magnitud de la respuesta en frecuencia
magH = np.abs(H)
# Fase de la respuesta en frecuencia
phiH = np.angle(H)
phiH_unw = np.unwrap(phiH)

# Para no dibujar la rectas verticales en los saltos de fase
phiH_der = phiH[1:] - phiH[:-1]
phase_jumps = np.argwhere(np.abs(phiH[1:] - phiH[:-1]) > np.pi).flatten()
phiH[phase_jumps] = 'nan'
rw = np.rint(phiH_unw / (2 * np.pi))
rw[phase_jumps] = 'nan'

# Para líneas de gráficas con colormap
cols = np.linspace(0, 1, len(H))

points1 = np.array([w, phiH]).T.reshape(-1, 1, 2)
segments1 = np.concatenate([points1[:-1], points1[1:]], axis=1)

points2 = np.array([w, phiH_unw]).T.reshape(-1, 1, 2)
segments2 = np.concatenate([points2[:-1], points2[1:]], axis=1)

points3 = np.array([H.real, H.imag]).T.reshape(-1, 1, 2)
segments3 = np.concatenate([points3[:-1], points3[1:]], axis=1)

cmap = 'viridis'

# Ticks
xticks = [0, np.pi/2, np.pi]
xtickslabels = ['$0$', '$\pi/2$', '$\pi$']

yticks1 = [-np.pi, 0, np.pi]
ytickslabels1 = ['$-\pi$', '$0$', '$\pi$']

yticks2 = np.pi * np.arange(-4, 1)
ytickslabels2 = ['$-4\pi$', '$-3\pi$','$-2\pi$', '$-\pi$', '$0$']

ylabel_xcoord = -0.07

fig = plt.figure(0, figsize=(7, 5), frameon=False)

ax1 = plt.subplot2grid((4, 1), (0, 0), rowspan=1, colspan=1)
lc1 = LineCollection(segments1, cmap=cmap)
lc1.set_array(cols)
lc1.set_linewidth(2)
line = ax1.add_collection(lc1)
ax1.set(xlim=(0, np.pi), ylim=(-np.pi, np.pi))
plt.xticks(xticks, '')
plt.yticks(yticks1, ytickslabels1)
plt.grid()
plt.ylabel(r'$\textrm{ARG}[H(e^{j\omega})]$')
ax1.yaxis.set_label_coords(ylabel_xcoord, 0.5)

ax2 = plt.subplot2grid((4, 1), (1, 0), rowspan=2, colspan=1)
lc2 = LineCollection(segments2, cmap=cmap)
lc2.set_array(cols)
lc2.set_linewidth(2)
line = ax2.add_collection(lc2)
ax2.set(xlim=(0, np.pi), ylim=(-4 * np.pi, 0))
plt.xticks(xticks, '')
plt.yticks(yticks2, ytickslabels2)
plt.grid()
plt.ylabel(r'$\textrm{arg}[H(e^{j\omega})]$')
ax2.yaxis.set_label_coords(ylabel_xcoord, 0.5)

ax3 = plt.subplot2grid((4, 1), (3, 0), rowspan=1, colspan=1)
plt.plot(w, rw, 'k', lw=2)
ax3.set(xlim=(0, np.pi))
plt.xticks(xticks, xtickslabels)
plt.grid()
plt.xlabel(r'$\omega\textrm{ (rad)}$')
plt.ylabel(r'$r(\omega)$')
ax3.yaxis.set_label_coords(ylabel_xcoord, 0.5)

plt.savefig('transform_analysis_phase_unwrap.pdf', bbox_inches='tight')

fig = plt.figure(1, figsize=(7, 6), frameon=False)
ax4 = fig.subplots()
ax4.set_aspect('equal', adjustable='box')
ax4.set(xlim=(-1, 1), ylim=(-1, 1))
lc3 = LineCollection(segments3, cmap=cmap)
lc3.set_array(cols)
lc3.set_linewidth(1.5)
line = ax4.add_collection(lc3)
plt.grid()
# puntos específicos
ms = 6
plt.plot(H[0].real, H[0].imag, 'k.', markersize=ms)
plt.annotate(r'$\omega=0$', xytext=(H[0].real - 0.06, H[0].imag + 0.1), xycoords='data', xy=(H[0].real, H[0].imag),
             textcoords='data', color='k', fontsize=10, va="baseline", ha="right",
             arrowprops=dict(arrowstyle="-|>, head_width=0.2, head_length=0.7", color='k', relpos=(1, 0),
                             patchA=None, patchB=None, shrinkA=1, shrinkB=1))
plt.plot(H[-1].real, H[-1].imag, 'k.', markersize=ms)
plt.annotate(r'$\omega=\pi$', xytext=(H[-1].real - 0.06, H[-1].imag - 0.12), xycoords='data', xy=(H[-1].real, H[-1].imag),
             textcoords='data', color='k', fontsize=10, va="baseline", ha="right",
             arrowprops=dict(arrowstyle="-|>, head_width=0.2, head_length=0.7", color='k', relpos=(1, 1),
                             patchA=None, patchB=None, shrinkA=1, shrinkB=1))
plt.xlabel(r'$\textrm{Re}[H(e^{j\omega})]$')
plt.ylabel(r'$\textrm{Im}[H(e^{j\omega})]$')
ax4.set(xlim=(-0.9, 1), ylim=(-1, 0.6))

plt.savefig('transform_analysis_phase_unwrap_z_plane.pdf', bbox_inches='tight')

plt.show()
