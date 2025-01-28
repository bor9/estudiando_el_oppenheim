import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from matplotlib.patches import ConnectionPatch


__author__ = 'ernesto'

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


# parámetros. Largo del truncamiento y frecuencia de corte del pasabajos ideal
M = 12
wc = np.pi / 2

# rango de frecuencia
wmin = -np.pi
wmax = np.pi

# número de muestras en frecuencia. se elige así para que los cruces por cero
# coincidan con un muestra de frecuencia: nw = K * (M + 1) + 1.
nw = 40 * (M + 1) + 1
w = np.linspace(wmin, wmax, nw)
# muestra correspondiente a la frecuencia w=0
w0 = round((nw - 1) / 2)
# espectro de la ventana rectangular
W = np.sin(w * (M + 1) / 2) / np.sin(w / 2)
W[w0] = M + 1
# separación entre muestras de frecuencia
Dw = 2 * np.pi / (nw - 1)
# número de muestras entre 0 y pi/2. w[w0 + wpi2] = pi/2
wpi2 = round((np.pi / 2) / Dw)
# número de muestras entre 0 y el primer cruce por cero.
wpl = round(2 * np.pi / (M + 1) / Dw)
# espectro del filtro ideal
Hi = np.zeros(w.shape)
Hi[w0 - wpi2: w0 + wpi2 + 1] = 1

# Respuesta en frecuencia del filtro enventanado
n = np.arange(M/2 + 1)
hn = np.sin(wc * n) / (np.pi * n)
hn[0] = wc / np.pi
hn = np.concatenate((np.flip(hn[1:]), hn))
_, H = signal.freqz(hn, 1, worN=w[w0:], whole=False, plot=None, fs=2*np.pi, include_nyquist=True)
# se elimina el componente de fase lineal para que Hz sea real
H = np.real(H * np.exp(1j * w[w0:] * M / 2))
H = np.concatenate((np.flip(H[1:]), H))

# Pruebas
plt.figure(0)
plt.plot(w, Hi, 'k-')
plt.plot(w, H, 'r-')
plt.grid()

plt.figure(1)
ax = plt.subplot2grid((4, 1), (0, 0), rowspan=4, colspan=1)
W1 = np.roll(W, -(wpi2+wpl))
plt.plot(w, Hi, 'k-')
plt.plot(w, W1, 'r-')
ax.fill_between(w[w0 - wpi2: w0 + wpi2 + 1], W1[w0 - wpi2: w0 + wpi2 + 1])
plt.grid()

plt.figure(2)
plt.plot(w, W, 'r-')
plt.grid()

# Gráfica

# valores maximos y minimos de los ejes
dx = 0.6
xmax_ax = wmax
xmin_ax = wmin
ymax_ax1 = 1.1 * (M + 1)
ymin_ax1 = -0.3 * (M + 1)
ymax_ax2 = 1.2
ymin_ax2 = -0.2

# length of the ticks for all subplot (6 pixels)
display_length = 7  # in pixels
# x y ticks labels margin
xtm = -0.9
xtm2 = -1.5
ytm = -0.2
# font size
fontsize = 11
fontsize2 = 10
labpadx1 = 0
labpadx2 = -5
ytit = ymax_ax1 + 4.5

xticks = np.linspace(wmin, wmax, 5)
xticks_labels = ['$-\pi$', '$-\dfrac{\pi}{2}$', '$0$', '$\dfrac{\pi}{2}$', '$\pi$']

fig = plt.figure(3, figsize=(8.5, 8), frameon=False)

ax1b = plt.subplot2grid((14, 6), (0, 0), rowspan=3, colspan=2)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax1, ymax_ax1)

Wi = np.roll(W, -(wpi2 + wpl))
plt.plot(w, Hi, 'k-', lw=1.5, label=r'$H_d(e^{j\theta})$', zorder=20)
plt.plot(w, Wi, 'r-', zorder=30, label=r'$W(e^{j(\omega-\theta)})$')
ax1b.fill_between(w[w0 - wpi2: w0 + wpi2 + 1], Wi[w0 - wpi2: w0 + wpi2 + 1], color='salmon', zorder=10)
plt.xticks(xticks, [], usetex=True)
plt.legend(loc='upper right', fontsize=fontsize, frameon=False, framealpha=1)
plt.text(xmin_ax, ytit, '$\omega=-\pi/2-\Delta\omega_m/2$', ha='left', va='baseline', fontsize=fontsize)
ax1b.xaxis.set_label_position('top')
ax1b.xaxis.tick_top()
plt.xlabel(r'$\theta$', labelpad=labpadx1, fontsize=fontsize)

ax1a = plt.subplot2grid((14, 6), (3, 0), rowspan=3, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax2, ymax_ax2)
plt.plot(w, Hi, 'k-', lw=1.5, label=r'$H_d(e^{j\omega})$')
wi = w0 - (wpi2 + wpl)
plt.plot(w[0: wi], H[0: wi], 'r', label=r'$H(e^{j\omega})$')
plt.legend(loc='center', fontsize=fontsize, frameon=False, framealpha=1)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel(r'$\omega$', labelpad=labpadx2, fontsize=fontsize)

con = ConnectionPatch(xyA=(w[wi], H[wi]), xyB=(w[wi], Wi[wi]), coordsA="data", coordsB="data",
                      axesA=ax1a, axesB=ax1b, color='b', lw=1, linestyle='dashed', zorder=10)
ax1a.add_artist(con)
##################################################
##################################################
ax2b = plt.subplot2grid((14, 6), (0, 2), rowspan=3, colspan=2)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax1, ymax_ax1)

Wi = np.roll(W, -wpi2)
plt.plot(w, Hi, 'k-', lw=1.5, label=r'$H_d(e^{j\theta})$', zorder=20)
plt.plot(w, Wi, 'r-', zorder=30, label=r'$W(e^{j(\omega-\theta)})$')
ax2b.fill_between(w[w0 - wpi2: w0 + wpi2 + 1], Wi[w0 - wpi2: w0 + wpi2 + 1], color='salmon', zorder=10)
plt.xticks(xticks, [], usetex=True)
plt.text(xmin_ax, ytit, '$\omega=-\pi/2$', ha='left', va='baseline', fontsize=fontsize)
ax2b.xaxis.set_label_position('top')
ax2b.xaxis.tick_top()
plt.xlabel(r'$\theta$', labelpad=labpadx1, fontsize=fontsize)
ax2b.set_yticklabels([])

ax2a = plt.subplot2grid((14, 6), (3, 2), rowspan=3, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax2, ymax_ax2)
plt.plot(w, Hi, 'k-', lw=1.5, label=r'$H_d(e^{j\omega})$')
wi = w0 - wpi2
plt.plot(w[0: wi], H[0: wi], 'r', label=r'$H(e^{j\omega})$')
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel(r'$\omega$', labelpad=labpadx2, fontsize=fontsize)
ax2a.set_yticklabels([])

con = ConnectionPatch(xyA=(w[wi], H[wi]), xyB=(w[wi], Wi[wi]), coordsA="data", coordsB="data",
                      axesA=ax2a, axesB=ax2b, color='b', lw=1, linestyle='dashed', zorder=10)
ax2a.add_artist(con)
##################################################
##################################################
ax3b = plt.subplot2grid((14, 6), (0, 4), rowspan=3, colspan=2)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax1, ymax_ax1)

Wi = np.roll(W, -wpi2 + wpl)
plt.plot(w, Hi, 'k-', lw=1.5, label=r'$H_d(e^{j\theta})$', zorder=20)
plt.plot(w, Wi, 'r-', zorder=30, label=r'$W(e^{j(\omega-\theta)})$')
ax3b.fill_between(w[w0 - wpi2: w0 + wpi2 + 1], Wi[w0 - wpi2: w0 + wpi2 + 1], color='salmon', zorder=10)
plt.xticks(xticks, [], usetex=True)
plt.text(xmin_ax, ytit, '$\omega=-\pi/2+\Delta\omega_m/2$', ha='left', va='baseline', fontsize=fontsize)
ax3b.xaxis.set_label_position('top')
ax3b.xaxis.tick_top()
plt.xlabel(r'$\theta$', labelpad=labpadx1, fontsize=fontsize)
ax3b.set_yticklabels([])

ax3a = plt.subplot2grid((14, 6), (3, 4), rowspan=3, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax2, ymax_ax2)
plt.plot(w, Hi, 'k-', lw=1.5, label=r'$H_d(e^{j\omega})$')
wi = w0 - wpi2 + wpl
plt.plot(w[0: wi], H[0: wi], 'r', label=r'$H(e^{j\omega})$')
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel(r'$\omega$', labelpad=labpadx2, fontsize=fontsize)
ax3a.set_yticklabels([])

con = ConnectionPatch(xyA=(w[wi], H[wi]), xyB=(w[wi], Wi[wi]), coordsA="data", coordsB="data",
                      axesA=ax3a, axesB=ax3b, color='b', lw=1, linestyle='dashed', zorder=10)
ax3a.add_artist(con)
##################################################
##################################################
ax4b = plt.subplot2grid((14, 6), (8, 0), rowspan=3, colspan=2)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax1, ymax_ax1)

Wi = np.roll(W, wpi2 - wpl)
plt.plot(w, Hi, 'k-', lw=1.5, label=r'$H_d(e^{j\theta})$', zorder=20)
plt.plot(w, Wi, 'r-', zorder=30, label=r'$W(e^{j(\omega-\theta)})$')
ax4b.fill_between(w[w0 - wpi2: w0 + wpi2 + 1], Wi[w0 - wpi2: w0 + wpi2 + 1], color='salmon', zorder=10)
plt.xticks(xticks, [], usetex=True)
plt.text(xmin_ax, ytit, '$\omega=\pi/2-\Delta\omega_m/2$', ha='left', va='baseline', fontsize=fontsize)
ax4b.xaxis.set_label_position('top')
ax4b.xaxis.tick_top()
plt.xlabel(r'$\theta$', labelpad=labpadx1, fontsize=fontsize)


ax4a = plt.subplot2grid((14, 6), (11, 0), rowspan=3, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax2, ymax_ax2)
plt.plot(w, Hi, 'k-', lw=1.5, label=r'$H_d(e^{j\omega})$')
wi = w0 + wpi2 - wpl
plt.plot(w[0: wi], H[0: wi], 'r', label=r'$H(e^{j\omega})$')
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel(r'$\omega$', labelpad=labpadx2, fontsize=fontsize)

con = ConnectionPatch(xyA=(w[wi], H[wi]), xyB=(w[wi], Wi[wi]), coordsA="data", coordsB="data",
                      axesA=ax4a, axesB=ax4b, color='b', lw=1, linestyle='dashed', zorder=10)
ax4a.add_artist(con)
##################################################
##################################################
ax5b = plt.subplot2grid((14, 6), (8, 2), rowspan=3, colspan=2)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax1, ymax_ax1)

Wi = np.roll(W, wpi2)
plt.plot(w, Hi, 'k-', lw=1.5, label=r'$H_d(e^{j\theta})$', zorder=20)
plt.plot(w, Wi, 'r-', zorder=30, label=r'$W(e^{j(\omega-\theta)})$')
ax5b.fill_between(w[w0 - wpi2: w0 + wpi2 + 1], Wi[w0 - wpi2: w0 + wpi2 + 1], color='salmon', zorder=10)
plt.xticks(xticks, [], usetex=True)
plt.text(xmin_ax, ytit, '$\omega=\pi/2$', ha='left', va='baseline', fontsize=fontsize)
ax5b.xaxis.set_label_position('top')
ax5b.xaxis.tick_top()
plt.xlabel(r'$\theta$', labelpad=labpadx1, fontsize=fontsize)
ax5b.set_yticklabels([])

ax5a = plt.subplot2grid((14, 6), (11, 2), rowspan=3, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax2, ymax_ax2)
plt.plot(w, Hi, 'k-', lw=1.5, label=r'$H_d(e^{j\omega})$')
wi = w0 + wpi2
plt.plot(w[0: wi], H[0: wi], 'r', label=r'$H(e^{j\omega})$')
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel(r'$\omega$', labelpad=labpadx2, fontsize=fontsize)
ax5a.set_yticklabels([])

con = ConnectionPatch(xyA=(w[wi], H[wi]), xyB=(w[wi], Wi[wi]), coordsA="data", coordsB="data",
                      axesA=ax5a, axesB=ax5b, color='b', lw=1, linestyle='dashed', zorder=10)
ax5a.add_artist(con)
##################################################
##################################################
ax6b = plt.subplot2grid((14, 6), (8, 4), rowspan=3, colspan=2)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax1, ymax_ax1)

Wi = np.roll(W, wpi2 + wpl)
plt.plot(w, Hi, 'k-', lw=1.5, label=r'$H_d(e^{j\theta})$', zorder=20)
plt.plot(w, Wi, 'r-', zorder=30, label=r'$W(e^{j(\omega-\theta)})$')
ax6b.fill_between(w[w0 - wpi2: w0 + wpi2 + 1], Wi[w0 - wpi2: w0 + wpi2 + 1], color='salmon', zorder=10)
plt.xticks(xticks, [], usetex=True)
plt.text(xmin_ax, ytit, '$\omega=\pi/2+\Delta\omega_m/2$', ha='left', va='baseline', fontsize=fontsize)
ax6b.xaxis.set_label_position('top')
ax6b.xaxis.tick_top()
plt.xlabel(r'$\theta$', labelpad=labpadx1, fontsize=fontsize)
ax6b.set_yticklabels([])

ax6a = plt.subplot2grid((14, 6), (11, 4), rowspan=3, colspan=2)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax2, ymax_ax2)
plt.plot(w, Hi, 'k-', lw=1.5, label=r'$H_d(e^{j\omega})$')
wi = w0 + wpi2 + wpl
plt.plot(w[0: wi], H[0: wi], 'r', label=r'$H(e^{j\omega})$')
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel(r'$\omega$', labelpad=labpadx2, fontsize=fontsize)
ax6a.set_yticklabels([])

con = ConnectionPatch(xyA=(w[wi], H[wi]), xyB=(w[wi], Wi[wi]), coordsA="data", coordsB="data",
                      axesA=ax6a, axesB=ax6b, color='b', lw=1, linestyle='dashed', zorder=10)
ax6a.add_artist(con)
##################################################
##################################################
# save as pdf image
plt.savefig('filter_design_windowing_convolution_ripple.pdf', bbox_inches='tight')

plt.show()
