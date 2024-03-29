import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# rango de n
nmax = 7
nmin = -1
n0 = -nmin

# vector n
n = np.arange(nmin, nmax + 1)

# tipo I
h1 = np.zeros(n.shape)
M1 = 4
h1[n0: n0 + M1 + 1] = signal.windows.triang(M1 + 1, sym=True)
# tipo II
h2 = np.zeros(n.shape)
M2 = 5
h2[n0: n0 + M2 + 1] = signal.windows.triang(M2 + 1, sym=True)
# tipo III y tipo IV
h3 = np.array([0, 0.5, 1, 0, -1, -0.5, 0, 0, 0])
h4 = np.array([0, 0.5, 1, -1, -0.5, 0, 0, 0, 0])


# rango de los ejes
ymax = 1.3
ymin = -1.1

delta_n = 1
nmin_ax = nmin - delta_n
nmax_ax = nmax + delta_n
delta_y = 0
ymax_ax = ymax + delta_y
ymin_ax = ymin - delta_y

baseline = -0.25
baseline_frac = -0.38
fontsize1 = 11
fontsize2 = 13
y_sep = 0.4

hmax = 1.2

fig = plt.figure(0, figsize=(8, 5), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, h1, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)
M = M1
plt.plot([M/2, M/2], [0, hmax], 'k--', lw=1)
plt.text(nmax_ax, baseline, '$n$', fontsize=fontsize1, ha='center', va='baseline')
plt.text(0, baseline, '$0$', fontsize=fontsize1, ha='center', va='baseline')
plt.text(M, baseline, '$M={}$'.format(M), fontsize=fontsize1, ha='center', va='baseline')
plt.text(M/2, baseline_frac, r'$\dfrac{M}{2}$', fontsize=fontsize1, ha='center', va='baseline')
plt.title(r'$\textrm{FIR tipo I: }M\textrm{ par, simetr\'ia par}$')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (0, 2), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, h2, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)
M = M2
plt.plot([M/2, M/2], [0, hmax], 'k--', lw=1)
plt.text(nmax_ax, baseline, '$n$', fontsize=fontsize1, ha='center', va='baseline')
plt.text(0, baseline, '$0$', fontsize=fontsize1, ha='center', va='baseline')
plt.text(M, baseline, '$M={}$'.format(M), fontsize=fontsize1, ha='center', va='baseline')
plt.text(M/2, baseline_frac, r'$\dfrac{M}{2}$', fontsize=fontsize1, ha='center', va='baseline')
plt.title(r'$\textrm{FIR tipo II: }M\textrm{ impar, simetr\'ia par}$')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 0), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, h3, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)
M = M1
plt.plot([M/2, M/2], [0, hmax], 'k--', lw=1)
plt.text(nmax_ax, baseline, '$n$', fontsize=fontsize1, ha='center', va='baseline')
plt.text(0, baseline, '$0$', fontsize=fontsize1, ha='center', va='baseline')
plt.text(M, 0.12, '$M={}$'.format(M), fontsize=fontsize1, ha='center', va='baseline')
plt.text(M/2, baseline_frac, r'$\dfrac{M}{2}$', fontsize=fontsize1, ha='center', va='baseline')
plt.title(r'$\textrm{FIR tipo III: }M\textrm{ par, simetr\'ia impar}$')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 2), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, h4, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)
M = 3
plt.plot([M/2, M/2], [0, hmax], 'k--', lw=1)
plt.text(nmax_ax, baseline, '$n$', fontsize=fontsize1, ha='center', va='baseline')
plt.text(0, baseline, '$0$', fontsize=fontsize1, ha='center', va='baseline')
plt.text(M, 0.12, '$M={}$'.format(M), fontsize=fontsize1, ha='center', va='baseline')
plt.text(M/2, baseline_frac, r'$\dfrac{M}{2}$', fontsize=fontsize1, ha='center', va='baseline')
plt.title(r'$\textrm{FIR tipo IV: }M\textrm{ impar, simetr\'ia impar}$')
plt.axis('off')

plt.savefig('transform_analysis_lineal_phase_fir_types.pdf', bbox_inches='tight')

## Respuestas en frecuencia

# número de frecuencias
nw = 512

w, H1 = signal.freqz(h1[n0:], 1, nw, whole=True)
magH1 = np.abs(H1)
_, H2 = signal.freqz(h2[n0:], 1, nw, whole=True)
magH2 = np.abs(H2)
_, H3 = signal.freqz(h3[n0:], 1, nw, whole=True)
magH3 = np.abs(H3)
_, H4 = signal.freqz(h4[n0:], 1, nw, whole=True)
magH4 = np.abs(H4)
# Se establecen los valores máximos y mínimos del eje y para todas las gráficas
dy = 0.4
ymax = np.amax(np.concatenate((magH1, magH2, magH3, magH3))) + dy
ymin = 0

y_label_coords = -0.1

xmin = 0
xmax = 2 * np.pi
xticks = np.linspace(xmin, xmax, 5)
xticks_labels = ['$0$', '$\dfrac{\pi}{2}$', '$\pi$', '$\dfrac{3\pi}{2}$', '$2\pi$']


fig = plt.figure(1, figsize=(8, 5), frameon=False)
ax = plt.subplot2grid((2, 4), (0, 0), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.plot(w, magH1, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$\textrm{FIR tipo I: }M\textrm{ par, simetr\'ia par}$')
plt.ylabel(r"$\textrm{Amplitud}$", fontsize=fontsize1)

ax = plt.subplot2grid((2, 4), (0, 2), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.plot(w, magH2, 'k-', lw=1.5)
plt.xticks(xticks, [], usetex=True)
ax.yaxis.set_ticklabels([])
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$\textrm{FIR tipo II: }M\textrm{ impar, simetr\'ia par}$')

ax = plt.subplot2grid((2, 4), (1, 0), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.plot(w, magH3, 'k-', lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$\textrm{FIR tipo III: }M\textrm{ par, simetr\'ia impar}$')
plt.xlabel(r"$\omega$", fontsize=fontsize1)
plt.ylabel(r"$\textrm{Amplitud}$", fontsize=fontsize1)

ax = plt.subplot2grid((2, 4), (1, 2), rowspan=1, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin, ymax)
plt.plot(w, magH4, 'k-', lw=1.5)
plt.xticks(xticks, xticks_labels, usetex=True)
ax.yaxis.set_ticklabels([])
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$\textrm{FIR tipo IV: }M\textrm{ impar, simetr\'ia impar}$')
plt.xlabel(r"$\omega$", fontsize=fontsize1)

plt.savefig('transform_analysis_lineal_phase_fir_types_freq_responses.pdf', bbox_inches='tight')

plt.show()

