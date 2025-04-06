import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack
from scipy import fft

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}', r'\usepackage{amssymb}']

# parámetros
# largo de la señal. número par.
N = 32
a = 0.95
w0 = 0.1 * np.pi
phi = 0

# construcción de la señal
n = np.arange(N)
x = (a ** n) * np.cos(w0 * n + phi)

# cálculo de la DFT y DCT-2
X = fft.fft(x)
Xc2 = fftpack.dct(x, type=2, norm=None)

# test DCT inversa. Hay que normalizar manualmente.
x_inv = fftpack.idct(Xc2, type=2, norm=None) / (2 * N)
print(np.sum(np.abs(x-x_inv)))

# Cálculo del error de truncamiento para la DFT
X_tr = X.copy()
ms_DFT = np.zeros(N // 2)
for i, m in enumerate(np.arange(1, N // 2 + 1)):
    X_tr[[N // 2 - m + 1, N // 2 + m - 1]] = 0
    x_tr = np.real(fft.ifft(X_tr))
    ms_DFT[i] = np.mean(np.square(x - x_tr))
m_DFT = np.arange(1, N // 2 + 1)
m_DFT[1:] = 2 * m_DFT[1:] - 1

# Cálculo del error de truncamiento para la DCT-2
Xc2_tr = Xc2.copy()
ms_DCT = np.zeros(N - 1)
for i, m in enumerate(np.arange(1, N)):
    Xc2_tr[-m] = 0
    x_tr = fftpack.idct(Xc2_tr, type=2) / (2 * N)
    ms_DCT[i] = np.mean(np.square(x - x_tr))
m_DCT = np.arange(1, N)

#
# Gráficas
#
fontsize = 11

# Señal en el tiempo

xmin_ax = n[0]
xmax_ax = n[-1]

fig = plt.figure(0, figsize=(8, 3), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
(markers, stemlines, bl) = plt.stem(n, x, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5, clip_on=False)
plt.setp(bl, visible=False)
plt.plot([xmin_ax, xmax_ax], [0, 0], 'k-', lw=1)

plt.xlabel('$n$', fontsize=fontsize)
plt.ylabel('$x[n]$', fontsize=fontsize)

plt.savefig('example_08_13_dct2_energy_compaction_x.pdf', bbox_inches='tight')

# DFT y DCT-2
K_dft = N // 2 + 1
k_dft = np.arange(K_dft)
k_dct = n

y_label_coords = -0.065
fig = plt.figure(1, figsize=(8, 6), frameon=False)
ax = plt.subplot2grid((3, 4), (0, 0), rowspan=1, colspan=4)
plt.xlim(k_dft[0], k_dft[-1])
(markers, stemlines, bl) = plt.stem(k_dft, np.real(X[:K_dft]), linefmt='k', markerfmt='sk',
                                    use_line_collection=True)
plt.setp(markers, markersize=3.5, clip_on=False)
plt.setp(bl, visible=False)
plt.plot([k_dft[0], k_dft[-1]], [0, 0], 'k-', lw=1)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.ylabel(r'$\operatorname{Re}\{X[k]\}$', fontsize=fontsize)

ax = plt.subplot2grid((3, 4), (1, 0), rowspan=1, colspan=4)
plt.xlim(k_dft[0], k_dft[-1])
(markers, stemlines, bl) = plt.stem(k_dft, np.imag(X[:K_dft]), linefmt='k', markerfmt='sk',
                                    use_line_collection=True)
plt.setp(markers, markersize=3.5, clip_on=False)
plt.setp(bl, visible=False)
plt.plot([k_dft[0], k_dft[-1]], [0, 0], 'k-', lw=1)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.ylabel(r'$\operatorname{Im}\{X[k]\}$', fontsize=fontsize)

ax = plt.subplot2grid((3, 4), (2, 0), rowspan=1, colspan=4)
plt.xlim(k_dct[0], k_dct[-1] + 1)
(markers, stemlines, bl) = plt.stem(k_dct, Xc2, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=3.5, clip_on=False)
plt.setp(bl, visible=False)
plt.plot([k_dct[0], k_dct[-1] + 1], [0, 0], 'k-', lw=1)
plt.xlabel('$k$', fontsize=fontsize)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.ylabel(r'$X^{c2}[k]$', fontsize=fontsize)

plt.savefig('example_08_13_dct2_energy_compaction_DFT_DCT2.pdf', bbox_inches='tight')

fig = plt.figure(2, figsize=(8, 3), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(m_DFT[0], m_DFT[-1])
plt.plot(m_DFT, ms_DFT, ls='-', marker='s', markersize=4, clip_on=False, zorder=10,
         label=r'$\textrm{Error de truncamiento de la DFT}$')
plt.plot(m_DCT, ms_DCT, ls='-', marker='s', markersize=4, clip_on=False, zorder=10,
         label=r'$\textrm{Error de truncamiento de la DCT-2}$')
leg = plt.legend(loc=2, frameon=False, fontsize=fontsize, framealpha=1)
plt.xlabel(r'$\textrm{N\'umero de coeficientes establecidos en cero}$', fontsize=fontsize)
plt.ylabel(r'$\textrm{Error MS}$', fontsize=fontsize)

plt.savefig('example_08_13_dct2_energy_compaction_MSE.pdf', bbox_inches='tight')

plt.show()



