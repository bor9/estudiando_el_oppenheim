import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Parámetros
N = 15
n01 = 0
n02 = N / 2

n = np.arange(N)
h1 = (1 + np.cos(2 * np.pi / 15 * (n - n01))) / 15
h2 = (1 + np.cos(2 * np.pi / 15 * (n - n02))) / 15

nw = 512
w, H1 = signal.freqz(h1, a=1, worN=512, whole=False)
_, H2 = signal.freqz(h2, a=1, worN=512, whole=False)

magH1 = np.abs(H1)
magH2 = np.abs(H2)

argH1 = np.angle(H1)
argH2 = np.angle(H2)

# expresiones analíticas
Ha1 = np.exp(-1j * w * 7) * np.sin(15 * w / 2) * (1 / np.sin(w / 2)
                             + (np.exp(-1j * np.pi / 15) / 2) / np.sin((w - 2 * np.pi / 15) / 2)
                             + (np.exp(1j * np.pi / 15) / 2) / np.sin((w + 2 * np.pi / 15) / 2)
                             ) / 15
Ha2 = np.exp(-1j * w * 7) * np.sin(15 * w / 2) * (1 / np.sin(w / 2)
                             - (np.exp(-1j * np.pi / 15) / 2) / np.sin((w - 2 * np.pi / 15) / 2)
                             - (np.exp(1j * np.pi / 15) / 2) / np.sin((w + 2 * np.pi / 15) / 2)
                             ) / 15

Ha2_alt = - np.exp(-1j * 15 * w / 2) * (np.sin(np.pi / 15) ** 2) * np.sin(15 * w / 2) * np.cos(w / 2)\
            / (15 * np.sin(w / 2) * np.sin((w - 2 * np.pi / 15) / 2) * np.sin((w + 2 * np.pi / 15) / 2))

print(np.real(Ha2) - np.real(Ha2_alt))
print(np.imag(Ha2) - np.imag(Ha2_alt))

magHa1 = np.abs(Ha1)
magHa2 = np.abs(Ha2_alt)

argHa1 = np.angle(Ha1)
argHa2 = np.angle(Ha2_alt)

# figuras de prueba para la verificación de las expresiones analíticas.
fig = plt.figure(0)
ax = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)
plt.plot(w, magH1, 'k-')
plt.plot(w, magHa1, 'r-')
ax = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)
plt.plot(w, magH2, 'k-')
plt.plot(w, magHa2, 'r-')

fig = plt.figure(1)
ax = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)
plt.plot(w, argH1, 'k-')
plt.plot(w, argHa1, 'r-')
ax = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)
plt.plot(w, argH2, 'k-')
plt.plot(w, argHa2, 'r-')

########## Respuesta al impulso ##########

# parámetros de las figuras de la respuesta el impulso
dn = 2
nmin = -dn
nmax = N + 1
nmin_ax = nmin - 0.5
nmax_ax = nmax + 0.5

n = np.arange(nmin, nmax + 1)
h1_p = np.zeros(n.shape)
h1_p[dn: dn + N] = h1
h2_p = np.zeros(n.shape)
h2_p[dn: dn + N] = h2

xticks = np.arange(nmin, nmax + 1, 2)
ms = 4
fs = 10

fig = plt.figure(2, figsize=(6.5, 4), frameon=False)
ax = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)
plt.plot(n, h1_p, 'ks', ms=ms)
plt.xlim(nmin_ax, nmax_ax)
plt.text(nmin, h1[0]/2, r'$n_0={}$'.format(n01), ha='left', va='baseline')
ax.set_xticks(xticks)
plt.ylabel(r'$\mathrm{Amplitud}$')

ax = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)
plt.plot(n, h2_p, 'ks', ms=ms)
plt.xlim(nmin_ax, nmax_ax)
plt.text(nmin, h1[0]/2, r'$n_0=\dfrac{15}{2}$', ha='left', va='baseline')
ax.set_xticks(xticks)
plt.xlabel(r'$n$')
plt.ylabel(r'$\mathrm{Amplitud}$')

# save as pdf image
plt.savefig('structures_problem_06_39_impulse_response.pdf', bbox_inches='tight')

########## Respuesta en frecuencia ##########

dy = 0.5
ymax_ph = np.amax(np.concatenate((argH1, argH2))) + dy
ymin_ph = np.amin(np.concatenate((argH1, argH2))) - dy

fs = 11  # fontsize
y_label_coords = -0.14

xmin = 0
xmax = np.pi
xticks = 2 * np.pi * np.arange((N + 1) / 2) / 15
xticks_labels = ['$\dfrac{{{:.0f}\pi}}{{15}}$'.format(i) for i in 2 * np.arange((N + 1) / 2)]
xticks_labels[0] = '$0$'

yticks_ph = [-np.pi, -2, 0, 2]
yticks_labels_ph = ['$-\pi$', '$-2$', '$0$', '$2$']

w0 = 2 * np.pi / 15
_, H10 = signal.freqz(h1, a=1, worN=[w0])
_, H20 = signal.freqz(h2, a=1, worN=[w0])

fig = plt.figure(3, figsize=(8, 4), frameon=False)
# H1
# Respuesta en magnitud
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=2, colspan=2)
plt.plot(w, magH1, 'k', lw=1.5)
plt.plot([w0, w0], [0, np.abs(H10)[0]], 'k--', lw=1)
plt.plot([0, w0], [np.abs(H10)[0], np.abs(H10)[0]], 'k--', lw=1)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(xmin, 1.05)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.title(r'$n_0=0$', fontsize=fs)
plt.ylabel(r"$\textrm{Magnitud}$", fontsize=fs)
# Respuesta en fase
ax = plt.subplot2grid((4, 4), (2, 0), rowspan=2, colspan=2)
plt.plot(w, argH1, 'k', lw=1.5)
plt.plot([w0, w0], [ymin_ph, np.angle(H10)[0]], 'k--', lw=1)
plt.plot([0, w0], [np.angle(H10)[0], np.angle(H10)[0]], 'k--', lw=1)
plt.xlim(xmin, xmax)
plt.ylim(ymin_ph, ymax_ph)
ax.yaxis.set_label_coords(y_label_coords, 0.5)
plt.ylabel(r"$\textrm{Fase (rad)}$", fontsize=fs)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.yticks(yticks_ph, yticks_labels_ph, usetex=True)
plt.xlabel(r"$\omega$", fontsize=fs)

# H2
# Respuesta en magnitud
ax = plt.subplot2grid((4, 4), (0, 2), rowspan=2, colspan=2)
plt.plot(w, magH2, 'k', lw=1.5)
plt.plot([w0, w0], [0, np.abs(H20)[0]], 'k--', lw=1)
plt.plot([0, w0], [np.abs(H20)[0], np.abs(H10)[0]], 'k--', lw=1)
plt.xticks(xticks, [], usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(0, 1.05)
ax.yaxis.set_ticklabels([])
plt.title(r'$n_0=15/2$', fontsize=fs)
# Respuesta en fase
ax = plt.subplot2grid((4, 4), (2, 2), rowspan=2, colspan=2)
plt.plot(w, argH2, 'k', lw=1.5)
plt.plot([w0, w0], [ymin_ph, np.angle(H20)[0]], 'k--', lw=1)
plt.plot([0, w0], [np.angle(H20)[0], np.angle(H20)[0]], 'k--', lw=1)
plt.xlim(xmin, xmax)
plt.ylim(ymin_ph, ymax_ph)
plt.yticks(yticks_ph, [], usetex=True)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlabel(r"$\omega$", fontsize=fs)

# save as pdf image
plt.savefig('structures_problem_06_39_freq_response.pdf', bbox_inches='tight')

plt.show()
