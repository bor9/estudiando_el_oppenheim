import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Filtro butter pasabajos
# Orden y frecuencia de corte
ord = 2
wc = 0.2 # frecuencia de corte, 0 < wc < 1
b1, a1 = signal.butter(ord, wc)
# Frecuencia específica
omega0 = 0.75

# Polos y ceros del filtro
zeros1 = np.roots(b1)
poles1 = np.roots(a1)

# Se agrega un par de polos y ceros de un filtro pasatodos para cambiar la respuesta en fase  y
# así hacer que el retardo de grupo sea mas grande.
pn = 0.9 * np.exp(1j * np.pi / 10)
zeros = np.append(zeros1, [1/pn, np.conj(1/pn)])
poles = np.append(poles1, [pn, np.conj(pn)])
b = np.poly(zeros)
a = np.poly(poles)

# Normalización para ganancia 1 en continua
G = np.polyval(b, 1) / np.polyval(a, 1)
a = G * a

# Respuesta en frecuencia y retardo de grupo del filtro
ntheta = 512
# respuesta en frecuencia
w_half, H = signal.freqz(b, a, ntheta)
# magnitud de la respuesta en frecuencia
magH_half = np.abs(H)
# fase de la respuesta en frecuencia
phiH_half = np.unwrap(np.angle(H))
# retardo de grupo de la respuesta en frecuencia
grdH = signal.group_delay((b, a), w_half)

# construcción de la respuesta en frecuencia entre -pi y pi
w = np.append(-np.flip(w_half[1:]), w_half)
magH = np.append(np.flip(magH_half[1:]), magH_half)
phiH = np.append(-np.flip(phiH_half[1:]), phiH_half)

# evaluación de la respuesta en frecuencia en theta_0
_, H_omega0 = signal.freqz(b, a, omega0)
G0 = np.abs(H_omega0[0])
# lo siguiente hay que hacerlo porque la fase esta desenvuelta
# indice del valor mas cercano a theta0 en w_half
idx_w0 = (np.abs(w_half - omega0)).argmin()
theta0 = np.angle(H_omega0[0]) + int(phiH_half[idx_w0] / np.pi) * np.pi

# Figuras
fs = 11
xmax = np.pi
xmin = -np.pi
xticks = [-np.pi, 0, omega0, np.pi]
xticks_labels = ['$-\pi$', '$0$', r'$\omega_0$', '$\pi$']
ind_theta0 = 2
ylabel_xcoord = -0.06
ylabel_coords = -0.05

fig = plt.figure(0, figsize=(7, 4), frameon=False)
ax = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)
ymax = 1.1
plt.xlim(xmin, xmax)
plt.ylim(0, ymax)
plt.plot(w, magH, lw=2, color='k')
plt.plot([0, 0], [0, ymax], lw=1, color='k')
plt.yticks([0, 0.5, 1])
plt.xticks(xticks, xticks_labels)
plt.ylabel(r"$\rm Magnitud$", fontsize=fs)
# se alinea la etiqueta del eje y de ambas graficas
ax.yaxis.set_label_coords(ylabel_xcoord, 0.5)
# se cambia el color de la etiqueta theta0
labels = plt.gca().get_xticklabels()
labels[ind_theta0].set_color('red')
# ganancia en theta0
plt.plot([omega0, omega0], [0, G0], 'r', lw=1)
plt.plot([0, omega0], [G0, G0], 'r', lw=1)
plt.text(ylabel_coords, G0, '$G_0$', fontsize=fs, ha='right', va='center', color='r')
# etiqueta de la gráfica
plt.text(0.98, 0.85, '$|H(e^{j\omega})|$', fontsize=fs, ha='right', va='baseline', transform = ax.transAxes)

ax = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)
plt.plot(w, phiH, lw=2, color='k')
plt.xlim(xmin, xmax)
ymin, ymax = ax.get_ylim()
plt.ylim(ymin, ymax)
plt.plot([xmin, xmax], [0, 0], lw=1, color='k')
plt.plot([0, 0], [ymin, ymax], lw=1, color='k')
plt.xticks(xticks, xticks_labels)
plt.xlabel(r"$\rm Frecuencia\;(rad)$", fontsize=fs)
plt.ylabel(r"$\rm Fase\;(rad)$", fontsize=fs)
# se alinea la etiqueta del eje y de ambas graficas
ax.yaxis.set_label_coords(ylabel_xcoord, 0.5)
# se cambia el color de la etiqueta theta_0
labels = plt.gca().get_xticklabels()
labels[ind_theta0].set_color('red')
# ganancia en theta0
plt.plot([omega0, omega0], [ymin, theta0], 'r', lw=1)
plt.plot([0, omega0], [theta0, theta0], 'r', lw=1)
plt.text(ylabel_coords, theta0, r'$\theta_0$', fontsize=fs, ha='right', va='center', color='r')
# etiqueta de la gráfica
plt.text(0.98, 0.85, r'$\angle H(e^{j\omega})$', fontsize=fs, ha='right', va='baseline', transform = ax.transAxes)

plt.savefig('seq_and_sys_frequency_response.pdf', bbox_inches='tight')
plt.show()
