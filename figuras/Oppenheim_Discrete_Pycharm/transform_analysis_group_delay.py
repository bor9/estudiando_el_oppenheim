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
ord = 4
wc = 0.5 # frecuencia de corte, 0 < wc < 1
b1, a1 = signal.butter(ord, wc)

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
nw = 1024
# respuesta en frecuencia
w_half, H = signal.freqz(b, a, nw)
# magnitud de la respuesta en frecuencia
magH_half = np.abs(H)
# fase de la respuesta en frecuencia
argH_half = np.unwrap(np.angle(H))
# retardo de grupo de la respuesta en frecuencia
_, grdH = signal.group_delay((b, a), w_half)

# Filtrado de sinusoide enventanada para ver el retardo de grupo
# Señal de entrada
nx = 200
# largo del relleno ceros a ambos lados
nz = 30
# Frecuencia de la portadora
w0 = np.pi / 10
x = np.cos(w0 * np.arange(nx)) * np.hanning(nx)
x = np.pad(x, (nz, nz), 'constant')
nmax = nx + 2 * nz
n = np.arange(nmax)
# Señal de salida
y = signal.lfilter(b, a, x)
# Espectro de la señal de entrada
_, X = signal.freqz(x, 1, w_half)
magX = np.abs(X)

# evaluación de la respuesta en frecuencia en las frecuencias de las sinusoides
_, H_theta_x = signal.freqz(b, a, w0)
# evaluación de la fase, el retardo de fase y el retardo de grupo en la frecuencia de la sinusoide
# retardo de fase
# lo siguiente hay que hacerlo porque la fase esta desenvuelta:
# indice del valor mas cercano a theta_x en w_half
idx = (np.abs(w_half - w0)).argmin()
# fase desenvuelta en esa frecuencia.
ph_0 = np.angle(H_theta_x) + np.rint(argH_half[idx] / (2 * np.pi)) * 2 * np.pi
# retardo de fase
tau_ph_0 = -ph_0 / w0
# Retardo de grupo
_, tau_gr_0 = signal.group_delay((b, a), w0)
# esto hay que hacerlo porque las funciones devuelven arreglos.
ph_0 = ph_0[0]
tau_ph_0 = tau_ph_0[0]
tau_gr_0 = tau_gr_0[0]

########## Gráficas ##########

fs = 11  # fontsize

# Señal de entrada y respuesta en fase con linealización

xmax = np.pi / 4
xmin = 0
xticks = [0, w0, np.pi/4]
xticks_labels = ['$0$', '$\omega_0$', '$\pi/4$']
yticks = [0]
yticks_labels = ['$0$']

fig = plt.figure(0, figsize=(8, 5), frameon=False)
ax = plt.subplot2grid((3, 1), (0, 0), rowspan=1, colspan=1)
plt.plot(w_half, magX, lw=1.5, color='k')
plt.xticks(xticks, xticks_labels, usetex=True)
plt.yticks(yticks, yticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(0, 60)
plt.title(r'${\rm Magnitud\;del\;espectro\;de\;la\;entrada\;de\;banda\;angosta}\;(|X(e^{j\omega})|)$', fontsize=fs)
plt.text(xmax / 2, -10, r"$\omega\;\rm(rad)$", fontsize=fs, ha='left', va='center')

# linealización en torno a omega_0
phi0 = w0 * (tau_ph_0 - tau_gr_0)
nd = tau_gr_0
dw = 0.4

ymin = -7
yticks = [0, ph_0]
yticks_labels = ['$0$', r'$\angle H(e^{j\omega_0})$']

ax = plt.subplot2grid((3, 1), (1, 0), rowspan=2, colspan=1)
plt.subplots_adjust(hspace=0.3)
plt.plot(w_half, argH_half, lw=1.5, color='k')
plt.plot([w0 - dw, w0 + dw], [-phi0 - nd * (w0 - dw), -phi0 - nd * (w0 + dw)], lw=1.5, color=[0.7, 0.7, 0.7])
plt.plot([0, w0], [ph_0, ph_0], 'k--', lw=1)
plt.plot([w0, w0], [0, ph_0], 'k--', lw=1)
plt.xticks(xticks, [])
plt.yticks(yticks, yticks_labels, usetex=True)
plt.xlim(xmin, xmax)
plt.ylim(ymin, 0)
ax.xaxis.tick_top()
plt.text(w0 + 0.03, ymin + 0.4, r"$\angle H(e^{j\omega})\approx-\phi_0-\omega_0n_d$", fontsize=fs, ha='center',
         va='bottom')
plt.text(xmax / 2, ymin - 0.4, r"${\rm Respuesta\;en\;fase\;}(\arg H(e^{j\omega}))\rm{\;y\;linealizaci\acute{o}n}$",
         fontsize=fs, ha='center', va='center')
plt.savefig('transform_analysis_group_delay_input_and_linealization.pdf', bbox_inches='tight')


# Respuesta en fase, retardo de fase y retardo de grupo
fs_t = 10
ylabel_xcoord = -0.05
xmax = np.pi
xmin = 0

xticks = [0, w0, np.pi/2, np.pi]
xticks_labels = ['$0$', '$\omega_0$', '$\pi/2$', '$\pi$']

# Retardo de fase
tauPhi_half = -argH_half / w_half

fig = plt.figure(1, figsize=(9, 5), frameon=False)
ax = plt.subplot2grid((6, 1), (0, 0), rowspan=2, colspan=1)
plt.plot(w_half, argH_half, lw=2, color='k')
plt.xlim(xmin, xmax)
ymin, ymax = ax.get_ylim()
plt.ylim(ymin, 0)
plt.xticks(xticks, [], usetex=True)
plt.ylabel(r"$\rm Fase\;(rad)$", fontsize=fs_t)
# se alinea la etiqueta del eje y de ambas graficas
ax.yaxis.set_label_coords(ylabel_xcoord, 0.5)
# fase en w0
plt.plot([w0, w0], [0, ph_0], 'k--', lw=1)
plt.plot([0, w0], [ph_0, ph_0], 'k--', lw=1)
# etiqueta de la gráfica
plt.text(0.98, 0.8, r'$\angle H(e^{j\omega})$', fontsize=12, ha='right', va='baseline', transform = ax.transAxes)

ax = plt.subplot2grid((6, 1), (2, 0), rowspan=2, colspan=1)
plt.plot(w_half, tauPhi_half, lw=2, color='k')
plt.xlim(xmin, xmax)
ymin, ymax = ax.get_ylim()
ymin = 0
plt.ylim(0, ymax)
plt.xticks(xticks, [], usetex=True)
plt.ylabel(r"$\rm Retardo\;(muestras)$", fontsize=fs_t)
# se alinea la etiqueta del eje y de ambas graficas
ax.yaxis.set_label_coords(ylabel_xcoord, 0.5)
# retardo de fase en w0
plt.plot([w0, w0], [0, tau_ph_0], 'k--', lw=1)
plt.plot([0, w0], [tau_ph_0, tau_ph_0], 'k--', lw=1)
plt.text(0.98, 0.8, r'$\tau_\textrm{ph}(\omega)=-\dfrac{\angle H(e^{j\omega})}{\omega}$', fontsize=12, ha='right',
         va='center', transform = ax.transAxes)

ax = plt.subplot2grid((6, 1), (4, 0), rowspan=2, colspan=1)
plt.plot(w_half, grdH, lw=2, color='k')
plt.xlim(xmin, xmax)
ymin, ymax = ax.get_ylim()
ymin = 0
plt.ylim(0, ymax)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.ylabel(r"$\rm Retardo\;(muestras)$", fontsize=fs_t)
# se alinea la etiqueta del eje y de ambas graficas
ax.yaxis.set_label_coords(ylabel_xcoord, 0.5)
# retardo de grupo en w0
plt.plot([w0, w0], [0, tau_gr_0], 'k--', lw=1)
plt.plot([0, w0], [tau_gr_0, tau_gr_0], 'k--', lw=1)
plt.xlabel(r"$\omega\;\rm(rad)$", fontsize=fs)
# etiqueta de la gráfica
plt.text(0.98, 0.8, r'$\tau_\textrm{gr}(\omega)=-\dfrac{d}{d\omega}\angle H(e^{j\omega})$', fontsize=12, ha='right',
         va='center', transform = ax.transAxes)
plt.savefig('transform_analysis_group_delay_phase_and_group_delay.pdf', bbox_inches='tight')


# Señales de entrada y salida del filtro
ms = 2.5
ymin = -1.1
ymax = 1.6

fig = plt.figure(2, figsize=(8, 5), frameon=False)
ax = fig.add_subplot(111)
plt.plot(n, x, ls='-', marker='s', color='k', markersize=ms, lw=1, label='${\\rm entrada}$')
plt.plot(np.arange(nx) + nz, np.hanning(nx), 'k--', lw=1)
plt.plot(n, y, ls='-', marker='s', color='r', markersize=ms, lw=1, label='${\\rm salida}$')
plt.plot(np.arange(nx) + nz + tau_gr_0, np.hanning(nx), 'r--', lw=1)
# Retardo de fase y retardo de grupo
hl1 = 1.1
dl = 0.03
nper = 5
ns = nz + np.pi * 2 * nper / w0
plt.plot([ns, ns + tau_ph_0], [hl1, hl1], 'k-', lw = 1)
plt.plot([ns, ns], [hl1 - dl, hl1 + dl], 'k-', lw = 1)
plt.plot([ns + tau_ph_0, ns + tau_ph_0], [hl1 - dl, hl1 + dl], 'k-', lw = 1)
plt.text(ns, hl1 + 0.05, r'$\tau_\textrm{{ph}}(\omega_0)\approx{:.2f}$'.format(tau_ph_0), fontsize=fs, ha='left',
         va='bottom')
hl2 = 1.35
plt.plot([ns, ns + tau_gr_0], [hl2, hl2], 'k-', lw = 1)
plt.plot([ns, ns], [hl2 - dl, hl2 + dl], 'k-', lw = 1)
plt.plot([ns + tau_gr_0, ns + tau_gr_0], [hl2 - dl, hl2 + dl], 'k-', lw = 1)
plt.text(ns, hl2 + 0.05, r'$\tau_\textrm{{gr}}(\omega_0)\approx{:.2f}$'.format(tau_gr_0), fontsize=12, ha='left',
         va='bottom')
plt.xlim(0, nmax)
plt.ylim(ymin, ymax)
plt.xlabel(r"$n$", fontsize=fs)
plt.ylabel(r"$\textrm{Amplitud}$", fontsize=fs)
plt.title(r"$\textrm{Entrada }x[n]=s[n]\cos\omega_0n\textrm{ y la correspondiente salida}$", fontsize=fs)

leg = plt.legend(loc=1, frameon=False, fontsize=fs, framealpha=1)

plt.savefig('transform_analysis_group_delay_input_output.pdf', bbox_inches='tight')

plt.show()
