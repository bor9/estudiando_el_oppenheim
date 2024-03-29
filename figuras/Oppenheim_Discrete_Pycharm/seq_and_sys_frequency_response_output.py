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

# frecuencias de sinusoides a filtrar
thetas = [0.25, 0.75, 1.5]
nt = len(thetas)

# evaluación de la respuesta en frecuencia en las frecuencias de las sinusoides
_, H_thetas = signal.freqz(b, a, thetas)
Gs = np.abs(H_thetas)
# lo siguiente hay que hacerlo porque la fase esta desenvuelta:
H_args = np.zeros(Gs.shape)
for i in np.arange(nt):
    # indice del valor mas cercano a thetas[i] en w_half
    idx = (np.abs(w_half - thetas[i])).argmin()
    # fase desenvuelta en esa frecuencia.
    H_args[i] = np.angle(H_thetas[i]) + int(phiH_half[idx] / np.pi) * np.pi

# GRÁFICAS DE RESPUESTA EN FRECUENCIA
fs = 11
xmax = np.pi
xmin = 0
xticks = np.concatenate(([0], thetas, [np.pi]), axis=0)
xticks_labels = ["${:.2f}$".format(xt) if xt != 0 else "${:.0f}$".format(xt) for xt in xticks]
xticks_labels[-1] = "$\pi$"
ylabel_xcoord = -0.1
ylabel_coords = -0.05

fig = plt.figure(0, figsize=(8, 5), frameon=False)
ax = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)
ymax = 1.1
plt.xlim(xmin, xmax)
plt.ylim(0, ymax)
plt.plot(w_half, magH_half, lw=2, color='k')
plt.plot([0, 0], [0, ymax], lw=1, color='k')
plt.xticks(xticks, xticks_labels, usetex=True)
yticks = np.concatenate(([0], Gs))
yticks_labels = ["${:.4f}$".format(yt) if yt != 0 else "${:.0f}$".format(yt) for yt in yticks]
plt.yticks(yticks, yticks_labels)
plt.ylabel(r"$\rm Magnitud$", fontsize=fs)
# se alinea la etiqueta del eje y de ambas graficas
ax.yaxis.set_label_coords(ylabel_xcoord, 0.5)
# ganancia en thetas
plt.plot(np.array([thetas, thetas]), np.array([np.zeros(Gs.shape), Gs]), 'k', lw=1)
plt.plot(np.array([np.zeros(Gs.shape), thetas]), np.array([Gs, Gs]), 'k', lw=1)
# etiqueta de la gráfica
plt.text(0.98, 0.85, '$|H(e^{j\omega})|$', fontsize=12, ha='right', va='baseline', transform = ax.transAxes)

ax = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)
plt.plot(w_half, phiH_half, lw=2, color='k')
plt.xlim(xmin, xmax)
ymin, ymax = ax.get_ylim()
plt.ylim(ymin, 0)
plt.xticks(xticks, xticks_labels, usetex=True)
yticks = np.concatenate(([0], H_args))
yticks_labels = ["${:.4f}$".format(yt) if yt != 0 else "${:.0f}$".format(yt) for yt in yticks]
plt.yticks(yticks, yticks_labels)
plt.xlabel(r"$\omega\;\rm(rad)$", fontsize=fs)
plt.ylabel(r"$\rm Fase\;(rad)$", fontsize=fs)
# se alinea la etiqueta del eje y de ambas graficas
ax.yaxis.set_label_coords(ylabel_xcoord, 0.5)
# se cambia el color de la etiqueta theta_0
# fase en thetas
plt.plot(np.array([thetas, thetas]), np.array([ymin * np.ones(H_args.shape), H_args]), 'k', lw=1)
plt.plot(np.array([np.zeros(H_args.shape), thetas]), np.array([H_args, H_args]), 'k', lw=1)
# etiqueta de la gráfica
plt.text(0.98, 0.85, r'$\angle H(e^{j\omega})$', fontsize=12, ha='right', va='baseline', transform = ax.transAxes)
plt.savefig('seq_and_sys_frequency_response_output_response.pdf', bbox_inches='tight')


# GRÁFICAS DE RESPUESTA EN FASE
# Retardo de fase
tauPhi_half = -phiH_half / w_half
# Retardo de fase en el vector thetas
tauPhis = -H_args / thetas

fig = plt.figure(1, figsize=(8, 5), frameon=False)
ax = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)
plt.plot(w_half, phiH_half, lw=2, color='k')
plt.xlim(xmin, xmax)
ymin, ymax = ax.get_ylim()
plt.ylim(ymin, 0)
plt.xticks(xticks, xticks_labels, usetex=True)
yticks = np.concatenate(([0], H_args))
yticks_labels = ["${:.4f}$".format(yt) if yt != 0 else "${:.0f}$".format(yt) for yt in yticks]
plt.yticks(yticks, yticks_labels)
plt.ylabel(r"$\rm Fase\;(rad)$", fontsize=fs)
# se alinea la etiqueta del eje y de ambas graficas
ax.yaxis.set_label_coords(ylabel_xcoord, 0.5)
# fase en thetas
plt.plot(np.array([thetas, thetas]), np.array([ymin * np.ones(H_args.shape), H_args]), 'k', lw=1)
plt.plot(np.array([np.zeros(H_args.shape), thetas]), np.array([H_args, H_args]), 'k', lw=1)
# etiqueta de la gráfica
plt.text(0.98, 0.85, r'$\angle H(e^{j\omega})$', fontsize=12, ha='right', va='baseline', transform = ax.transAxes)

ax = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)
plt.plot(w_half, tauPhi_half, lw=2, color='k')
plt.xlim(xmin, xmax)
ymin, ymax = ax.get_ylim()
ymin = 0
plt.ylim(0, ymax)
plt.xticks(xticks, xticks_labels, usetex=True)
yticks = np.concatenate(([0], tauPhis))
yticks_labels = ["${:.4f}$".format(yt) if yt != 0 else "${:.0f}$".format(yt) for yt in yticks]
plt.yticks(yticks, yticks_labels)
plt.ylabel(r"$\rm Retardo\;de\;fase\;(muestras)$", fontsize=fs)
# se alinea la etiqueta del eje y de ambas graficas
ax.yaxis.set_label_coords(ylabel_xcoord, 0.5)
# retardo de fase en thetas
plt.plot(np.array([thetas, thetas]), np.array([ymin * np.ones(tauPhis.shape), tauPhis]), 'k', lw=1)
plt.plot(np.array([np.zeros(tauPhis.shape), thetas]), np.array([tauPhis, tauPhis]), 'k', lw=1)
plt.xlabel(r"$\omega\;\rm(rad)$", fontsize=fs)
# etiqueta de la gráfica
plt.text(0.98, 0.85, r'$-\dfrac{\angle H(e^{j\omega})}{\omega}$', fontsize=12, ha='right', va='center', transform = ax.transAxes)
plt.savefig('seq_and_sys_frequency_response_output_phase_response.pdf', bbox_inches='tight')

# FILTRADO DE SINUSOIDES
nx = 200
n = np.arange(nx)
# entradas
x0 = np.sin(thetas[0] * n)
x1 = np.sin(thetas[1] * n)
x2 = np.sin(thetas[2] * n)
# salidas
y0 = signal.lfilter(b, a, x0)
y1 = signal.lfilter(b, a, x1)
y2 = signal.lfilter(b, a, x2)

n1 = 120
di = 0.01
ni = np.arange(n1, nx-1, di)
# señales interpoladas para mejor visualización
xi0 = np.sin(thetas[0] * ni)
xi1 = np.sin(thetas[1] * ni)
xi2 = np.sin(thetas[2] * ni)
yi0 = Gs[0] * np.sin(thetas[0] * ni + H_args[0])
yi1 = Gs[1] * np.sin(thetas[1] * ni + H_args[1])
yi2 = Gs[2] * np.sin(thetas[2] * ni + H_args[2])

# GRAFICAS DE RESPUESTA EN MAGNITUD DE SINUSOIDES
ms = 3.5
xmin = n1
xmax = nx-1
ymin = -1.1
ymax = 1.1
fig = plt.figure(2, figsize=(8, 6), frameon=False)
ax = plt.subplot2grid((3, 1), (0, 0), rowspan=1, colspan=1)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.plot(n[n1:], x0[n1:], ls='', marker='s', color='k', markersize=ms, label='${\\rm entrada}$')
plt.plot(ni, xi0, 'k', lw=1)
plt.plot(n[n1:], y0[n1:], ls='', marker='s', color='r', markersize=ms, label='${\\rm salida}$')
plt.plot(ni, yi0, 'r', lw=1)
yticks = [-1, 0, Gs[0]]
yticks_labels = ["${:.4f}$".format(yt) if yt % 1 != 0 else "${:.0f}$".format(yt) for yt in yticks]
plt.yticks(yticks, yticks_labels)
plt.plot([xmin, xmax], [Gs[0], Gs[0]], 'k', lw=1, zorder=-1)
leg = plt.legend(loc=1, frameon=True, fontsize=fs, framealpha=1)
leg.get_frame().set_edgecolor('k')
leg.get_frame().set_linewidth(0.5)
plt.title('${{\\rm Entrada\;de\;frecuencia\;}}\omega={:.2f}{{\;\\rm radianes}}$'.format(thetas[0]), fontsize=fs)

ax = plt.subplot2grid((3, 1), (1, 0), rowspan=1, colspan=1)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.plot(n[n1:], x1[n1:], ls='', marker='s', color='k', markersize=ms)
plt.plot(ni, xi1, 'k', lw=1)
plt.plot(n[n1:], y1[n1:], ls='', marker='s', color='r', markersize=ms)
plt.plot(ni, yi1, 'r', lw=1)
yticks = [-1, 0, Gs[1], 1]
yticks_labels = ["${:.4f}$".format(yt) if yt % 1 != 0 else "${:.0f}$".format(yt) for yt in yticks]
plt.yticks(yticks, yticks_labels)
plt.plot([xmin, xmax], [Gs[1], Gs[1]], 'k', lw=1, zorder=-1)
plt.title('${{\\rm Entrada\;de\;frecuencia\;}}\omega={:.2f}{{\;\\rm radianes}}$'.format(thetas[1]), fontsize=fs)

ax = plt.subplot2grid((3, 1), (2, 0), rowspan=1, colspan=1)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.plot(n[n1:], x2[n1:], ls='', marker='s', color='k', markersize=ms)
plt.plot(ni, xi2, 'k', lw=1)
plt.plot(n[n1:], y2[n1:], ls='', marker='s', color='r', markersize=ms)
plt.plot(ni, yi2, 'r', lw=1)
yticks = [-1, 0, Gs[2], 1]
yticks_labels = ["${:.4f}$".format(yt) if yt % 1 != 0 else "${:.0f}$".format(yt) for yt in yticks]
yticks_labels[1] = ''
plt.yticks(yticks, yticks_labels)
plt.plot([xmin, xmax], [Gs[2], Gs[2]], 'k', lw=1, zorder=-1)
plt.title('${{\\rm Entrada\;de\;frecuencia\;}}\omega={:.2f}{{\;\\rm radianes}}$'.format(thetas[2]), fontsize=fs)
plt.xlabel('$n$', fontsize=fs)

plt.tight_layout()
plt.savefig('seq_and_sys_frequency_response_output_response_gain.pdf', bbox_inches='tight')

# GRAFICAS DE RESPUESTA EN FASE DE SINUSOIDES
nx = 160
n = np.arange(nx)
di = 0.01
ni = np.arange(n1, nx-1, di)
# señales interpoladas para mejor visualización
xi0 = np.sin(thetas[0] * ni)
xi1 = np.sin(thetas[1] * ni)
xi2 = np.sin(thetas[2] * ni)
yi0 = Gs[0] * np.sin(thetas[0] * ni + H_args[0])
yi1 = Gs[1] * np.sin(thetas[1] * ni + H_args[1])
yi2 = Gs[2] * np.sin(thetas[2] * ni + H_args[2])

# máximo de la sinusoide de entrada mayor a n1
thetas = np.array(thetas)
kmax = np.ceil(thetas * n1 / (2 * np.pi) - 1 / 4)
nmax = np.pi * (1 + 4 * kmax) / (2 * thetas)

# parametros de la gráfica
xmin = n1
xmax = nx-1
ymin = -1.1
ymax = 2
# parámetros de la etiqueta de la cantidad de muestras
h = 1.2
ht = 1.4
bt = 1.35

fig = plt.figure(3, figsize=(8, 6), frameon=False)
ax = plt.subplot2grid((3, 1), (0, 0), rowspan=1, colspan=1)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.plot(n[n1:nx], x0[n1:nx], ls='', marker='s', color='k', markersize=ms, label='${\\rm entrada}$')
plt.plot(ni, xi0, 'k', lw=1)
plt.plot(n[n1:nx], y0[n1:nx], ls='', marker='s', color='r', markersize=ms, label='${\\rm salida}$')
plt.plot(ni, yi0, 'r', lw=1)
# dibujo de la etuqueta
i = 0
plt.plot([nmax[i], nmax[i] + tauPhis[i]], [h, h], 'k', lw=1)
plt.plot([nmax[i], nmax[i]], [1, ht], 'k', lw=1)
plt.plot([nmax[i] + tauPhis[i], nmax[i] + tauPhis[i]], [Gs[i], ht], 'k', lw=1, zorder=-1)
plt.text(nmax[i] + tauPhis[i] / 2, bt, '${:.4f}{{\;\\rm muestras}}$'.format(tauPhis[i]),
         fontsize=fs, ha='center', va='baseline')
leg = plt.legend(loc=1, frameon=True, fontsize=fs, framealpha=1)
leg.get_frame().set_edgecolor('k')
leg.get_frame().set_linewidth(0.5)
plt.title('${{\\rm Entrada\;de\;frecuencia\;}}\omega={:.2f}{{\;\\rm radianes}}$'.format(thetas[0]), fontsize=fs)

ax = plt.subplot2grid((3, 1), (1, 0), rowspan=1, colspan=1)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.plot(n[n1:nx], x1[n1:nx], ls='', marker='s', color='k', markersize=ms)
plt.plot(ni, xi1, 'k', lw=1)
plt.plot(n[n1:nx], y1[n1:nx], ls='', marker='s', color='r', markersize=ms)
plt.plot(ni, yi1, 'r', lw=1)
# dibujo de la etuqueta
i = 1
plt.plot([nmax[i], nmax[i] + tauPhis[i]], [h, h], 'k', lw=1)
plt.plot([nmax[i], nmax[i]], [1, ht], 'k', lw=1)
plt.plot([nmax[i] + tauPhis[i], nmax[i] + tauPhis[i]], [Gs[i], ht], 'k', lw=1, zorder=-1)
plt.text(nmax[i] + tauPhis[i] / 2, bt, '${:.4f}{{\;\\rm muestras}}$'.format(tauPhis[i]),
         fontsize=fs, ha='center', va='baseline')
plt.title('${{\\rm Entrada\;de\;frecuencia\;}}\omega={:.2f}{{\;\\rm radianes}}$'.format(thetas[1]), fontsize=fs)

ax = plt.subplot2grid((3, 1), (2, 0), rowspan=1, colspan=1)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.plot(n[n1:nx], x2[n1:nx], ls='', marker='s', color='k', markersize=ms)
plt.plot(ni, xi2, 'k', lw=1)
plt.plot(n[n1:nx], y2[n1:nx], ls='', marker='s', color='r', markersize=ms)
plt.plot(ni, yi2, 'r', lw=1)
# dibujo de la etuqueta
i = 2
plt.plot([nmax[i], nmax[i] + tauPhis[i]], [h, h], 'k', lw=1)
plt.plot([nmax[i], nmax[i]], [1, ht], 'k', lw=1)
plt.plot([nmax[i] + tauPhis[i], nmax[i] + tauPhis[i]], [Gs[i], ht], 'k', lw=1, zorder=-1)
plt.text(nmax[i] + tauPhis[i] / 2, bt, '${:.4f}{{\;\\rm muestras}}$'.format(tauPhis[i]),
         fontsize=fs, ha='center', va='baseline')
plt.title('${{\\rm Entrada\;de\;frecuencia\;}}\omega={:.2f}{{\;\\rm radianes}}$'.format(thetas[2]), fontsize=fs)
plt.xlabel('$n$', fontsize=fs)

plt.tight_layout()
plt.savefig('seq_and_sys_frequency_response_output_response_phase.pdf', bbox_inches='tight')


plt.show()
