import matplotlib.pyplot as plt
import numpy as np
from scipy import signal, fft

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

#
# parámetros del filtro pasabajos
#
# Frecuencias de los bordes de las bandas pasante y suprimida
wp = 0.32 * np.pi
ws = 0.40 * np.pi
# ripple de las bandas pasante y suprimida
delta1 = 0.1
delta2 = 0.015

# taplength: largo del filtro: M = taplength - 1, L = M / 2
# debe ser un número impar
taplength = 39

#
# parámetros del algoritmo
#
# densidad de la grilla de frecuencias de Ae(e^{jw}): nomegas = griddensity * nalternations
griddensity = 20
# error relativo mínimo para terminar el algoritmo
min_error = 1e-7
# número máximo de itraciones del algoritmo
maxiterationcount = 20

# nalternations = L + 2
nalternations = ((taplength - 1) // 2) + 2
# listas para almacenar el resultado de cada iteración
alternations_freq_idx_list = []
new_alternations_freq_idx_list = []
error_list = []
amplitude_list = []
delta_list = []

# construcción de la grilla densa de frecuencias equiespaciadas
# cantidad de slots de frecuencia
nomega = griddensity * nalternations
# indices del vector de frecuencias: omega_idx = 0: nomega (largo = nomega + 1)
omega_idx = np.arange(nomega + 1)
# vector de frecuencias en radianes
omegas = omega_idx * np.pi / nomega

# construcción de la respuesta en frecuencia deseada en la grilla de frecuencias (ecuación 7.97)
Hd = np.zeros(nomega + 1)
Hd[omegas <= (wp + ws) / 2] = 1

# construcción del vector de pesos en la grilla de frecuencia (ecuación 7.98)
W = np.zeros(nomega + 1)
K = delta1 / delta2
W[omegas <= wp] = 1 / K
W[omegas >= ws] = 1

# indices de las frecuencias de los bordes de las bandas pasante y suprimida
# indices correspondientes a la banda de transición
trans_band_freq_idx = np.where(W == 0)[0]
# be1, be2: índice de la última frecuencia de la banda pasante y primer frecuencia de la banda atenuada
wp_idx = trans_band_freq_idx[0] - 1
ws_idx = trans_band_freq_idx[-1] + 1
# índices de las frecuencias de la banda pasante y suprimida
non_trans_band_freq_idx = np.concatenate([np.arange(wp_idx + 1), np.arange(ws_idx, nomega + 1)])

# selección inicial de los índices de frecuencia de los puntos de alternancia (extremos)
# se eligen elementos de non_trans_band_freq_idx equiespaciados
alternations_idx = np.linspace(0, len(non_trans_band_freq_idx) - 1, nalternations)
alternations_idx = np.round(alternations_idx)
alternations_idx = alternations_idx.astype(int)
freq_idx_alternations = np.take(non_trans_band_freq_idx, alternations_idx)

sign_array = (-1) ** np.arange(nalternations)

# respuesta deseada en las alternancias
Hd_alternations = np.array(nalternations)
# respuesta en frecuencia en la grilla densa de frecuencias
Ae = np.zeros(nomega + 1)
iterationcount = 0
# diferencia relativa entre el máximo del error y delta para la decisión de
# terminación del algoritmo
relativedifference = np.zeros(maxiterationcount)
# indices de frecuencia de las nuevas alternancias
new_freq_idx_alternations = np.empty(nalternations, dtype=np.int32)

while 1:
    iterationcount += 1
    # cos(omega_k), donde omega_k es la frecuencia de los extremos, k = 1: numberofextrema
    cos_freq_aternations = np.cos(freq_idx_alternations * np.pi / nomega)
    # resuesta deseada y ponderación en las alternancias
    Hd_alternations = np.take(Hd, freq_idx_alternations)
    W_alternations = np.take(W, freq_idx_alternations)
    # x = cos(omega) en la grilla densa de frecuencias
    cos_omega = np.cos(omegas)
    # inicialización de los coeficientes b_k de Lagrange
    bk = np.zeros(nalternations)

    # cálculo de los coeficientes bk (ecuación 7.115)
    for i in np.arange(nalternations):
        bk[i] = 1
        for j in np.arange(nalternations):
            if j != i:
                bk[i] = bk[i] * (cos_freq_aternations[i] - cos_freq_aternations[j])
        bk[i] = 1 / bk[i]

    # cálculo de delta (ecuación 7.114)
    delta = np.sum(bk * Hd_alternations) / np.sum(bk * sign_array / W_alternations)
    # lo que sigue es para el cálculo de la amplitud de los extremos (ecuación 7.116a)
    # cálculo de C_k (ecuación 7.116b)
    Ck = Hd_alternations - delta * sign_array / W_alternations
    # cálculo de d_k (ecuación 7.116c)
    dk = bk * (cos_freq_aternations - cos_freq_aternations[-1])

    # cálculo de A(omega) (ecuación 7.116a)
    for i in omega_idx:
        # hay que distinguir el caso en que en que la frecuencia de la grilla densa
        # es una alternancia o no
        if not (np.isin(i, freq_idx_alternations)):
            Ae_num = np.sum(dk[:-1] * Ck[:-1] / (cos_omega[i] - cos_freq_aternations[:-1]))
            Ae_den = np.sum(dk[:-1] / (cos_omega[i] - cos_freq_aternations[:-1]))
            Ae[i] = Ae_num / Ae_den
        else:
            j = np.nonzero(freq_idx_alternations == i)[0]
            Ae[i] = Ck[j]
    # cálculo del error de aproximación (ecuación 7.96)
    E = W * (Hd - Ae)

    # búsqueda de las nuevas alternancias
    j = 0
    for i in np.arange(1, nomega):
        if abs(E[i]) > abs(E[i - 1]) and abs(E[i]) > abs(E[i + 1]):
            new_freq_idx_alternations[j] = i
            j = j + 1
    # se agregan los extremos de las bandas de frecuencias si no están
    if wp_idx not in new_freq_idx_alternations:
        new_freq_idx_alternations[j] = wp_idx
        j = j + 1
    if ws_idx not in new_freq_idx_alternations:
        new_freq_idx_alternations[j] = ws_idx
        j = j + 1
    # se agrega el extremo con mayor error o ambos si el número de
    # alternancias no es nalternations
    if abs(E[0]) > abs(E[-1]):
        new_freq_idx_alternations[j] = 0
        j = j + 1
        if j < nalternations:
            new_freq_idx_alternations[j] = nomega
    else:
        new_freq_idx_alternations[j] = nomega
        j = j + 1
        if j < nalternations:
            new_freq_idx_alternations[j] = 0

    # se almacenan los resultados para graficar
    new_freq_idx_alternations = np.sort(new_freq_idx_alternations)
    relativedifference[iterationcount] = (max(abs(E)) - abs(delta)) / abs(delta)
    alternations_freq_idx_list.append(freq_idx_alternations.copy())
    new_alternations_freq_idx_list.append(new_freq_idx_alternations.copy())
    error_list.append(E)
    amplitude_list.append(Ae.copy())
    delta_list.append(delta)
    if relativedifference[iterationcount] < min_error:
        print("I have converged at iteration: ", iterationcount)
        flag = 0
        break
    freq_idx_alternations = new_freq_idx_alternations.copy()
    if iterationcount == maxiterationcount:
        break

# cálculo de la respuesta al impulso h mediante la DFT inversa de Ae
HH = np.concatenate((Ae, np.flip(Ae[1:])))
h = fft.ifft(HH)
M = taplength - 1
L = M // 2
h = np.real(h[: L + 1])
h = np.concatenate((np.flip(h[1:]), h))

### comparación con el resultado de la función remez de scipy
hsci = signal.remez(taplength, [0, wp, ws, np.pi], [1, 0], [1 / K, 1], fs=2*np.pi)
print('sum_n|h[n] - h_remez[n]| = {}'.format(np.sum(np.abs(h - hsci))))

#
# Gráficas
#
# Error de aproximación y alternancias en cada iteración

fs = 12
subplotheight = 2

xticks = np.linspace(0, 1, 11)
xticks_labels = ['${:.1f}\pi$'.format(xt) for xt in xticks]
xticks_labels[0] = '$0$'
xticks_labels[-1] = '$\pi$'
xticks *= np.pi

new_alt_marker_style = dict(marker='o', linestyle='', markersize=5, markerfacecolor='w',
                            markeredgewidth=1.5, markeredgecolor='r')
alt_marker_style = dict(marker='s', linestyle='', markersize=5, markerfacecolor='w',
                        markeredgewidth=1.5, markeredgecolor='b')

fig = plt.figure(0, figsize=(8, 7), frameon=False)
for k in np.arange(iterationcount):
    ax = plt.subplot2grid((subplotheight * iterationcount, 1), (subplotheight * k, 0), rowspan=subplotheight)
    plt.xlim(omegas[0], omegas[-1])
    Ek = error_list[k]
    Ek[trans_band_freq_idx] = 'nan'
    alt_freq_idx_k = alternations_freq_idx_list[k]
    new_alt_freq_idx_k = new_alternations_freq_idx_list[k]
    plt.plot(omegas, Ek, 'k-', lw=1.5, label='$E(\omega)$')
    plt.plot(omegas[alt_freq_idx_k], Ek[alt_freq_idx_k], **alt_marker_style, zorder=10,
             label=r'$\textrm{Alternancias}$')
    plt.plot(omegas[new_alt_freq_idx_k], Ek[new_alt_freq_idx_k], **new_alt_marker_style, zorder=15,
             label=r'$\textrm{Nuevas alternancias}$')
    if k == 0:
        plt.legend(loc='upper right', fontsize=fs, frameon=False, framealpha=1)
    if k != iterationcount - 1:
        plt.xticks(xticks, [])
    else:
        plt.xlabel(r'$\omega$', fontsize=fs)
        plt.xticks(xticks, xticks_labels, usetex=True)

plt.savefig('filter_design_fir_optimal_pm_implementation_plots_error.pdf', bbox_inches='tight')


# Respuesta en frecuencia y respuesta al impulso

fig = plt.figure(1, figsize=(8, 5), frameon=False)
# Respuesta en frecuencia
ax = plt.subplot2grid((11, 1), (0, 0), rowspan=6, colspan=1)
plt.xlim(omegas[0], omegas[-1])
plt.plot(omegas, Ae, 'k-', lw=1.5)
plt.text(0.99, 0.9, '$A_e(e^{j\omega})$', fontsize=fs, ha='right', va='baseline', transform=ax.transAxes)
plt.xlabel(r'$\omega$', fontsize=fs)
plt.xticks(xticks, xticks_labels, usetex=True)

# Respuesta al impulso
n = np.arange(taplength)
ax = plt.subplot2grid((11, 1), (7, 0), rowspan=5, colspan=1)
plt.xlim(n[0], n[-1])
(markers, stemlines, bl) = plt.stem(n, h, linefmt='k', markerfmt='s', use_line_collection=True)
plt.setp(markers, markersize=4, markeredgecolor='k', markerfacecolor='k')
plt.setp(bl, visible=False)
plt.plot([n[0], n[-1]], [0, 0], 'k-', lw=1, zorder=-1)
plt.text(0.99, 0.85, '$h[n]$', fontsize=fs, ha='right', va='baseline', transform=ax.transAxes)
plt.xlabel('$n$')

plt.savefig('filter_design_fir_optimal_pm_implementation_plots_Ae_h.pdf', bbox_inches='tight')

# Amplitud en la iteración k
k = 1
alt_freq_idx_k = alternations_freq_idx_list[k]
new_alt_freq_idx_k = new_alternations_freq_idx_list[k]
Ak = amplitude_list[k]
deltak = delta_list[k]

# banda pasante
xmin = 0
dw = 0.1
dy = 0.2
xmax = wp + dw
xmax_ax = xmax
xmin_ax = xmin
ymax_ax = 1 + dy
ymin_ax = 1 - dy

fig = plt.figure(2, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=2)

plt.plot(omegas, Ak, 'k-', lw=1.5)
plt.plot(omegas[alt_freq_idx_k], Ak[alt_freq_idx_k], **alt_marker_style, zorder=10,
             label=r'$\textrm{Alternancias previas}$')
plt.plot(omegas[new_alt_freq_idx_k], Ak[new_alt_freq_idx_k], **new_alt_marker_style, zorder=15,
             label=r'$\textrm{Nuevas alternancias}$')
plt.plot([xmin_ax, wp], [1 + deltak * K, 1 + deltak * K], 'k--', lw=1)
plt.plot([xmin_ax, wp], [1 - deltak * K, 1 - deltak * K], 'k--', lw=1)
plt.xlabel(r'$\omega$', fontsize=fs)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.title(r'$\textrm{Banda pasante de }A_e(e^{j\omega})$', fontsize=fs)

# banda atenuada

xmin = ws - dw
xmax = np.pi
xmax_ax = xmax
xmin_ax = xmin
ymax_ax = dy / K
ymin_ax = -dy / K

ax = plt.subplot2grid((4, 4), (0, 2), rowspan=4, colspan=2)
plt.plot(omegas, Ak, 'k-', lw=1.5)
plt.plot(omegas[alt_freq_idx_k], Ak[alt_freq_idx_k], **alt_marker_style, zorder=10,
             label=r'$\textrm{Alternancias}$')
plt.plot(omegas[new_alt_freq_idx_k], Ak[new_alt_freq_idx_k], **new_alt_marker_style, zorder=15,
             label=r'$\textrm{Nuevas alternancias}$')
plt.plot([ws, xmax_ax], [deltak, deltak], 'k--', lw=1)
plt.plot([ws, xmax_ax], [-deltak,-deltak], 'k--', lw=1)

ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
plt.xlabel(r'$\omega$', fontsize=fs)
plt.xticks(xticks, xticks_labels, usetex=True)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.title(r'$\textrm{Banda suprimida de }A_e(e^{j\omega})$', fontsize=fs)
plt.legend(loc='upper right', fontsize=fs, frameon=False, framealpha=1)

plt.savefig('filter_design_fir_optimal_pm_implementation_plots_Ae_iter_2.pdf', bbox_inches='tight')

plt.show()
