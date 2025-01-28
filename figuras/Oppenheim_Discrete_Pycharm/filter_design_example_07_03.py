import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from matplotlib.patches import Polygon

__author__ = 'ernesto'

# use latex
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

# Diseño de filtro Butterworth pasabajos mediante la transformación bilineal
# Ejemplo 7.3 del Oppenheim

def butterworth_filter_parameters(x):
    return [1 + (wp / (x[0] * Td)) ** (2 * x[1]) - (1 / gp) ** 2,
            1 + (ws / (x[0] * Td)) ** (2 * x[1]) - (1 / gs) ** 2]

# Parametros
# Especificaciones del filtro en tiempo
# Ganancia mínima en banda pasante (dB) y frecuencia (rad)
gp_db = -1
wp = 0.2 * np.pi
# Ganancia máxima en banda suprimida (dB) y frecuencia (rad)
gs_db = -15
ws = 0.3 * np.pi
# Td: período de muestreo en invarinza al impulso
Td = 1

# Conversión de decibeles a escala lineal
gp = 10 ** (gp_db / 20)
gs = 10 ** (gs_db / 20)
print('Ganancia mínima en banda pasante: {:.5f}, Antenuación mínima en banda suprimida: {:.5f}'.format(gp, gs))
# Conversion a las frecuencias en tiempo continuo
Wp = (2 / Td) * np.tan(wp / 2)
Ws = (2 / Td) * np.tan(ws / 2)
# Calculo de del orden N y frecuencia de corte W_c del filtro en tiempo continuo
N = np.log10((1 - (1 / gs) ** 2) / (1 - (1 / gp) ** 2)) / (2 * np.log10(Ws / Wp))
Wc = Wp * ((1 / ((1 / gp) ** 2 - 1)) ** (1 / (2 * N)))
print('N: {:.4f}'.format(N))
print('Wc: {:.5f}'.format(Wc))
# N debe ser el entero superior
N = int(np.ceil(N))
# Se calcula W_c para que las especificaciones se cumplan exactas en la banda suprimida.
Wc = Ws / ((1 / gs) ** 2 - 1) ** (1 / (2 * N))
print('N entero: {:.4f}'.format(N))
print('Wc: {:.5f}'.format(Wc))
# Cálculo de los polos del sistema en el semiplano derecho
ks = np.arange(N)
sk = np.exp(1j * np.pi * (2 * ks + N + 1) / (2 * N)) * Wc
np.set_printoptions(precision=4, suppress=True, floatmode='fixed')
print('Polos del sistema en tiempo continuo: {}'.format(sk))
# Numerador y denominador
Hs_num = Wc ** N
print('Numerador: {:.5f}'.format(Hs_num))
print('Factores de segundo grado del denominador: ')
# Hay que distinguir los casos en que N es par o impar.
if N % 2 == 0:
    for i in np.arange(N // 2):
        print(np.real(np.polymul([1, -sk[i]], [1, -np.conj(sk[i])])))
else:
    for i in np.arange((N - 1) // 2):
        print(np.real(np.polymul([1, -sk[i]], [1, -np.conj(sk[i])])))
    print(np.real([1, sk[(N + 1) / 2]]))
# Aplicación de la transformación bilineal
num = Hs_num
zk = np.zeros(sk.shape, dtype=complex)
for i in np.arange(N):
    num /= (2 / Td - sk[i])
    zk[i] = (2 / Td + sk[i]) / (2 / Td - sk[i])
print('Numerador de H(z): {}'.format(num))
print('Factor multiplicativo de z^{{-1}} del denominador de H(z) (1-z^{{-1}}z_k): {}'.format(zk))
# Agrupación de términos complejos conjugados
print('Factor del numerador de H(z): {:.7f}'.format(np.real(num)))
print('Factores de segundo grado del denominador de H(z): ')
# Hay que distinguir los casos en que N es par o impar.
if N % 2 == 0:
    for i in np.arange(N // 2):
        print(np.real(np.polymul([1, -zk[i]], [1, -np.conj(zk[i])])))
else:
    for i in np.arange((N - 1) // 2):
        print(np.real(np.polymul([1, -zk[i]], [1, -np.conj(zk[i])])))
    print(np.real([1, zk[(N + 1) / 2]]))
# Desarrollo del numerador y denominador de H(z)
Hz_num = 1
Hz_den = 1
for i in np.arange(N):
    Hz_num = np.polymul(Hz_num, [1, 1])
    Hz_den = np.polymul(Hz_den, [1, -zk[i]])
np.set_printoptions(precision=6, suppress=True, floatmode='fixed')
print('Numerador de H(z): {}'.format(np.real(Hz_num)))
print('Denominador de H(z): {}'.format(np.real(Hz_den)))

# H(z) en forma directa
b_Hz = num * np.real(Hz_num)
a_Hz = np.real(Hz_den)
print('H(z) en forma directa. Numerador: {}'.format(np.real(b_Hz)))
print('H(z) en forma directa. Denominador: {}'.format(np.real(a_Hz)))

# Respuesta en frecuencia
w, Hz = signal.freqz(b_Hz, a_Hz, worN=512, whole=False, plot=None, fs=2*np.pi, include_nyquist=False)
_, Hz_grp = signal.group_delay((b_Hz, a_Hz), w=512, whole=False, fs=2*np.pi)
Hz_grp[-1] = 'nan'


# Graficas

### Parámetros de la gráfica
fontsize = 11
fontsize2 = 12
# x y ticks labels margin
xtm = -0.07
ytm = -0.04
display_length = 6

xmin = 0
xmax = np.pi
dx = 0.15
xmax_ax = xmax + dx
xmin_ax = xmin - dx
ymin_ax = -0.1
ymax_ax = 1.2

grey = [0.9, 0.9, 0.9]

fig = plt.figure(0, figsize=(9, 6), frameon=False)
ax = plt.subplot2grid((8, 4), (0, 0), rowspan=5, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
# magnitud de las respuestas en frecuencia
plt.plot(w, np.absolute(Hz), 'r-', lw=2)

# mascara
plt.plot([0, ws], [1, 1], 'k-')
plt.plot([0, wp], [gp, gp], 'k-')
plt.plot([ws, np.pi], [gs, gs], 'k-')
plt.plot([wp, wp], [0, gp], 'k-')
plt.plot([ws, ws], [gs, 1], 'k-')
# región pintada
vert = np.vstack(([0, wp, wp, 0], [0, 0, gp, gp]))
p1 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p1)
mask_max = 1.1
vert = np.vstack(([0, ws, ws, 0], [1, 1, mask_max, mask_max]))
p2 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p2)
vert = np.vstack(([ws, xmax, xmax, ws], [gs, gs, mask_max, mask_max]))
p3 = Polygon(vert.T, fc=grey, alpha=1, zorder=-1)
ax.add_artist(p3)

# etiquetas
# eje x
xtm2 = -0.11
plt.plot([wp, wp], [0, xtl], 'k-', lw=1)
plt.text(wp, xtm, '${:.1f}\pi$'.format(wp / np.pi), fontsize=fontsize, ha='center', va='baseline')
plt.plot([ws, ws], [0, xtl], 'k-', lw=1)
plt.text(ws, xtm, '${:.1f}\pi$'.format(ws / np.pi), fontsize=fontsize, ha='center', va='baseline')
plt.plot([np.pi, np.pi], [0, xtl], 'k-', lw=1)
plt.text(np.pi, xtm, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytm, xtm, '$0$', fontsize=fontsize, ha='right', va='baseline')
# eje y
plt.plot([0, xtl], [gp, gp], 'k-', lw=1)
plt.text(ytm, gp, '${:.4f}$'.format(gp), fontsize=fontsize, ha='right', va='center')
plt.plot([0, xtl], [1, 1], 'k-', lw=1)
plt.text(ytm, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, xtl], [gs, gs], 'k-', lw=1)
plt.text(ytm, gs, '${:.4f}$'.format(gs), fontsize=fontsize, ha='right', va='center')
# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.text(0.06, ymax_ax, '$|H(e^{j\omega})|$', fontsize=fontsize2, ha='left', va='center')
ylab = -0.2
plt.text(ylab, (ymax_ax+ymin_ax)/2, r'$\textrm{Amplitud}$', fontsize=fontsize, ha='center', va='center', rotation=90)

plt.axis('off')

# Retardo de grupo
# x y ticks labels margin
xtm = -1.2
ytm = -0.04

xmin = 0
xmax = np.pi
dx = 0.15
xmax_ax = xmax + dx
xmin_ax = xmin - dx
ymin_ax = -1
ymax_ax = 12

xticks = np.linspace(0.2, 1, 5) * np.pi
xticks_labels = ['${:.1f}\pi$'.format(t) for t in np.linspace(0.2, 1, 5)]
xticks_labels[-1] = '$\pi$'
yticks = np.linspace(2, 10, 5)
yticks_labels = ['${}$'.format(int(t)) for t in np.linspace(2, 10, 5)]

ax = plt.subplot2grid((8, 4), (5, 0), rowspan=3, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))

plt.plot(w, Hz_grp, 'k-', lw=2)

for i, xt in enumerate(xticks):
    plt.plot([xt, xt], [0, xtl], 'k-', lw=1)
    plt.text(xt, xtm, xticks_labels[i], fontsize=fontsize, ha='center', va='baseline')

for i, yt in enumerate(yticks):
    plt.plot([0, ytl], [yt, yt], 'k-', lw=1)
    plt.text(ytm, yt, yticks_labels[i], fontsize=fontsize, ha='right', va='center')

plt.text(ytm, xtm, '$0$', fontsize=fontsize, ha='right', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.text(0.06, ymax_ax, r'$\textrm{grd}[H(e^{j\omega})]$', fontsize=fontsize2, ha='left', va='center')

plt.text(ylab, (ymax_ax+ymin_ax)/2, r'$\textrm{Muestras}$', fontsize=fontsize, ha='center', va='center', rotation=90)

plt.axis('off')

plt.savefig('filter_design_example_07_03_freq_response.pdf', bbox_inches='tight')

#
# Diagrama  de polos yceros
#
# Polos de H(s)H(-s)
ks = np.arange(2 * N)
sk = np.exp(1j * np.pi * (2 * ks + N + 1) / (2 * N)) * Wc

# círculo donde están los polos
M = 200
thetas = np.linspace(0, 2 * np.pi, M)
xz = Wc * np.cos(thetas)
yz = Wc * np.sin(thetas)

# valores maximos y minimos de los ejes
max_ax = 1.1
xmin = -max_ax
xmax = max_ax
ymin = -max_ax
ymax = max_ax
# axis parameters
xmin_ax = xmin
xmax_ax = xmax
ymin_ax = ymin
ymax_ax = ymax

# x ticks labels margin
xtm = -0.18
ytm = -0.09
# font size
fontsize = 14
fontsize2 = 10

fig = plt.figure(1, figsize=(4, 4), frameon=False)
ax = fig.add_subplot(111)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.gca().set_aspect('equal', adjustable='box')

# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

# ceros
zeros_marker_style = dict(marker='o', linestyle='', markersize=8, markerfacecolor='w', markeredgecolor='k',
                          markeredgewidth=1.5)
# polos
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgecolor='k', markeredgewidth=1.5)
plt.plot(sk.real, sk.imag, **polos_marker_style)

# angulo
k = 9
plt.plot([0, sk[k].real], [0, sk[k].imag], 'k--', lw=1)
plt.plot([0, sk[k + 1].real], [0, sk[k + 1].imag], 'k--', lw=1)
r = 0.5
nn = 20
ths = np.linspace(np.angle(sk[k]), np.angle(sk[k + 1]), nn)
xa = r * np.cos(ths)
ya = r * np.sin(ths)
plt.plot(xa, ya, 'k-', lw=1)
plt.annotate('$\dfrac{{\pi}}{{{}}}$'.format(N), xytext=(xmax_ax, 0.6), xycoords='data', xy=(xa[int(nn/2)], ya[int(nn/2)]),
             textcoords='data', fontsize=fontsize, va="center", ha="right",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0, 0.4),
                             patchA=None, patchB=None, shrinkA=4, shrinkB=0,
                             connectionstyle="angle3,angleA={:d},angleB={:d}".format(180, 30)))

# longitud
k = 3
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(sk[k].real, sk[k].imag), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=3, headlength=8, facecolor='black', shrink=0.002))
plt.annotate('$\Omega_c={:.5f}$'.format(Wc), xytext=(xmin_ax, -0.4), xycoords='data',
             xy=(sk[k].real / 2, sk[k].imag / 2),
             textcoords='data', fontsize=fontsize, va="center", ha="left",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.8, 1),
                             patchA=None, patchB=None, shrinkA=4, shrinkB=0,
                             connectionstyle="angle3,angleA={:d},angleB={:d}".format(20, 70)))

# etiquetas de los ejes
plt.text(xmax_ax, xtm, r'$\textrm{Re}(s)$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, r'$\textrm{Im}(s)$', fontsize=fontsize, ha='left', va='center')

plt.axis('off')

# save as pdf image
plt.savefig('filter_design_example_07_03_poles.pdf', bbox_inches='tight')

plt.show()
