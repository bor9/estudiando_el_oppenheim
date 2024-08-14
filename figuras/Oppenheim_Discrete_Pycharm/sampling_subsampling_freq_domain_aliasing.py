import matplotlib.pyplot as plt
import numpy as np
import math

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


def create_signal_spectrum(NomegaN):
    X = np.sqrt(np.hanning(2 * NomegaN + 1))
    Lh = round(NomegaN / 2)
    X[NomegaN - Lh: NomegaN + Lh + 1] = X[NomegaN - Lh: NomegaN + Lh + 1] + 1.5 * np.hanning(2 * Lh + 1) ** 2
    return X / np.amax(X)


def create_impulse_train_spectrum(Nomega, X, NomegaN, Nshift=0):
    # relleno de ceros para simplificar la construcción
    zpad = Nomega + 2 * NomegaN + 1
    Xs = np.zeros((Nfrec + 2 * zpad,))
    if0_aux = if0 + zpad
    # número de copias del espectro
    nc = math.ceil(((Nfrec - 1) / 2) / Nomega)
    for i in np.arange(-nc, nc + 1):
        # indice de la frecuencia central
        ic = if0_aux + i * Nomega + Nshift
        Xs[ic - NomegaN: ic + NomegaN + 1] += X
    # se elimina el relleno de ceros
    return Xs[zpad: Nfrec + zpad]

#
# Parámetros
#
# Períodos de muestreo
T = 0.05
# Frecuencia de muestreo
omega1 = 2 * math.pi / T
# rango de frecuencias del eje horizontal
fmin = -2.6 * omega1
fmax = 2.6 * omega1
# frecuencia máxima del espectro de la señal continua
omegaN = 0.25 * omega1
# cantidad de frecuencias en el eje horizontal
Nfrec = 2001
omega = np.linspace(fmin, fmax, Nfrec, endpoint=True)
# índice del elemento de frecuencia 0
if0 = int((Nfrec - 1) / 2)
# intervalo de frecuencia entre dos muestras en el eje horizontal
Domega = (fmax - fmin) / (Nfrec - 1)
# se redefine omega1, omega2 y omegaN para que sean
# un número entero de muestras de frecuencia
Nomega1 = 2 * round(omega1 / (2 * Domega))
omega1 = Nomega1 * Domega
NomegaN = round(omegaN / Domega)
omegaN = NomegaN * Domega

# Construcción del espectro
X = create_signal_spectrum(NomegaN)

# Construcción de Xc
Xc = np.zeros(omega.shape)
Xc[if0 - NomegaN: if0 + NomegaN + 1] = X

# Construcción de Xs
Xs = create_impulse_train_spectrum(Nomega1, X, NomegaN)

# Construcción del espectro de la señal submuestreada
M = 3
X = create_signal_spectrum(M * NomegaN)
Xsub1 = create_impulse_train_spectrum(M * Nomega1, X, M * NomegaN)
Xsub2 = create_impulse_train_spectrum(M * Nomega1, X, M * NomegaN, Nomega1)
Xsub3 = create_impulse_train_spectrum(M * Nomega1, X, M * NomegaN, -Nomega1)
# señal submuestreada
Xsub = create_impulse_train_spectrum(Nomega1, X, M * NomegaN)

# para las etiquetas de los ejes
nc = math.floor(((Nfrec - 1) / 2) / (Nomega1 / 2))
# altura. esto es 1/T
d_heigh = 1.3

# Gráficas

df = 15
xmin_ax = fmin - df
xmax_ax = fmax + df
ymin_ax =-0.1
ymax_ax = 1.8

display_length = 6
fontsize = 10
fontsize2 = 12
# x y ticks labels margin
xtm = -0.27
xtm2 = -0.35
ytm = 6
ytm2 = -6
grey = [0.7, 0.7, 0.7]
dymax = 0.15

fig = plt.figure(0, figsize=(9.5, 11 * 5 / 6), frameon=False)

### X(e^{jw}) con T1
ax = plt.subplot2grid((5, 3), (0, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax - dymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Xs * d_heigh, 'k-', lw=1.5)
# etiquetas de los ejes
for i in np.arange(-nc, nc + 1):
    omega_i = i * omega1 / 2
    plt.plot([omega_i, omega_i], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(omega_i, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(omega_i, xtm, '${}{}\pi$'.format(sg, np.absolute(i)),
                 fontsize=fontsize, ha='center', va='baseline')
plt.text(-omegaN, xtm, '$-\omega_N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(omegaN, xtm, '$\omega_N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, ytl], [d_heigh, d_heigh], 'k-', lw=1)
plt.text(ytm, ymax_ax - dymax, '$X(e^{j\omega})$', fontsize=fontsize2, ha='left', va='center')
plt.annotate('$\dfrac{\pi}{2}$', xytext=(2.1 * omegaN, 0.5), xycoords='data',
             xy=(omegaN, 0), textcoords='data', fontsize=fontsize, va="center", ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0, 0.5),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=0))
plt.text(xmin_ax + 4, ymax_ax - 0.2, '$\omega_N=\dfrac{\pi}{2}$', fontsize=fontsize2, ha='left', va='top')
plt.axis('off')

### X(e^{jw/2})
ax = plt.subplot2grid((5, 3), (1, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax - dymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Xsub1 * d_heigh, 'k-', lw=1.5)
# etiquetas de los ejes
for i in np.arange(-nc, nc + 1):
    omega_i = i * omega1 / 2
    plt.plot([omega_i, omega_i], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(omega_i, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(omega_i, xtm, '${}{}\pi$'.format(sg, np.absolute(i)),
                 fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh + 0.1, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, ytl], [d_heigh, d_heigh], 'k-', lw=1)
plt.text(ytm, ymax_ax - dymax, '$X(e^{j\omega/3})$', fontsize=fontsize2, ha='left', va='center')
plt.annotate('$\dfrac{3\pi}{2}$', xytext=(4.1 * omegaN, 0.5), xycoords='data',
             xy=(M * omegaN, 0), textcoords='data', fontsize=fontsize, va="center", ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0, 0.5),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=0))
plt.text(xmin_ax + 4, ymax_ax - 0.2, '$M=3$', fontsize=fontsize2, ha='left', va='top')
plt.axis('off')

### X(e^{j(w-2pi)/2})
ax = plt.subplot2grid((5, 3), (2, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax - dymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Xsub2 * d_heigh, 'k-', lw=1.5)
# etiquetas de los ejes
for i in np.arange(-nc, nc + 1):
    omega_i = i * omega1 / 2
    plt.plot([omega_i, omega_i], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(omega_i, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(omega_i, xtm, '${}{}\pi$'.format(sg, np.absolute(i)),
                 fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, ytl], [d_heigh, d_heigh], 'k-', lw=1)
plt.text(ytm, ymax_ax - dymax, '$X(e^{j(\omega-2\pi)/3})$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

### X(e^{j(w-2pi)/2})
ax = plt.subplot2grid((5, 3), (3, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax - dymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Xsub3 * d_heigh, 'k-', lw=1.5)
# etiquetas de los ejes
for i in np.arange(-nc, nc + 1):
    omega_i = i * omega1 / 2
    plt.plot([omega_i, omega_i], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(omega_i, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(omega_i, xtm, '${}{}\pi$'.format(sg, np.absolute(i)),
                 fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, ytl], [d_heigh, d_heigh], 'k-', lw=1)
plt.text(ytm, ymax_ax - dymax, '$X(e^{j(\omega-4\pi)/3})$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

d_heigh = d_heigh / 1.5
### Xd(e^{jw})
ax = plt.subplot2grid((5, 3), (4, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax - dymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Xsub * d_heigh, 'k-', lw=1.5)
plt.plot(omega, Xsub1 * d_heigh, color=grey, lw=1, zorder=1)
plt.plot(omega, Xsub2 * d_heigh, color=grey, lw=1, zorder=1)
plt.plot(omega, Xsub3 * d_heigh, color=grey, lw=1, zorder=1)
# etiquetas de los ejes
for i in np.arange(-nc, nc + 1):
    omega_i = i * omega1 / 2
    plt.plot([omega_i, omega_i], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(omega_i, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(omega_i, xtm, '${}{}\pi$'.format(sg, np.absolute(i)),
                 fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh + 0.15, '$\dfrac{1}{MT}=\dfrac{1}{3T}$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, ytl], [d_heigh, d_heigh], 'k-', lw=1)
plt.text(ytm, ymax_ax - dymax, '$X_d(e^{j\omega})$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

plt.savefig('sampling_subsampling_freq_domain_aliasing.pdf', bbox_inches='tight')

plt.show()