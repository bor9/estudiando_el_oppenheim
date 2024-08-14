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
fmin = -1.6 * omega1
fmax = 1.6 * omega1
# frecuencia máxima del espectro de la señal continua
omegaN = 0.5 * omega1
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

# Construcción del espectro de la señal sobremuestreada
L = 2
X = create_signal_spectrum(NomegaN // L)
# señal submuestreada
Xe = create_impulse_train_spectrum(NomegaN, X, NomegaN // L)
Xup = create_impulse_train_spectrum(Nomega1, X, NomegaN // L)

# Construcción del espectro del pasabajos
H = np.ones((NomegaN + 1, ))
H = create_impulse_train_spectrum(Nomega1, H, NomegaN // L)

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
### Xc(jw)
ax = plt.subplot2grid((5, 3), (0, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))

plt.plot(omega, Xc, 'k-', lw=1.5)

# etiquetas
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omegaN, xtm, '$-\Omega_N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(omegaN, xtm, '$\Omega_N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.plot([0, ytl], [1, 1], 'k-', lw=1)
plt.text(ytm2, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_c(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

### X_s(jw) con T1
Nomega = Nomega1
omega_a = omega1

nc = math.floor(((Nfrec - 1) / 2) / (Nomega / 2))
d_heigh = 0.8

### X(e^{jw})
ax = plt.subplot2grid((5, 3), (1, 0), rowspan=1, colspan=3)
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
    omega_i = i * omega_a / 2
    plt.plot([omega_i, omega_i], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        if i == -1: # esto es para omitir el valor en pi
            plt.text(omega_i, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(omega_i, xtm, '${}{}\pi$'.format(sg, np.absolute(i)), fontsize=fontsize, ha='center', va='baseline')
plt.text(omegaN, xtm, '$\omega_N=\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh + 0.1, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, ytl], [d_heigh, d_heigh], 'k-', lw=1)
plt.text(ytm, ymax_ax - dymax, '$X(e^{j\omega})$', fontsize=fontsize2, ha='left', va='center')
plt.text(xmin_ax + 4, ymax_ax - 0.2, '$\Omega_s=2\Omega_N$', fontsize=fontsize2, ha='left', va='baseline')
plt.axis('off')

### Xe(e^{jw})
ax = plt.subplot2grid((5, 3), (2, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax - dymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Xe * d_heigh, 'k-', lw=1.5)
# etiquetas de los ejes
for i in np.arange(-nc, nc + 1):
    omega_i = i * omega_a / 2
    plt.plot([omega_i, omega_i], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(omega_i, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(omega_i, xtm, '${}{}\pi$'.format(sg, np.absolute(i)), fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, ytl], [d_heigh, d_heigh], 'k-', lw=1)
plt.text(ytm, ymax_ax - dymax, '$X_e(e^{j\omega})=X(e^{j\omega L})$', fontsize=fontsize2, ha='left', va='center')
plt.annotate('$\dfrac{\omega_N}{L}=\dfrac{\pi}{2}$', xytext=(omegaN / 2, 0.8), xycoords='data',
             xy=(omegaN / 2, 0), textcoords='data', fontsize=fontsize, va="center", ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 0),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=0))
plt.text(xmin_ax + 4, ymax_ax - 0.2, '$L=2$', fontsize=fontsize2, ha='left', va='baseline')
plt.axis('off')


### H(e^{jw})
d_heigh_H = 1.2
ax = plt.subplot2grid((5, 3), (3, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax - dymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, H * d_heigh_H, 'k-', lw=1.5)
# etiquetas de los ejes
for i in np.arange(-nc, nc + 1):
    omega_i = i * omega_a / 2
    plt.plot([omega_i, omega_i], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(omega_i, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(omega_i, xtm, '${}{}\pi$'.format(sg, np.absolute(i)), fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh_H, '$L$', fontsize=fontsize, ha='right', va='bottom')
plt.text(ytm, ymax_ax - dymax, '$H(e^{j\omega})$', fontsize=fontsize2, ha='left', va='center')
plt.annotate('$\omega_c=\dfrac{\pi}{L}$', xytext=(omegaN, 0.8), xycoords='data',
             xy=(omegaN / 2, 0), textcoords='data', fontsize=fontsize, va="center", ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0, 0),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=0))
plt.axis('off')

### Xd(e^{jw})
d_heigh = 1.3

ax = plt.subplot2grid((5, 3), (4, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax - dymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Xup * d_heigh, 'k-', lw=1.5)
# etiquetas de los ejes
for i in np.arange(-nc, nc + 1):
    omega_i = i * omega_a / 2
    plt.plot([omega_i, omega_i], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(omega_i, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(omega_i, xtm, '${}{}\pi$'.format(sg, np.absolute(i)), fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{L}{T}=\dfrac{1}{T_i}$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, ytl], [d_heigh, d_heigh], 'k-', lw=1)
plt.text(ytm, ymax_ax - dymax, '$X_i(e^{j\omega})$', fontsize=fontsize2, ha='left', va='center')

plt.text(omegaN / L, xtm2, '$\omega_N=\dfrac{\pi}{L}$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omegaN / L, xtm2, '$-\dfrac{\pi}{L}$', fontsize=fontsize, ha='center', va='baseline')

plt.axis('off')

plt.savefig('sampling_upsampling_freq_domain.pdf', bbox_inches='tight')

plt.show()