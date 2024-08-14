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


def create_impulse_train_spectrum(Nomega):
    # relleno de ceros para simplificar la construcción
    zpad = Nomega + 2 * NomegaN + 1
    Xs = np.zeros((Nfrec + 2 * zpad,))
    if0_aux = if0 + zpad
    # número de copias del espectro
    nc = math.ceil(((Nfrec - 1) / 2) / Nomega)
    for i in np.arange(-nc, nc + 1):
        # indice de la frecuencia central
        ic = if0_aux + i * Nomega
        Xs[ic - NomegaN: ic + NomegaN + 1] += X
    # se elimina el relleno de ceros
    return Xs[zpad: Nfrec + zpad]

#
# Parámetros
#
# Períodos de muestreo
T1 = 0.75
# Correspondientes frecuencias de muestreo
omega1 = 2 * math.pi / T1
# rango de frecuencias del eje horizontal
fmin = -2.5 * omega1
fmax = 2.5 * omega1
# frecuencia máxima del espectro de la señal continua
omegaN = 0.8 * omega1 / 2
# cantidad de frecuencias en el eje horizontal
Nfrec = 1001
omega = np.linspace(fmin, fmax, Nfrec, endpoint=True)
# índice del elemento de frecuencia 0
if0 = int((Nfrec - 1) / 2)
# intervalo de frecuencia entre dos muestras en el eje horizontal
Domega = (fmax - fmin) / (Nfrec - 1)
# se redefine omega1, omega2 y omegaN para que sean
# un número entero de muestras de frecuencia
Nomega1 = round(omega1 / Domega)
omega1 = Nomega1 * Domega
NomegaN = round(omegaN / Domega)
omegaN = NomegaN * Domega

# Construcción del espectro
X = create_signal_spectrum(NomegaN)

# Construcción de Xc
Xc = np.zeros(omega.shape)
Xc[if0 - NomegaN: if0 + NomegaN + 1] = X

# Construcción de Xs1 y Xs2
Xs1 = create_impulse_train_spectrum(Nomega1)

# Gráficas

df = 1
xmin_ax = fmin - df
xmax_ax = fmax + df
ymin_ax = -0.1
ymax_ax = 1.8

display_length = 6
fontsize = 10
fontsize2 = 12
# x y ticks labels margin
xtm = -0.28
ytm = 0.5
ytm2 = -0.5
grey = [0.5, 0.5, 0.5]

dy = 0.2

fig = plt.figure(0, figsize=(8, 7), frameon=False)
### Xc(jw)
ax = plt.subplot2grid((4, 6), (0, 0), rowspan=1, colspan=6)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + dy)
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
plt.text(ytm2, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_c(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')


### X_s(jw)
Nomega = Nomega1
omega_a = omega1
nc = math.floor(((Nfrec - 1) / 2) / Nomega)
d_heigh = 1 / T1

ax = plt.subplot2grid((4, 6), (1, 0), rowspan=1, colspan=6)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + dy)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Xs1 * d_heigh, 'k-', lw=1.5)
for i in np.arange(-nc, nc + 1):
    plt.plot([i * omega_a, i * omega_a], [0, xtl], 'k-', lw=1)
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(i * omega_a, xtm, '$-\Omega_s$' if i < 0 else '$\Omega_s$', fontsize=fontsize, ha='center',
                 va='baseline')
    else:
        plt.text(i * omega_a, xtm, '${}\Omega_s$'.format(i), fontsize=fontsize, ha='center', va='baseline')
# etiquetas
plt.text(omegaN, xtm, '$\Omega_N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omegaN, xtm, '$-\Omega_N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_s(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.annotate('$\Omega_s-\Omega_N$', xytext=(omega_a - omegaN + 0.5, -0.5), xycoords='data', xy=(omega_a - omegaN, 0),
             textcoords='data', fontsize=fontsize, va="center", ha="left",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.08, 0.9),
                             patchA=None, patchB=None, shrinkA=4, shrinkB=0))

plt.axis('off')

### Hr(jw)
omega_c = omega1 / 2
ax = plt.subplot2grid((4, 6), (2, 0), rowspan=1, colspan=6)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + dy)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))

plt.plot([fmin, -omega_c], [0, 0], 'k-', lw=1.5)
plt.plot([omega_c, fmax], [0, 0], 'k-', lw=1.5)
plt.plot([-omega_c, omega_c], [T1, T1], 'k-', lw=1.5)
plt.plot([-omega_c, -omega_c], [0, T1], 'k-', lw=1.5)
plt.plot([omega_c, omega_c], [0, T1], 'k-', lw=1.5)

# etiquetas
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omega_c, xtm, '$-\Omega_c$', fontsize=fontsize, ha='center', va='baseline')
plt.text(omega_c, xtm, '$\Omega_c$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, 1 / d_heigh + 0.05, '$T$', fontsize=fontsize, ha='right', va='bottom')
plt.text(ytm, ymax_ax, '$H_r(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.text(omega_c, 1 / d_heigh + 0.05, '$\Omega_N<\Omega_c<\Omega_s-\Omega_N$',
         fontsize=fontsize, ha='left', va='bottom')
plt.axis('off')

### Xr(jw)
ax = plt.subplot2grid((4, 6), (3, 0), rowspan=1, colspan=6)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + dy)
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
plt.text(ytm2, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_r(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

plt.savefig('sampling_impulse_modulator_frequency_representation_reconstruction.pdf', bbox_inches='tight')

plt.show()
