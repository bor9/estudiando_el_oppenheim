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
T1 = 0.05
T2 = 0.083
# Correpondientes frecuencias de muestreo
omega1 = 2 * math.pi / T1
omega2 = 2 *math. pi / T2
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
Nomega2 = round(omega2 / Domega)
omega2 = Nomega2 * Domega
NomegaN = round(omegaN / Domega)
omegaN = NomegaN * Domega

# Construcción del espectro
X = create_signal_spectrum(NomegaN)

# Construcción de Xc
Xc = np.zeros(omega.shape)
Xc[if0 - NomegaN: if0 + NomegaN + 1] = X

# Construcción de Xs1 y Xs2
Xs1 = create_impulse_train_spectrum(Nomega1)
Xs2 = create_impulse_train_spectrum(Nomega2)

# Gráficas

df = 30
xmin_ax = fmin - df
xmax_ax = fmax + df
ymin_ax =-0.1
ymax_ax = 1.8

display_length = 6
fontsize = 10
fontsize2 = 12
# x y ticks labels margin
xtm = -0.2
ytm = 14
ytm2 = -10
grey = [0.5, 0.5, 0.5]

fig = plt.figure(0, figsize=(10, 7), frameon=False)
### Xc(jw)
ax = plt.subplot2grid((3, 6), (0, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
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

### S(jw) con T1
ax = plt.subplot2grid((3, 6), (1, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
Nomega = Nomega1
omega_a = omega1
nc = math.floor(((Nfrec - 1) / 2) / Nomega)
d_heigh = T2 / T1 * 0.8
for i in np.arange(-nc, nc + 1):
    plt.annotate(s='', xytext=(i * omega_a, 0), xy=(i * omega_a, d_heigh),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', facecolor='black',
                                 shrinkA=0, shrinkB=0))
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(i * omega_a, xtm, '$-\Omega_s$' if i < 0 else '$\Omega_s$', fontsize=fontsize, ha='center',
                 va='baseline')
    else:
        plt.text(i * omega_a, xtm, '${}\Omega_s$'.format(i), fontsize=fontsize, ha='center', va='baseline')

# etiquetas
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{2\pi}{T}$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$S(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.text(xmin_ax, ymax_ax, '$T=T_1$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

ax = plt.subplot2grid((3, 6), (2, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
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
plt.annotate('$\Omega_s-\Omega_N$', xytext=(omega_a - omegaN + 6, -0.35), xycoords='data', xy=(omega_a - omegaN, 0),
             textcoords='data', fontsize=fontsize, va="center", ha="left",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.08, 0.9),
                             patchA=None, patchB=None, shrinkA=4, shrinkB=0))


plt.axis('off')

##########
# T>T1
##########

### Xc(jw)
ax = plt.subplot2grid((3, 6), (0, 3), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
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

### S(jw) con T>T1
ax = plt.subplot2grid((3, 6), (1, 3), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
Nomega = Nomega2
omega_a = omega2
nc = math.floor(((Nfrec - 1) / 2) / Nomega)
d_heigh = 0.8
for i in np.arange(-nc, nc + 1):
    plt.annotate(s='', xytext=(i * omega_a, 0), xy=(i * omega_a, d_heigh),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', facecolor='black',
                                 shrinkA=0, shrinkB=0))
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(i * omega_a, xtm, '$-\Omega_s$' if i < 0 else '$\Omega_s$', fontsize=fontsize, ha='center',
                 va='baseline')
    else:
        plt.text(i * omega_a, xtm, '${}\Omega_s$'.format(i), fontsize=fontsize, ha='center', va='baseline')

# etiquetas
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{2\pi}{T}$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$S(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.text(xmin_ax, ymax_ax, '$T>T_1$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

ax = plt.subplot2grid((3, 6), (2, 3), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Xs2 * d_heigh, 'k-', lw=1.5)
for i in np.arange(-nc, nc + 1):
    plt.plot([i * omega_a, i * omega_a], [0, xtl], 'k-', lw=1)
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(i * omega_a, xtm, '$-\Omega_s$' if i < 0 else '$\Omega_s$', fontsize=fontsize, ha='center',
                 va='baseline')
    else:
        plt.text(i * omega_a, xtm, '${}\Omega_s$'.format(i), fontsize=fontsize, ha='center', va='baseline')
# graficas de los  espectros sin solaparse
for i in np.arange(-nc, nc + 1):
    # indice de la muestra central del espectro
    ic = if0 + i * Nomega
    # indices de la muestra inicial y final del espectro
    ii = ic - NomegaN
    ie = ic + NomegaN + 1
    if ii < 0:
        plt.plot(omega[0: ie], d_heigh * X[-ii:], color=grey, lw=1, zorder=1)
    elif ie > Nfrec + 1:
        plt.plot(omega[ii: Nfrec + 1], d_heigh * X[: -(ie - Nfrec)], color=grey, lw=1, zorder=1)
    else:
        plt.plot(omega[ii: ie], d_heigh * X, color=grey, lw=1, zorder=1)

# etiquetas
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_s(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.annotate('$\Omega_N$', xytext=(omegaN + 4, -0.35), xycoords='data', xy=(omegaN, 0),
             textcoords='data', fontsize=fontsize, va="center", ha="left",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.08, 0.9),
                             patchA=None, patchB=None, shrinkA=4, shrinkB=0))
plt.annotate('$\Omega_s-\Omega_N$', xytext=(omega_a - omegaN + 9, -0.35), xycoords='data', xy=(omega_a - omegaN, 0),
             textcoords='data', fontsize=fontsize, va="center", ha="right",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.8, 0.9),
                             patchA=None, patchB=None, shrinkA=4, shrinkB=0))
plt.axis('off')

plt.savefig('sampling_impulse_modulator_frequency_representation.pdf', bbox_inches='tight')

plt.show()
