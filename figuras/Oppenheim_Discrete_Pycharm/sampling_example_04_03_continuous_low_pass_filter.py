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


def create_impulse_train_spectrum(Nomega, X, NomegaN):
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
T2 = 0.08
# Correpondientes frecuencias de muestreo
omega1 = 2 * math.pi / T1
omega2 = 2 * math. pi / T2
# rango de frecuencias del eje horizontal
fmin = -1.5 * omega1
fmax = 1.5 * omega1
# frecuencia máxima del espectro de la señal continua
omegaN = 0.4 * omega1
# cantidad de frecuencias en el eje horizontal
Nfrec = 2001
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
Xs1 = create_impulse_train_spectrum(Nomega1, X, NomegaN)
Xs2 = create_impulse_train_spectrum(Nomega2, X, NomegaN)

# Filtro pasabajos discreto con T1
# frecuencia de corte discreta. equivale a pi/2
omegac1 = omega1 / 4
# se redefine que sea
# un número entero de muestras de frecuencia
Nomegac1 = round(omegac1 / Domega)
omegac1 = Nomegac1 * Domega
H1 = np.ones((2 * Nomegac1 + 1,))
H1 = create_impulse_train_spectrum(Nomega1, H1, Nomegac1)
Yr1 = np.zeros(omega.shape)
Yr1[if0 - Nomegac1: if0 + Nomegac1 + 1] = X[NomegaN - Nomegac1: NomegaN + Nomegac1 + 1]

# Filtro pasabajos discreto con T2
omegac2 = omegac1 * T1 / T2
# se redefine que sea
# un número entero de muestras de frecuencia
Nomegac2 = round(omegac2 / Domega)
omegac2 = Nomegac2 * Domega
H2 = np.ones((2 * Nomegac2 + 1,))
H2 = create_impulse_train_spectrum(Nomega2, H2, Nomegac2)
# Construcción de Yr
Yr2 = np.zeros(omega.shape)
Yr2[if0 - Nomegac2: if0 + Nomegac2 + 1] = X[NomegaN - Nomegac2: NomegaN + Nomegac2 + 1]

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
ytm = 4
ytm2 = -6
grey = [0.7, 0.7, 0.7]

fig = plt.figure(0, figsize=(9.5, 11), frameon=False)
### Xc(jw)
ax = plt.subplot2grid((6, 3), (0, 0), rowspan=1, colspan=3)
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
plt.text(xmin_ax + 4, ymax_ax - 0.2, '$T=T_1$\n$\omega_c=\dfrac{\pi}{2}$', fontsize=fontsize2, ha='left', va='top')
plt.axis('off')

### X_s(jw) con T1
Nomega = Nomega1
omega_a = omega1
nc = math.floor(((Nfrec - 1) / 2) / Nomega)
d_heigh = T2 / T1 * 0.8

ax = plt.subplot2grid((6, 3), (1, 0), rowspan=1, colspan=3)
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
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(i * omega_a, xtm2, '${}\dfrac{{{}\pi}}{{T}}$'.format(sg, 2 * np.absolute(i)),
                 fontsize=fontsize, ha='center', va='baseline')
        omega_aa = (i / np.absolute(i)) * (np.absolute(i) - 0.5) * omega_a
        plt.plot([omega_aa, omega_aa], [0, xtl], 'k-', lw=1)
        if np.absolute(i) == 1:
            plt.text(omega_aa, xtm2, '${}\dfrac{{\pi}}{{T}}$'.format(sg),
                     fontsize=fontsize, ha='center', va='baseline')
        else:
            plt.text(omega_aa, xtm2, '${}\dfrac{{{}\pi}}{{T}}$'.format(sg, 2 * np.absolute(i)),
                     fontsize=fontsize, ha='center', va='baseline')

# etiquetas
plt.text(omegaN, xtm, '$\Omega_N$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_s(j\Omega)=X(e^{j\Omega T})$', fontsize=fontsize2, ha='left', va='center')
plt.annotate('$\dfrac{2\pi}{T}-\Omega_N$', xytext=(omega_a - omegaN - 3, 0.5), xycoords='data',
             xy=(omega_a - omegaN, 0), textcoords='data', fontsize=fontsize, va="center", ha="right",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.9, 0.1),
                             patchA=None, patchB=None, shrinkA=1, shrinkB=0))
plt.axis('off')

### X(e^{jw}) con T1
ax = plt.subplot2grid((6, 3), (2, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Xs1 * d_heigh, 'k-', lw=1.5)
plt.plot(omega, H1, 'k-', lw=1)
for i in np.arange(-nc, nc + 1):
    plt.plot([i * omega_a, i * omega_a], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(i * omega_a, xtm, '${}{}\pi$'.format(sg, 2 * np.absolute(i)),
                 fontsize=fontsize, ha='center', va='baseline')
        omega_aa = (i / np.absolute(i)) * (np.absolute(i) - 0.5) * omega_a
        plt.plot([omega_aa, omega_aa], [0, xtl], 'k-', lw=1)
        if np.absolute(i) == 1:
            plt.text(omega_aa, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
        else:
            plt.text(omega_aa, xtm, '${}{}\pi$'.format(sg, 2 * np.absolute(i)),
                     fontsize=fontsize, ha='center', va='baseline')

# etiquetas
plt.text(omegaN, xtm, '$\Omega_N T$', fontsize=fontsize, ha='center', va='baseline')
plt.text(omegac1, xtm, '$\omega_c$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omegac1, xtm, '$-\omega_c$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh + 0.1, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X(e^{j\omega})$', fontsize=fontsize2, ha='left', va='center')
plt.annotate('$2\pi-\Omega_N T$', xytext=(omega_a - omegaN - 3, 0.5), xycoords='data',
             xy=(omega_a - omegaN, 0), textcoords='data', fontsize=fontsize, va="center", ha="right",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.9, 0.1),
                             patchA=None, patchB=None, shrinkA=1, shrinkB=0))
plt.annotate('$-\Omega_N T$', xytext=(-omegaN - 3, 0.5), xycoords='data',
             xy=(-omegaN, 0), textcoords='data', fontsize=fontsize, va="center", ha="right",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.9, 0.1),
                             patchA=None, patchB=None, shrinkA=1, shrinkB=0))
plt.annotate('$H(e^{j\omega})$', xytext=(omegac1 + 8, 1.3), xycoords='data', xy=(omegac1, 1),
             textcoords='data', fontsize=fontsize, va="center", ha="left",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0, 0.1),
                             patchA=None, patchB=None, shrinkA=1, shrinkB=0))
plt.axis('off')

### Y(e^{jw}) con T1
ax = plt.subplot2grid((6, 3), (3, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, H1 * Xs1 * d_heigh, 'k-', lw=1.5)
for i in np.arange(-nc, nc + 1):
    plt.plot([i * omega_a, i * omega_a], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(i * omega_a, xtm, '${}{}\pi$'.format(sg, 2 * np.absolute(i)),
                 fontsize=fontsize, ha='center', va='baseline')
        omega_aa = (i / np.absolute(i)) * (np.absolute(i) - 0.5) * omega_a
        plt.plot([omega_aa, omega_aa], [0, xtl], 'k-', lw=1)
        if np.absolute(i) == 1:
            plt.text(omega_aa, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
        else:
            plt.text(omega_aa, xtm, '${}{}\pi$'.format(sg, 2 * np.absolute(i)),
                     fontsize=fontsize, ha='center', va='baseline')

# etiquetas
plt.text(omegac1, xtm, '$\omega_c$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omegac1, xtm, '$-\omega_c$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$Y(e^{j\omega})$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

### Ys(jw) con T1
ax = plt.subplot2grid((6, 3), (4, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, H1 * Xs1 * d_heigh, 'k-', lw=1.5)
for i in np.arange(-nc, nc + 1):
    plt.plot([i * omega_a, i * omega_a], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(i * omega_a, xtm2, '${}\dfrac{{{}\pi}}{{T}}$'.format(sg, 2 * np.absolute(i)),
                 fontsize=fontsize, ha='center', va='baseline')
        omega_aa = (i / np.absolute(i)) * (np.absolute(i) - 0.5) * omega_a
        plt.plot([omega_aa, omega_aa], [0, xtl], 'k-', lw=1)
        if np.absolute(i) == 1:
            plt.text(omega_aa, xtm2, '${}\dfrac{{\pi}}{{T}}$'.format(sg),
                     fontsize=fontsize, ha='center', va='baseline')
        else:
            plt.text(omega_aa, xtm2, '${}\dfrac{{{}\pi}}{{T}}$'.format(sg, 2 * np.absolute(i)),
                     fontsize=fontsize, ha='center', va='baseline')

# pasabajos ideal
plt.plot([-omega_a / 2, -omega_a / 2], [0, 1 / d_heigh], 'k--', lw=1)
plt.plot([omega_a / 2, omega_a / 2], [0, 1 / d_heigh], 'k--', lw=1)
plt.plot([-omega_a / 2, omega_a / 2], [1 / d_heigh, 1 / d_heigh], 'k--', lw=1)

# etiquetas
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$Y(e^{j\Omega T})$', fontsize=fontsize2, ha='left', va='center')
plt.annotate('$H_r(j\Omega)$', xytext=(omega_a / 2 + 8, 1), xycoords='data', xy=(omega_a / 2, 1 / d_heigh),
             textcoords='data', fontsize=fontsize, va="center", ha="left",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0, 0.1),
                             patchA=None, patchB=None, shrinkA=1, shrinkB=0))
plt.text(omegac1, xtm2, '$\dfrac{\omega_c}{T}$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omegac1, xtm2, '$-\dfrac{\omega_c}{T}$', fontsize=fontsize, ha='center', va='baseline')
plt.annotate('$T$', xytext=(18, 1.2), xycoords='data', xy=(0, 1 / d_heigh),
             textcoords='data', fontsize=fontsize, va="center", ha="left",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0, 0.1),
                             patchA=None, patchB=None, shrinkA=1, shrinkB=0))
plt.axis('off')

### Ys(jw) con T1
ax = plt.subplot2grid((6, 3), (5, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
dymax = 0.2
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax - dymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Yr1, 'k-', lw=1.5)

# etiquetas
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax - dymax, '$Y_r(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.text(omegac1, xtm2, '$\dfrac{\omega_c}{T}$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omegac1, xtm2, '$-\dfrac{\omega_c}{T}$', fontsize=fontsize, ha='center', va='baseline')
plt.axis('off')

plt.savefig('sampling_example_04_03_continuous_low_pass_filter.pdf', bbox_inches='tight')

#############################################
#############################################
# Con aliasing
#############################################
#############################################

fig = plt.figure(1, figsize=(9.5, 11 * 5 / 6), frameon=False)
ax = plt.subplot2grid((5, 3), (0, 0), rowspan=1, colspan=3)
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
plt.text(xmin_ax + 4, ymax_ax - 0.2, '$T>T_1$\n$\omega_c=\dfrac{\pi}{2}$', fontsize=fontsize2, ha='left', va='top')
plt.axis('off')

### X_s(jw) con T2
Nomega = Nomega2
omega_a = omega2
nc = math.floor(((Nfrec - 1) / 2) / Nomega)
d_heigh = 0.8

ax = plt.subplot2grid((5, 3), (1, 0), rowspan=1, colspan=3)
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
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(i * omega_a, xtm2, '${}\dfrac{{{}\pi}}{{T}}$'.format(sg, 2 * np.absolute(i)),
                 fontsize=fontsize, ha='center', va='baseline')
        omega_aa = (i / np.absolute(i)) * (np.absolute(i) - 0.5) * omega_a
        plt.plot([omega_aa, omega_aa], [0, xtl], 'k-', lw=1)
        if np.absolute(i) == 1:
            plt.text(omega_aa, xtm2, '${}\dfrac{{\pi}}{{T}}$'.format(sg),
                     fontsize=fontsize, ha='center', va='baseline')
        else:
            plt.text(omega_aa, xtm2, '${}\dfrac{{{}\pi}}{{T}}$'.format(sg, 2 * np.absolute(i) - 1),
                     fontsize=fontsize, ha='center', va='baseline')
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
plt.text(ytm, ymax_ax, '$X_s(j\Omega)=X(e^{j\Omega T})$', fontsize=fontsize2, ha='left', va='center')
plt.annotate('$\dfrac{2\pi}{T}-\Omega_N$', xytext=(omega_a - omegaN, 0.8), xycoords='data',
             xy=(omega_a - omegaN, 0), textcoords='data', fontsize=fontsize, va="center", ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.5, 0.1),
                             patchA=None, patchB=None, shrinkA=1, shrinkB=0))
plt.annotate('$\Omega_N$', xytext=(omegaN + 10, -0.4), xycoords='data',
             xy=(omegaN, 0), textcoords='data', fontsize=fontsize, va="center", ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.1, 0.9),
                             patchA=None, patchB=None, shrinkA=3, shrinkB=0))

plt.axis('off')

### X(e^{jw}) con T2
ax = plt.subplot2grid((5, 3), (2, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Xs2 * d_heigh, 'k-', lw=1.5)
plt.plot(omega, H2, 'k-', lw=1)
for i in np.arange(-nc, nc + 1):
    plt.plot([i * omega_a, i * omega_a], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(i * omega_a, xtm, '${}{}\pi$'.format(sg, 2 * np.absolute(i)),
                 fontsize=fontsize, ha='center', va='baseline')
        omega_aa = (i / np.absolute(i)) * (np.absolute(i) - 0.5) * omega_a
        plt.plot([omega_aa, omega_aa], [0, xtl], 'k-', lw=1)
        if np.absolute(i) == 1:
            plt.text(omega_aa, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
        else:
            plt.text(omega_aa, xtm, '${}{}\pi$'.format(sg, 2 * np.absolute(i) - 1),
                     fontsize=fontsize, ha='center', va='baseline')
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
doc = 4
plt.text(omegac2, xtm, '$\omega_c$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omegac2, xtm, '$-\omega_c$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(-ytm, 1, '$1$', fontsize=fontsize, ha='right', va='bottom')
plt.text(ytm, ymax_ax, '$X(e^{j\omega})$', fontsize=fontsize2, ha='left', va='center')
plt.annotate('$H(e^{j\omega})$', xytext=(omegac2 + 8, 1.3), xycoords='data', xy=(omegac2, 1),
             textcoords='data', fontsize=fontsize, va="center", ha="left",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0, 0.1),
                             patchA=None, patchB=None, shrinkA=1, shrinkB=0))
plt.annotate('$\Omega_N T$', xytext=(omegaN + 10, -0.4), xycoords='data',
             xy=(omegaN, 0), textcoords='data', fontsize=fontsize, va="center", ha="center",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.1, 0.9),
                             patchA=None, patchB=None, shrinkA=4, shrinkB=0))
plt.annotate('$2\pi-\Omega_N T$', xytext=(omega_a - omegaN + 4, -0.7), xycoords='data',
             xy=(omega_a - omegaN, 0), textcoords='data', fontsize=fontsize, va="center", ha="left",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.1, 0.9),
                             patchA=None, patchB=None, shrinkA=4, shrinkB=0))
plt.axis('off')

### Y(e^{jw}) con T2
ax = plt.subplot2grid((5, 3), (3, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, H2 * Xs2 * d_heigh, 'k-', lw=1.5)
for i in np.arange(-nc, nc + 1):
    plt.plot([i * omega_a, i * omega_a], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(i * omega_a, xtm, '${}{}\pi$'.format(sg, 2 * np.absolute(i)),
                 fontsize=fontsize, ha='center', va='baseline')
        omega_aa = (i / np.absolute(i)) * (np.absolute(i) - 0.5) * omega_a
        plt.plot([omega_aa, omega_aa], [0, xtl], 'k-', lw=1)
        if np.absolute(i) == 1:
            plt.text(omega_aa, xtm, '${}\pi$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
        else:
            plt.text(omega_aa, xtm, '${}{}\pi$'.format(sg, 2 * np.absolute(i) - 1),
                     fontsize=fontsize, ha='center', va='baseline')

# etiquetas
plt.text(omegac2, xtm, '$\omega_c$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omegac2, xtm, '$-\omega_c$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, d_heigh, '$\dfrac{1}{T}$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$Y(e^{j\omega})$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

### Ys(jw) con T2
ax = plt.subplot2grid((5, 3), (4, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax - dymax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(omega, Yr2, 'k-', lw=1.5)

# etiquetas
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytm2, 1, '$1$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax - dymax, '$Y_r(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.text(omegac2, xtm2, '$\dfrac{\omega_c}{T}$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omegac2, xtm2, '$-\dfrac{\omega_c}{T}$', fontsize=fontsize, ha='center', va='baseline')
plt.axis('off')

plt.savefig('sampling_example_04_03_continuous_low_pass_filter_aliasing.pdf', bbox_inches='tight')

plt.show()
