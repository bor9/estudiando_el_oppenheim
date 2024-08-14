import matplotlib.pyplot as plt
import numpy as np

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


#
# Parámetros
#
# Esto es mejor dejarlo fijo
# frecuencia y tiempo máximo de los ejes
N = 4
wmax = N * np.pi
tmax = N
# rango de frecuencias
fmin = -wmax
fmax = wmax
# rango de tiempos
dt = 0.5
tmin = -tmax - dt
tmax = tmax + dt
#
# Fin de parámentros
#

# parámetros del pasabajos
# frecuencia de corte
W = np.pi
# ganancia
G = 1
# construcción de x(t)
Nt = 501
t = np.linspace(tmin, tmax, Nt)
xt = G * np.sin(W * t) / (np.pi * t)
xt[(Nt - 1) // 2] = W / np.pi

dax = 0.4
df = dax * np.pi
fmin_ax = fmin - df
fmax_ax = fmax + df

dt = dax
tmin_ax = tmin - dt
tmax_ax = tmax + dt
ymin_ax = -0.3
ymax_ax = 1.6

#
# Gráficas
#

display_length = 5
fontsize = 11
lw = 2
# x y ticks labels margin
xtm = -0.25
ytmf = 0.5
ytmt = ytmf / np.pi
ytmf2 = 0.3
ytmt2 = ytmf2 / np.pi

d_ymax = 0

fig = plt.figure(0, figsize=(7, 4), frameon=False)
### X(jw)
ax = plt.subplot2grid((2, 6), (0, 0), rowspan=1, colspan=6)
plt.xlim(fmin_ax, fmax_ax)
plt.ylim(ymin_ax, ymax_ax + d_ymax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(fmin_ax, 0), xycoords='data', xy=(fmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot([W, W], [0, G], 'k-', lw=lw)
plt.plot([-W, -W], [0, G], 'k-', lw=lw)
plt.plot([-W, W], [G, G], 'k-', lw=lw)
plt.plot([fmin, -W], [0, 0], 'k-', lw=lw)
plt.plot([W, fmax], [0, 0], 'k-', lw=lw)

# etiquetas
plt.text(-ytmf2, G, '$1$', fontsize=fontsize, ha='right', va='bottom')
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(fmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(-W, xtm, '$-W$', fontsize=fontsize, ha='center', va='baseline')
plt.text(W, xtm, '$W$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytmf, ymax_ax, '$X(j\omega)$', fontsize=fontsize, ha='left', va='center')
plt.axis('off')

### x(t)
ax = plt.subplot2grid((2, 6), (1, 0), rowspan=1, colspan=6)
plt.xlim(tmin_ax, tmax_ax)
plt.ylim(ymin_ax, ymax_ax + d_ymax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(tmin_ax, 0), xycoords='data', xy=(tmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(t, xt, 'k-', lw=lw)

xtm2 = 1.6 * xtm
for i in np.arange(-N, N + 1):
    plt.plot([i, i], [0, xtl], 'k-', lw=1)
    sg = '-' if i < 0 else ''
    if i == 0:
        plt.text(i, xtm, '$0$'.format(-i), fontsize=fontsize, ha='center', va='baseline')
    elif np.abs(i) == 1:
        plt.text(i, xtm2, '${:s}\dfrac{{\pi}}{{W}}$'.format(sg), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(i, xtm2, '${:s}\dfrac{{{:d}\pi}}{{W}}$'.format(sg, np.abs(i)), fontsize=fontsize, ha='center',
                 va='baseline')

plt.plot([0, ytl], [G * W / np.pi, G * W / np.pi], 'k-', lw=1)
plt.text(-0.2, G * W / np.pi + 0.1, '$\dfrac{W}{\pi}$', fontsize=fontsize, ha='right', va='center')
plt.text(tmax_ax, xtm, '$t$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytmt, ymax_ax, '$x(t)$', fontsize=fontsize, ha='left', va='center')
plt.axis('off')

plt.savefig('sampling_lowpass_continuous_inverse_fourier_transform_1.pdf', bbox_inches='tight')


fig = plt.figure(1, figsize=(7, 4), frameon=False)
### X(jw)
ax = plt.subplot2grid((2, 6), (0, 0), rowspan=1, colspan=6)
plt.xlim(fmin_ax, fmax_ax)
plt.ylim(ymin_ax, ymax_ax + d_ymax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(fmin_ax, 0), xycoords='data', xy=(fmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot([W, W], [0, G], 'k-', lw=lw)
plt.plot([-W, -W], [0, G], 'k-', lw=lw)
plt.plot([-W, W], [G, G], 'k-', lw=lw)
plt.plot([fmin, -W], [0, 0], 'k-', lw=lw)
plt.plot([W, fmax], [0, 0], 'k-', lw=lw)

# etiquetas
xtm2 = 1.3 * xtm
plt.text(-ytmf2, G, '$T$', fontsize=fontsize, ha='right', va='bottom')
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(fmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(-W, xtm2, '$-\dfrac{\pi}{T}$', fontsize=fontsize, ha='center', va='baseline')
plt.text(W, xtm2, '$\dfrac{\pi}{T}$', fontsize=fontsize, ha='center', va='baseline')
plt.text(ytmf, ymax_ax, '$H_r(j\Omega)$', fontsize=fontsize, ha='left', va='center')
plt.axis('off')

### x(t)
ax = plt.subplot2grid((2, 6), (1, 0), rowspan=1, colspan=6)
plt.xlim(tmin_ax, tmax_ax)
plt.ylim(ymin_ax, ymax_ax + d_ymax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(tmin_ax, 0), xycoords='data', xy=(tmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(t, xt, 'k-', lw=lw)

for i in np.arange(-N, N + 1):
    plt.plot([i, i], [0, xtl], 'k-', lw=1)
    if i == 0:
        plt.text(i, xtm, '$0$'.format(-i), fontsize=fontsize, ha='center', va='baseline')
    elif i == -1:
        plt.text(i + 0.1, xtm, '$-T$', fontsize=fontsize, ha='center', va='baseline')
    elif i == 1:
        plt.text(i - 0.1, xtm, '$T$', fontsize=fontsize, ha='center', va='baseline')
    elif i == -2:
        plt.text(i - 0.2, xtm, '$-{:d}T$'.format(-i), fontsize=fontsize, ha='center', va='baseline')
    elif i == 2:
        plt.text(i + 0.1, xtm, '${:d}T$'.format(i), fontsize=fontsize, ha='center', va='baseline')
    elif i < -1:
        plt.text(i, xtm, '$-{:d}T$'.format(-i), fontsize=fontsize, ha='center', va='baseline')
    else:
        plt.text(i, xtm, '${:d}T$'.format(i), fontsize=fontsize, ha='center', va='baseline')

plt.plot([0, ytl], [G * W / np.pi, G * W / np.pi], 'k-', lw=1)
plt.text(-ytmt2, G * W / np.pi, '$1$', fontsize=fontsize, ha='right', va='bottom')
plt.text(tmax_ax, xtm, '$t$', fontsize=fontsize, ha='right', va='baseline')
plt.text(ytmt, ymax_ax, '$h_r(t)$', fontsize=fontsize, ha='left', va='center')
plt.axis('off')

plt.savefig('sampling_lowpass_continuous_inverse_fourier_transform_2.pdf', bbox_inches='tight')


plt.show()
