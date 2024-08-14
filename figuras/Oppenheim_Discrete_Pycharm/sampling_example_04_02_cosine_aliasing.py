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
# frecuencia de muestreo y frecuencia de corte del pasabajos
fs = 6000
Ts = 1 / fs
omega_s = 2 * np.pi / Ts
omega_c = omega_s / 2

f0 = 8000
omega_0 = 2 * f0 * np.pi

fmin = -1.5 * omega_s
fmax = 1.5 * omega_s

# Gráficas

# Altura de las flechas. Esto es pi.
a_height = 1
T = 1.3

df = 30
xmin_ax = fmin - df
xmax_ax = fmax + df
ymin_ax =-0.1
ymax_ax = 1.8

display_length = 5
fontsize = 10
fontsize2 = 12
# x y ticks labels margin
xtm = -0.23
ytm = 1300
ytm2 = -900
grey = [0.5, 0.5, 0.5]

d_ymax = 0.15

fig = plt.figure(0, figsize=(8, 6), frameon=False)
### Xc(jw)
ax = plt.subplot2grid((3, 6), (0, 0), rowspan=1, colspan=6)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + d_ymax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))

# deltas
plt.annotate(s='', xytext=(omega_0, 0), xy=(omega_0, a_height),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='black',
                                 shrinkA=0, shrinkB=0))
plt.annotate(s='', xytext=(-omega_0, 0), xy=(-omega_0, a_height),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='r',
                                 shrinkA=0, shrinkB=0))

# etiquetas
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omega_0, xtm, '$-{:d}\pi$'.format(2 * f0), fontsize=fontsize, ha='center', va='baseline')
plt.text(omega_0, xtm, '${:d}\pi$'.format(2 * f0), fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.plot([0, ytl], [a_height, a_height], 'k-', lw=1)
plt.text(ytm2, a_height, '$\pi$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_c(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

### X(e^{jw})
ax = plt.subplot2grid((3, 6), (1, 0), rowspan=1, colspan=6)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + d_ymax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
# cantidad de periodos en el eje
nper = int(np.ceil(fmax / omega_s))
for i in np.arange(-nper, nper + 1):
    omega = i * omega_s + omega_0
    if (omega > fmin) & (omega < fmax):
        plt.annotate(s='', xytext=(omega, 0), xy=(omega, a_height / T),
                     arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='black',
                                     shrinkA=0, shrinkB=0))
        plt.text(i * omega_s + omega_0, xtm, '${:d}\pi$'.format(2 * (i * fs + f0)), fontsize=fontsize,
                 ha='center', va='baseline')
    omega = i * omega_s - omega_0
    if (omega > fmin) & (omega < fmax):
        plt.annotate(s='', xytext=(i * omega_s - omega_0, 0), xy=(i * omega_s - omega_0, a_height / T),
                    arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='r',
                                    shrinkA=0, shrinkB=0))
        plt.text(i * omega_s - omega_0, xtm, '${:d}\pi$'.format(2 * (i * fs - f0)), fontsize=fontsize,
                 ha='center', va='baseline')
    omega = i * omega_s
    if (omega > fmin) & (omega < fmax):
        plt.plot([i * omega_s, i * omega_s], [0, xtl], 'k-', lw=1)
        plt.text(i * omega_s, xtm, '${:d}\pi$'.format(i * 2 * fs) if i != 0 else '$0$', fontsize=fontsize,
                ha='center', va='baseline')

# pasabajos ideal
plt.plot([omega_c, omega_c], [0, T], 'k--', lw=1)
plt.plot([-omega_c, -omega_c], [0, T], 'k--', lw=1)
plt.plot([-omega_c, omega_c], [T, T], 'k--', lw=1)
dx = 1.8
plt.text(omega_c, dx * xtm, '${:d}\pi$'.format(fs), fontsize=fontsize, ha='center', va='baseline')
plt.text(-omega_c, dx * xtm, '${:d}\pi$'.format(-fs), fontsize=fontsize, ha='center', va='baseline')

dy = 0.08
plt.text(ytm2, T + dy, '$T$', fontsize=fontsize, ha='right', va='baseline')
plt.text(omega_c, T + dy, '$H_r(j\Omega)$', fontsize=fontsize, ha='center', va='baseline')

plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.plot([0, ytl], [a_height / T, a_height / T], 'k', lw=1)
plt.text(ytm2, a_height / T, '$\dfrac{\pi}{T}$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_s(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

xtm2 = -0.36
ax = plt.subplot2grid((3, 6), (2, 0), rowspan=1, colspan=6)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + d_ymax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
# cantidad de periodos en el eje
nper = int(np.ceil(fmax / omega_s))
for i in np.arange(-nper, nper + 1):
    omega = i * omega_s + omega_0
    if (omega > fmin) & (omega < fmax):
        plt.annotate(s='', xytext=(omega, 0), xy=(omega, a_height),
                    arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='black',
                                    shrinkA=0, shrinkB=0))
        k = int(omega / fs)
        plt.text(omega, xtm2,
                 '$\dfrac{{{:d}\pi}}{{3}}$'.format(k) if k > 0 else '$-\dfrac{{{:d}\pi}}{{3}}$'.format(-k),
                 fontsize=fontsize, ha='center', va='baseline')
    omega = i * omega_s - omega_0
    if (omega > fmin) & (omega < fmax):
        plt.annotate(s='', xytext=(omega, 0), xy=(omega, a_height),
                    arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='r',
                                    shrinkA=0, shrinkB=0))
        k = int(omega / fs)
        plt.text(i * omega_s - omega_0, xtm2,
                 '$\dfrac{{{:d}\pi}}{{3}}$'.format(k) if k > 0 else '$-\dfrac{{{:d}\pi}}{{3}}$'.format(-k),
                 fontsize=fontsize, ha='center', va='baseline')
    omega = i * omega_s
    if (omega > fmin) & (omega < fmax):
        plt.plot([i * omega_s, i * omega_s], [0, xtl], 'k-', lw=1)
        plt.text(i * omega_s, xtm, '${:d}\pi$'.format(i * 2) if i != 0 else '$0$', fontsize=fontsize,
                 ha='center', va='baseline')


plt.text(omega_c, xtm, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omega_c, xtm, '$-\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.plot([omega_s / 2, omega_s / 2], [0, xtl], 'k-', lw=1)
plt.plot([-omega_s / 2, -omega_s / 2], [0, xtl], 'k-', lw=1)

plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='right', va='baseline')
plt.plot([0, ytl], [a_height, a_height], 'k', lw=1)
plt.text(ytm2, a_height, '$\pi$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X(e^{j\omega})=X_s(j\omega/T)$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

plt.savefig('sampling_example_04_02_cosine_aliasing.pdf', bbox_inches='tight')

plt.show()
