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
omega_s = 120
omega_c = omega_s / 2
# frecuencia de las sinusoides
omega1 = 40
omega2 = 80

fmin = -1.8 * omega_s
fmax = 1.8 * omega_s

# Altura de las flechas. Esto es pi.
a_height = 1.25
T = 1.3

# Gráficas

df = 30
xmin_ax = fmin - df
xmax_ax = fmax + df
ymin_ax =-0.1
ymax_ax = 1.8

display_length = 5
fontsize = 10
fontsize2 = 12
# x y ticks labels margin
xtm = -0.2
ytm = 12
ytm2 = -8
grey = [0.5, 0.5, 0.5]

fig = plt.figure(0, figsize=(10, 7), frameon=False)
### Xc(jw) con omega1
ax = plt.subplot2grid((3, 6), (0, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + 0.2)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))

# deltas
plt.annotate(s='', xytext=(omega1, 0), xy=(omega1, a_height),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='black',
                                 shrinkA=0, shrinkB=0))
plt.annotate(s='', xytext=(-omega1, 0), xy=(-omega1, a_height),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='r',
                                 shrinkA=0, shrinkB=0))

# etiquetas
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omega1, xtm, '$-\Omega_0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(omega1, xtm, '$\Omega_0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.plot([0, ytl], [a_height, a_height], 'k-', lw=1)
plt.text(ytm2, a_height, '$\pi$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_c(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

### X_s(jw) con omega1
ax = plt.subplot2grid((3, 6), (1, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + 0.2)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
# cantidad de periodos en el eje
nper = np.floor(fmax / omega_s)
for i in np.arange(-nper, nper + 1):
    plt.annotate(s='', xytext=(i * omega_s + omega1, 0), xy=(i * omega_s + omega1, a_height / T),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='black',
                                 shrinkA=0, shrinkB=0))
    plt.annotate(s='', xytext=(i * omega_s - omega1, 0), xy=(i * omega_s - omega1, a_height / T),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='r',
                                 shrinkA=0, shrinkB=0))
    plt.plot([i * omega_s, i * omega_s], [0, xtl], 'k-', lw=1)
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(i * omega_s, xtm, '$-\Omega_s$' if i < 0 else '$\Omega_s$', fontsize=fontsize, ha='center',
                 va='baseline')
    else:
        plt.text(i * omega_s, xtm, '${}\Omega_s$'.format(i), fontsize=fontsize, ha='center', va='baseline')

# pasabajos ideal
plt.plot([omega_c, omega_c], [0, T], 'k--', lw=1)
plt.plot([-omega_c, -omega_c], [0, T], 'k--', lw=1)
plt.plot([-omega_c, omega_c], [T, T], 'k--', lw=1)
dy = 0.08
plt.text(ytm2, T + dy, '$T$', fontsize=fontsize, ha='right', va='baseline')
plt.text(omega_c, T + dy, '$H_r(j\Omega)$', fontsize=fontsize, ha='center', va='baseline')

# etiquetas
plt.text(-omega1, xtm, '$-\Omega_0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(omega1, xtm, '$\Omega_0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(omega_c * 1.1, -0.3, r'$\dfrac{\Omega_s}{2}$', fontsize=fontsize, ha='center', va='baseline')
plt.annotate('$\Omega_s-\Omega_0$', xytext=(omega_s - omega1 + 12, -0.35), xycoords='data', xy=(omega_s - omega1, 0),
             textcoords='data', fontsize=fontsize, va="center", ha="left",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.08, 0.9),
                             patchA=None, patchB=None, shrinkA=4, shrinkB=0))

plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.plot([0, ytl], [a_height / T, a_height / T], 'k', lw=1)
plt.text(ytm2, a_height / T, '$\dfrac{\pi}{T}$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_s(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.text(xmin_ax, ymax_ax, '$\Omega_0<\dfrac{\Omega_s}{2}$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

ax = plt.subplot2grid((3, 6), (2, 0), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + 0.2)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))

# deltas
plt.annotate(s='', xytext=(omega1, 0), xy=(omega1, a_height),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='black',
                                 shrinkA=0, shrinkB=0))
plt.annotate(s='', xytext=(-omega1, 0), xy=(-omega1, a_height),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='r',
                                 shrinkA=0, shrinkB=0))

# etiquetas
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omega1, xtm, '$-\Omega_0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(omega1, xtm, '$\Omega_0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.plot([0, ytl], [a_height, a_height], 'k-', lw=1)
plt.text(ytm2, a_height, '$\pi$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_r(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

##########
# omega2
##########

ax = plt.subplot2grid((3, 6), (0, 3), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + 0.2)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))

# deltas
plt.annotate(s='', xytext=(omega2, 0), xy=(omega2, a_height),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='black',
                                 shrinkA=0, shrinkB=0))
plt.annotate(s='', xytext=(-omega2, 0), xy=(-omega2, a_height),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='r',
                                 shrinkA=0, shrinkB=0))

# etiquetas
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-omega2, xtm, '$-\Omega_0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(omega2, xtm, '$\Omega_0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.plot([0, ytl], [a_height, a_height], 'k-', lw=1)
plt.text(ytm2, a_height, '$\pi$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_c(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

### X_s(jw) con omega2
ax = plt.subplot2grid((3, 6), (1, 3), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + 0.2)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
# cantidad de periodos en el eje
nper = np.floor(fmax / omega_s)
for i in np.arange(-nper, nper + 1):
    plt.annotate(s='', xytext=(i * omega_s + omega2, 0), xy=(i * omega_s + omega2, a_height / T),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='black',
                                 shrinkA=0, shrinkB=0))
    plt.annotate(s='', xytext=(i * omega_s - omega2, 0), xy=(i * omega_s - omega2, a_height / T),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='r',
                                 shrinkA=0, shrinkB=0))
    plt.plot([i * omega_s, i * omega_s], [0, xtl], 'k-', lw=1)
    if i == 0:
        plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
    elif np.absolute(i) == 1:
        plt.text(i * omega_s, xtm, '$-\Omega_s$' if i < 0 else '$\Omega_s$', fontsize=fontsize, ha='center',
                 va='baseline')
    else:
        plt.text(i * omega_s, xtm, '${}\Omega_s$'.format(i), fontsize=fontsize, ha='center', va='baseline')

# pasabajos ideal
plt.plot([omega_c, omega_c], [0, T], 'k--', lw=1)
plt.plot([-omega_c, -omega_c], [0, T], 'k--', lw=1)
plt.plot([-omega_c, omega_c], [T, T], 'k--', lw=1)
dy = 0.08
plt.text(ytm2, T + dy, '$T$', fontsize=fontsize, ha='right', va='baseline')
plt.text(omega_c, T + dy, '$H_r(j\Omega)$', fontsize=fontsize, ha='center', va='baseline')

# etiquetas
plt.text(-omega2, xtm, '$-\Omega_0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(omega2 * 1.1, xtm, '$\Omega_0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(omega_c, -0.3, r'$\dfrac{\Omega_s}{2}$', fontsize=fontsize, ha='center', va='baseline')
plt.annotate('$\Omega_s-\Omega_0$', xytext=(omega_s - omega2 - 8, -0.35), xycoords='data', xy=(omega_s - omega2, 0),
             textcoords='data', fontsize=fontsize, va="center", ha="right",
             arrowprops=dict(arrowstyle="-|>, head_width=0.1, head_length=0.4", facecolor='black', relpos=(0.85, 0.9),
                             patchA=None, patchB=None, shrinkA=4, shrinkB=0))

plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.plot([0, ytl], [a_height / T, a_height / T], 'k', lw=1)
plt.text(ytm2, a_height / T, '$\dfrac{\pi}{T}$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_s(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.text(xmin_ax, ymax_ax, '$\Omega_0>\dfrac{\Omega_s}{2}$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

ax = plt.subplot2grid((3, 6), (2, 3), rowspan=1, colspan=3)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax + 0.2)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))

# deltas
plt.annotate(s='', xytext=(omega_s - omega2, 0), xy=(omega_s - omega2, a_height),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='r',
                                 shrinkA=0, shrinkB=0))
plt.annotate(s='', xytext=(-(omega_s - omega2), 0), xy=(-(omega_s - omega2), a_height),
                 arrowprops=dict(arrowstyle='->, head_width=0.2, head_length=0.6', color='k',
                                 shrinkA=0, shrinkB=0))

# etiquetas
plt.text(0, xtm, '$0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(omega_s - omega2 + 10, xtm, '$\Omega_s-\Omega_0$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-(omega_s - omega2) - 15, xtm, '$-(\Omega_s-\Omega_0)$', fontsize=fontsize, ha='center', va='baseline')
plt.text(xmax_ax, xtm, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.plot([0, ytl], [a_height, a_height], 'k-', lw=1)
plt.text(ytm2, a_height, '$\pi$', fontsize=fontsize, ha='right', va='center')
plt.text(ytm, ymax_ax, '$X_r(j\Omega)$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')


plt.savefig('sampling_aliasing_sinusoid.pdf', bbox_inches='tight')

plt.show()
