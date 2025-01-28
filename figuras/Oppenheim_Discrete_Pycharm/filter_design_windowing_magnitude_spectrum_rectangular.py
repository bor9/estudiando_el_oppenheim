import matplotlib.pyplot as plt
import numpy as np

__author__ = 'ernesto'

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


# parámetros
M = 7

# rango de frecuencia
domega = 4 * np.pi / (M + 1)
omega_max = 2 * np.pi + domega
omega_min = -domega
nomega = 500
omega = np.linspace(omega_min, omega_max, nomega)

# repuesta en frecuencia
H = np.sin(omega * (M + 1) / 2) / np.sin(omega / 2) * np.exp(-1j * omega * M / 2)
# magnitud y fase de la respuesta en frecuencia
magH = np.abs(H)
phH = np.angle(H)

# valores maximos y minimos de los ejes
dx = 0.6
xmax_ax = omega_max + dx
xmin_ax = omega_min - dx
ymax_ax = 1.3 * (M + 1)
ymin_ax = -0.3 * (M + 1)
ymin_plot = -0.5 * (M + 1)

# length of the ticks for all subplot (6 pixels)
display_length = 7  # in pixels
# x y ticks labels margin
xtm = -0.9
xtm2 = -1.5
ytm = -0.2
# font size
fontsize = 12

xticks1 = [np.pi, 2 * np.pi]
xticks_labels1 = ['$\pi$', '$2\pi$']
xticks2 = [-2 * np.pi / (M + 1), 2 * np.pi / (M + 1)]
xticks_labels2 = [r'$-\dfrac{2\pi}{M+1}$', r'$\dfrac{2\pi}{M+1}$']

fig = plt.figure(0, figsize=(7, 4), frameon=False)

ax = plt.subplot2grid((4, 1), (0, 0), rowspan=4, colspan=1)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_plot, ymax_ax)

# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

# magnitud de la resuesta en frecuencia
plt.plot(omega, magH, 'k-', lw=2)

for i in np.arange(len(xticks1)):
    plt.plot([xticks1[i], xticks1[i]], [0, xtl], 'k-', lw=1)
    plt.text(xticks1[i], xtm, xticks_labels1[i], fontsize=fontsize, ha='center', va='baseline')
for i in np.arange(len(xticks2)):
    plt.plot([xticks2[i], xticks2[i]], [0, xtl], 'k-', lw=1)
    plt.text(xticks2[i], xtm2, xticks_labels2[i], fontsize=fontsize, ha='center', va='baseline')

# yticks y ytickslabels
plt.plot([0, ytl], [M + 1, M + 1], 'k-', lw=1)
plt.text(ytm, M + 1, '$M+1$', fontsize=fontsize, ha='right', va='center')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, r'$\left|\dfrac{\operatorname{sen}[\omega(M+1)/2]}{\operatorname{sen}(\omega/2)}\right|$',
         fontsize=fontsize, ha='left', va='center')

plt.text(xmax_ax, ymax_ax, r'$M={}$'.format(M), fontsize=fontsize, ha='right', va='center')

# ancho del lóbulo principal
ym1 = -0.3 * (M + 1)
ym2 = -0.5 * (M + 1)
i = 0
plt.plot([xticks2[i], xticks2[i]], [ym1, ym2], 'k-', lw=1)
i = 1
plt.plot([xticks2[i], xticks2[i]], [ym1, ym2], 'k-', lw=1)
plt.text(0, (ym1 + ym2) / 2, r'$\Delta\omega_m$', fontsize=fontsize, ha='center', va='center')
plt.annotate("", xytext=(0.4, (ym1 + ym2) / 2), xycoords='data', xy=(xticks2[1], (ym1 + ym2) / 2), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=4, headlength=6, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(-0.4, (ym1 + ym2) / 2), xycoords='data', xy=(xticks2[0], (ym1 + ym2) / 2), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=4, headlength=6, facecolor='black', shrink=0.002))
plt.axis('off')

# save as pdf image
plt.savefig('filter_design_windowing_magnitude_spectrum_rectangular.pdf', bbox_inches='tight')

plt.show()
