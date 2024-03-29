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
M2 = 4

# rango de frecuencia
omega_max = 2 * np.pi + 0.6
omega_min = -omega_max
nomega = 500
omega = np.linspace(omega_min, omega_max, nomega)

# repuesta en frecuencia
H = np.sin(omega * (M2 + 1) / 2)/((M2 + 1) * np.sin(omega / 2)) * np.exp(-1j * omega * M2 / 2)
# magnitud y fase de la respuesta en frecuencia
magH = np.abs(H)
phH = np.angle(H)

# valores maximos y minimos de los ejes
xmax_ax = omega_max + 0.6
xmin_ax = -xmax_ax
ymax_ax = 1.3
ymin_ax = -0.4

# length of the ticks for all subplot (6 pixels)
display_length = 7  # in pixels
# x y ticks labels margin
xtm = -0.18
xtm2 = -0.25
ytm = -0.3
# font size
fontsize = 12

xticks1 = [-2 * np.pi, -np.pi, np.pi, 2 * np.pi]
xticks_labels1 = ['$-2\pi$', '$-\pi$', '$\pi$', '$2\pi$']
xticks2 = [-2 * np.pi / (M2 + 1), 2 * np.pi / (M2 + 1)]
xticks_labels2 = [r'$-\dfrac{{2\pi}}{{{:.0f}}}$'.format(M2 + 1), r'$\dfrac{{2\pi}}{{{:.0f}}}$'.format(M2 + 1)]

fig = plt.figure(0, figsize=(9, 5), frameon=False)

ax = plt.subplot2grid((4, 1), (0, 0), rowspan=2, colspan=1)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)

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
plt.plot([0, ytl], [1, 1], 'k-', lw=1)
plt.text(ytm, 1, '$1$', fontsize=fontsize, ha='left', va='bottom')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, r'$|H(e^{j\omega})|$', fontsize=fontsize, ha='left', va='center')

plt.axis('off')

############################

# x y ticks labels margin
xtm = -0.9
xtm2 = -1.3
ytm = -0.3

ymax_ax = np.pi + 1
ymin_ax = -ymax_ax

ax = plt.subplot2grid((4, 1), (2, 0), rowspan=2, colspan=1)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)

# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

# respuesta en fase
# magnitud de la resuesta en frecuencia
plt.plot(omega, phH, 'k-', lw=2)

for i in np.arange(len(xticks1)):
    plt.plot([xticks1[i], xticks1[i]], [0, xtl], 'k-', lw=1)
    if xticks1[i] < 0:
        plt.text(xticks1[i], xtm, xticks_labels1[i], fontsize=fontsize, ha='right', va='baseline')
    else:
        plt.text(xticks1[i], xtm, xticks_labels1[i], fontsize=fontsize, ha='center', va='baseline')

# yticks y ytickslabels
plt.plot([0, ytl], [np.pi, np.pi], 'k-', lw=1)
plt.text(-0.2, np.pi, '$\pi$', fontsize=fontsize, ha='right', va='center')
plt.plot([0, ytl], [-np.pi, -np.pi], 'k-', lw=1)
plt.text(-0.2, -np.pi, '$-\pi$', fontsize=fontsize, ha='right', va='center')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, r'$\angle H(e^{j\omega})$', fontsize=fontsize, ha='left', va='center')

# plt.text(np.pi + 0.9, -phi0 - np.pi * nd - 0.7, '${\\rm pendiente\,}=-n_d$', fontsize=fontsize, ha='right', va='top')
plt.text(xmax_ax, ymax_ax, r'$(b)$', fontsize=fontsize, ha='right', va='top')

plt.axis('off')

# save as pdf image
plt.savefig('example_02_16_moving_average_freq_response.pdf', bbox_inches='tight')

plt.show()
