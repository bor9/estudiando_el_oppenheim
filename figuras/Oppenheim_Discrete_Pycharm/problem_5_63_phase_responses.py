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
phi0 = 3.5
nd = 2
omega0 = 2 * np.pi / 3

# valores maximos y minimos de los ejes
xmax_ax = 4
ymax_ax = 12
xmin_ax = -xmax_ax
ymin_ax = -ymax_ax

# length of the ticks for all subplot (6 pixels)
display_length = 7  # in pixels
# x ticks labels margin
xtm = -1.5
ytm = -0.3
# font size
fontsize = 14
fontsize2 = 10

fig = plt.figure(0, figsize=(9, 4), frameon=False)

ax = plt.subplot2grid((1, 4), (0, 0), rowspan=1, colspan=2)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)

# horizontal and vertical ticks length
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

# resuesta en fase
plt.plot([0, np.pi], [-phi0, -phi0], 'k-', lw=2)
plt.plot([-np.pi, 0], [phi0, phi0], 'k-', lw=2)

# xticks y xtickslabels
plt.plot([np.pi, np.pi], [0, xtl], 'k-', lw=1)
plt.text(np.pi, xtm, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.plot([-np.pi, -np.pi], [0, xtl], 'k-', lw=1)
plt.text(-np.pi, xtm, '$-\pi$', fontsize=fontsize, ha='center', va='baseline')

# yticks y ytickslabels
plt.plot([0, ytl], [phi0, phi0], 'k-', lw=1)
plt.text(-ytm, phi0, r'$\phi_0$', fontsize=fontsize, ha='left', va='center')
plt.text(ytm, -phi0, r'$-\phi_0$', fontsize=fontsize, ha='right', va='center')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, r'$\angle H(e^{j\omega})$', fontsize=fontsize, ha='left', va='center')

plt.text(xmax_ax, ymax_ax, r'$(a)$', fontsize=fontsize, ha='right', va='top')

plt.axis('off')

############################

ax = plt.subplot2grid((1, 4), (0, 2), rowspan=1, colspan=2)

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
plt.plot([0, np.pi], [-phi0, -phi0-np.pi*nd], 'k-', lw=2)
plt.plot([-np.pi, 0], [phi0+np.pi*nd, phi0], 'k-', lw=2)

# xticks y xtickslabels
plt.plot([np.pi, np.pi], [0, xtl], 'k-', lw=1)
plt.text(np.pi, xtm, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
plt.plot([-np.pi, -np.pi], [0, xtl], 'k-', lw=1)
plt.text(-np.pi, xtm, '$-\pi$', fontsize=fontsize, ha='center', va='baseline')

# yticks y ytickslabels
plt.plot([0, ytl], [phi0, phi0], 'k-', lw=1)
plt.text(-ytm, phi0, r'$\phi_0$', fontsize=fontsize, ha='left', va='center')
plt.plot([0, ytl], [-phi0, -phi0], 'k-', lw=1)
plt.text(ytm, -phi0, r'$-\phi_0$', fontsize=fontsize, ha='right', va='center')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, r'$\angle H(e^{j\omega})$', fontsize=fontsize, ha='left', va='center')

# omega_0
plt.plot([omega0, omega0], [0, xtl], 'k-', lw=1)
plt.text(omega0, -xtm, '$\omega_0$', fontsize=fontsize, ha='center', va='center')
plt.plot([omega0, omega0], [0, -phi0 - omega0 * nd], 'k--', lw=1)
plt.plot([0, omega0], [-phi0 - omega0 * nd, -phi0 - omega0 * nd], 'k--', lw=1)
plt.text(ytm, -phi0 - omega0 * nd, r'$-\phi_1$', fontsize=fontsize, ha='right', va='center')

plt.text(np.pi + 0.9, -phi0 - np.pi * nd - 0.7, '${\\rm pendiente\,}=-n_d$', fontsize=fontsize, ha='right', va='top')
plt.text(xmax_ax, ymax_ax, r'$(b)$', fontsize=fontsize, ha='right', va='top')

plt.axis('off')

# save as pdf image
plt.savefig('problem_5_63_phase_responses.pdf', bbox_inches='tight')

plt.show()
