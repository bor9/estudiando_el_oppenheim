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
# bordes de frecuencias y ganancia de cada banda
wi = np.array([0, 0.2, 0.35, 0.6, 1]) * np.pi
Gi = np.array([1, 0.3, 0, 0.6])
Nmb = len(Gi)

#
# Gráfica
#

wmin = 0
wmax = np.pi
dw = 0.2
wmax_ax = wmax + dw
wmin_ax = wmin - dw

dy = 0.1
ymin_ax = -dy
ymax_ax = 1 + 2 * dy

display_length = 6
fontsize = 11
lw = 2
# x y ticks labels margin
xtm = -0.09
ytm = -0.07

d_ymax = 0

fig = plt.figure(0, figsize=(5, 3), frameon=False)
### X(jw)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(wmin_ax, wmax_ax)
plt.ylim(ymin_ax, ymax_ax + d_ymax)
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
plt.annotate("", xytext=(wmin_ax, 0), xycoords='data', xy=(wmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
for i in np.arange(Nmb):
    plt.plot([wi[i], wi[i + 1]], [Gi[i], Gi[i]], 'k-', lw=lw)
    plt.plot([wi[i + 1], wi[i + 1]], [0, xtl], 'k-', lw=1)
    plt.plot([0, ytl], [Gi[i], Gi[i]], 'k-', lw=1)
plt.plot([wi[i + 1], wi[i + 1]], [0, xtl], 'k-', lw=1)
for i in np.arange(Nmb - 1):
    plt.plot([wi[i + 1], wi[i + 1]], [Gi[i], Gi[i + 1]], 'k-', lw=lw)
    plt.text(wi[i + 1], xtm, '$\omega_{}$'.format(i + 1), fontsize=fontsize, ha='center', va='baseline')
i += 1
plt.plot([wi[i + 1], wi[i + 1]], [0, Gi[i]], 'k--', lw=1)

plt.text(-0.05, xtm, '$0$', fontsize=fontsize, ha='right', va='baseline')
plt.text(np.pi, xtm, '$\pi$', fontsize=fontsize, ha='center', va='baseline')
i = 0
plt.text(ytm, Gi[i], '$G_{}$'.format(i + 1), fontsize=fontsize, ha='right', va='center')
i = 1
plt.text(ytm, Gi[i], '$G_{}$'.format(i + 1), fontsize=fontsize, ha='right', va='center')
i = 3
plt.text(ytm, Gi[i], '$G_{}$'.format(i + 1), fontsize=fontsize, ha='right', va='center')

plt.text(wmax_ax, xtm, '$\omega$', fontsize=fontsize, ha='center', va='baseline')
plt.text(0.12, ymax_ax, r'$|H_\textrm{mb}(e^{j\omega})|$', fontsize=fontsize, ha='left', va='center')
plt.text(1.5, 0.9, r'$N_\textrm{mb}=4$', fontsize=fontsize, ha='center', va='center')

plt.axis('off')

plt.savefig('filter_design_windowing_multiband.pdf', bbox_inches='tight')

plt.show()
