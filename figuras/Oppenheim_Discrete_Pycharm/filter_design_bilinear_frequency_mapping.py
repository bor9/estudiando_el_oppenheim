import matplotlib.pyplot as plt
import numpy as np

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


Wmax = 20
Td = 1
W = np.linspace(-Wmax, Wmax, 200)
w = 2 * np.arctan(W * Td / 2)


fontsize = 12
xmax = Wmax
xmin = -xmax
ymax = np.pi
ymin = -ymax


xmax_ax = xmax + 2
xmin_ax = -xmax_ax
ymax_ax = ymax + 1.5
ymin_ax = -ymax_ax

fig = plt.figure(0, figsize=(7, 3), frameon=False)

ax1 = plt.subplot2grid((4, 4), (0, 0), rowspan=4, colspan=4)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=5, headlength=7, facecolor='black', shrink=0.002))
plt.plot(W, w, 'k-', lw=1.5)
plt.plot([xmin, xmax], [np.pi, np.pi], 'k--', lw=1)
plt.plot([xmin, xmax], [-np.pi, -np.pi], 'k--', lw=1)
# etiquetas
plt.text(xmax_ax, -0.8, '$\Omega$', fontsize=fontsize, ha='right', va='baseline')
plt.text(1, ymax_ax - 0.3, r'$\omega=2\arctan\left(\dfrac{\Omega T_d}{2}\right)$',
         fontsize=fontsize, ha='left', va='center')
plt.text(-0.5, np.pi + 0.1, r'$\pi$', fontsize=fontsize, ha='right', va='bottom')
plt.text(-0.5, -np.pi - 0.1, r'$-\pi$', fontsize=fontsize, ha='right', va='top')

plt.axis('off')

plt.savefig('filter_design_bilinear_frequency_mapping.pdf', bbox_inches='tight')

plt.show()
