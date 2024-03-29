import matplotlib.pyplot as plt
import numpy as np

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Diagrama de ceros y polos

# magnitud y ángulo del cero
r = 0.9
theta = 0.8 * np.pi / 2
# cero
c = r * np.exp(1j * theta)
# valor específico de omega
w0 = 0.3 * np.pi / 2

# circulo |z|=1
w = np.linspace(0, 2 * np.pi, 200)
# circulo unidad
xd = np.cos(w)
yd = np.sin(w)

# valores maximos y minimos de los ejes
max_ax = 1.2
xmin = -max_ax
xmax = max_ax
ymin = -max_ax
ymax = max_ax

# axis parameters
xmin_ax = xmin
xmax_ax = xmax
ymin_ax = ymin
ymax_ax = ymax

# length of the ticks for all subplot (6 pixels)
display_length = 6  # in pixels
# x ticks labels margin
xtm = -0.11
ytm = -0.06
# font size
fontsize = 14
fontsize2 = 10

fig = plt.figure(3, figsize=(5, 5), frameon=False)
ax = fig.add_subplot(111)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.gca().set_aspect('equal', adjustable='box')

# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

# círculo unidad
plt.plot(xd, yd, 'k', lw=2)
# ceros
zeros_marker_style = dict(marker='o', linestyle='', markersize=7, markerfacecolor='w', markeredgecolor='k',
                          markeredgewidth=1.5)
plt.plot(c.real, c.imag, **zeros_marker_style, zorder=10)
# polos
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgecolor='k', markeredgewidth=1.5)
plt.plot(0, 0, **polos_marker_style)

# vectores
plt.annotate('', xytext=(0, 0), xycoords='data', xy=(np.cos(w0), np.sin(w0)), textcoords='data',
             arrowprops=dict(arrowstyle="-|>, head_width=0.25, head_length=1", facecolor='black',
                             patchA=None, patchB=None, shrinkA=0, shrinkB=0))
plt.annotate('', xytext=(0, 0), xycoords='data', xy=(c.real, c.imag), textcoords='data',
             arrowprops=dict(arrowstyle="-|>, head_width=0.25, head_length=1", facecolor='black',
                             patchA=None, patchB=None, shrinkA=0, shrinkB=2))
plt.annotate('', xytext=(c.real, c.imag), xycoords='data', xy=(np.cos(w0), np.sin(w0)), textcoords='data',
             arrowprops=dict(arrowstyle="-|>, head_width=0.25, head_length=1", facecolor='black',
                             patchA=None, patchB=None, shrinkA=0, shrinkB=0))
# etiquetas de los vectores
k = 0.6
plt.text(k * np.cos(w0), k * np.sin(w0) + 0.05, '$v_1$', fontsize=fontsize, ha='center', va='bottom')
k = 0.5
plt.text(k * c.real + 0.04, k * c.imag, '$v_2$', fontsize=fontsize, ha='left', va='center')
plt.text((c.real + np.cos(w0)) / 2, (c.imag + np.sin(w0)) / 2 - 0.05, '$v_3$', fontsize=fontsize, ha='left', va='top')

# ángulos
nw1 = 30
# omega
w1 = np.linspace(0, w0, 30)
ra = 0.4
plt.plot(ra * np.cos(w1), ra * np.sin(w1), 'k-', lw=1)
plt.text(0.47 * np.cos(w1[int(nw1/2)]), 0.47 * np.sin(w1[int(nw1/2)]), '$\omega$', fontsize=fontsize, ha='center',
         va='center')
# theta
w1 = np.linspace(0, theta, 30)
ra = 0.2
plt.plot(ra * np.cos(w1), ra * np.sin(w1), 'k-', lw=1)
plt.text(0.26 * np.cos(w1[18]), 0.26 * np.sin(w1[18]), r'$\theta$', fontsize=fontsize, ha='center',
         va='center')
# phi
plt.plot([c.real, c.real + 0.5], [c.imag, c.imag], 'k--', lw=1)
phi = np.angle(np.exp(1j * w0) - c) + 2 * np.pi
w1 = np.linspace(0, phi, 80)
ra = 0.25
plt.plot(c.real + ra * np.cos(w1), c.imag + ra * np.sin(w1), 'k-', lw=1)
plt.text(c.real + ra * np.cos(w1[10]), c.imag + ra * np.sin(w1[10]) + 0.02, r'$\phi_3$', fontsize=fontsize, ha='left',
         va='bottom')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, '$x$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, '$y$', fontsize=fontsize, ha='left', va='center')

plt.axis('off')

# save as pdf image
plt.savefig('transform_analysis_first_order_zplane.pdf', bbox_inches='tight')

plt.show()
