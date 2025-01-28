import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import polynomial as P

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


N = np.arange(2, 6)
x = np.linspace(-2, 3, 500)

x0 = np.where(x < -1)
x1 = np.where((x >= -1) & (x <= 1))
x2 = np.where(x > 1)

cheb_polys = np.zeros((N.shape[0], x.shape[0]))
for i, Ni in enumerate(N):
    cheb_polys[i, x0] = ((-1) ** Ni) * np.cosh(Ni * np.arccosh(-x[x0]))
    cheb_polys[i, x1] = np.cos(Ni * np.arccos(x[x1]))
    cheb_polys[i, x2] = np.cosh(Ni * np.arccosh(x[x2]))

# Ganancia de filtro Chebyshev
# amplitud del ripple
r = 0.8
# parámetro epsilon
eps = np.sqrt(1 - r ** 2) / r
idx = np.where((x >= 0))
w = x[idx]
H = 1 / np.sqrt(1 + (eps * cheb_polys[:, idx[0]]) ** 2)

# Gráficas

# Polinomios de Chebyshev
xmax = 1.2
idx = np.where((x >= -1) & (x <= xmax))
leg = ['$N={}$'.format(Ni) for Ni in N]

xmin_ax = -1
xmax_ax = xmax
ymin_ax = -1
ymax_ax = np.amax(cheb_polys[:, idx[0]])

fig = plt.figure(0, figsize=(6, 3), frameon=False)
ax = fig.add_subplot(111)
plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.plot(x[idx], (cheb_polys[:, idx[0]]).T)

plt.legend(leg, loc='upper left', frameon=False, framealpha=1)
plt.xlabel(r'$x$')
plt.ylabel(r'$V_N(x)$')

plt.savefig('continuous_filters_chebyshev_polynomials.pdf', bbox_inches='tight')

# Filtros de Chebyshev

fig = plt.figure(1, figsize=(6, 3), frameon=False)
ax = fig.add_subplot(111)
plt.xlim(w[0], w[-1])
plt.ylim(0, 1)
plt.plot(w, H.T)
# frecuencia de corte
plt.plot([1, 1], [0, r], 'k--', lw=1)
# ripple
plt.plot([0, 1], [r, r], 'k--', lw=1)

plt.legend(leg, loc='upper right', frameon=False, framealpha=1)
plt.xlabel(r'$\Omega$')
plt.ylabel(r'$|H_c(j\Omega)|$')

plt.savefig('continuous_filters_chebyshev_responses.pdf', bbox_inches='tight')

#
# Filtro de Chebyshev tipo II
#
r = 0.2 # maxima ganancia en la banda suprimida
# parámetro epsilon
eps = r / np.sqrt(1 - r ** 2)
w = np.linspace(0, 3, 500)
w0 = np.where((w >= 0) & (w <= 1))  # w^{-1} > 1
w1 = np.where(w > 1)   # w^{-1} < 1

cheb_polys_inv = np.zeros((N.shape[0], w.shape[0]))
for i, Ni in enumerate(N):
    cheb_polys_inv[i, w0] = np.cosh(Ni * np.arccosh(1 / w[w0]))
    cheb_polys_inv[i, w1] = np.cos(Ni * np.arccos(1 / w[w1]))

H_II = np.sqrt(((eps * cheb_polys_inv) ** 2) / (1 + ((eps * cheb_polys_inv) ** 2)))

fig = plt.figure(2, figsize=(6, 3), frameon=False)
ax = fig.add_subplot(111)
plt.xlim(w[0], w[-1])
plt.ylim(0, 1)
plt.plot(w, H_II.T)
# frecuencia de corte
plt.plot([1, 1], [r, 1], 'k--', lw=1)
# ripple
plt.plot([1, w[-1]], [r, r], 'k--', lw=1)

plt.legend(leg, loc='upper right', frameon=False, framealpha=1)
plt.xlabel(r'$\Omega$')
plt.ylabel(r'$|H_c(j\Omega)|$')

plt.savefig('continuous_filters_chebyshev_II_responses.pdf', bbox_inches='tight')




#
# Polos
#
N = 5
k = np.arange(2 * N)
eps = 0.2
sk = np.sin(np.pi / 2 * (1 + 2 * k) / N) * np.sinh(np.arcsinh(1 / eps) / N)
wk = np.cos(np.pi / 2 * (1 + 2 * k) / N) * np.cosh(np.arcsinh(1 / eps) / N)

# Expresión alternativa de los polos
alpha = 1 / eps + np.sqrt(1 / (eps ** 2) + 1)
sk_alt = 0.5 * (alpha ** (1 / N) - alpha ** (-1 / N)) * np.sin(np.pi / 2 * (1 + 2 * k) / N)
wk_alt = 0.5 * (alpha ** (1 / N) + alpha ** (-1 / N)) * np.cos(np.pi / 2 * (1 + 2 * k) / N)

# Polinomio de denominador de H(s)
Hs_den = [1]
for ki in k:
    Hs_den = P.polymul(Hs_den, [-sk_alt[ki] - 1j * wk_alt[ki], 1])
Hs_den = np.real(Hs_den)
Hs_den = (eps ** 2) * (2 ** (2 * (N - 1))) * Hs_den
print('Denominador de H(s)H(-s): {}'.format(Hs_den))

# Verificación del resultado comparando con el polinomio de Chebyshev
a = np.zeros(N + 1)
a[-1] = 1
cheb = np.polynomial.chebyshev.Chebyshev(a)
VN = np.polynomial.chebyshev.cheb2poly(cheb.coef)
VN = VN * ((-1j) ** np.arange(N + 1)) # esto es la evaluación en s/j
Hs_den_theo = np.real(((-1) ** N) * P.polyadd([1], (eps ** 2) * P.polymul(VN, VN)))
print('Denominador de H(s)H(-s) teórico: {}'.format(Hs_den_theo))

# valores maximos y minimos de los ejes
max_ax = 1.5
xmin = -max_ax
xmax = max_ax
ymin = -max_ax
ymax = max_ax
# axis parameters
xmin_ax = xmin
xmax_ax = xmax
ymin_ax = ymin
ymax_ax = ymax

display_length = 6

# x ticks labels margin
xtm = -0.23
ytm = -0.1
# font size
fontsize = 14
fontsize2 = 10

fig = plt.figure(3, figsize=(4, 4), frameon=False)
ax = fig.add_subplot(111)

plt.xlim(xmin_ax, xmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.gca().set_aspect('equal', adjustable='box')
xtl, ytl = convert_display_to_data_coordinates(ax.transData, length=display_length)
# axis arrows
plt.annotate("", xytext=(xmin_ax, 0), xycoords='data', xy=(xmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
plt.annotate("", xytext=(0, ymin_ax), xycoords='data', xy=(0, ymax_ax), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))

# elipse
# semieje menor
a = (alpha ** (1 / N) - alpha ** (-1 / N)) / 2
# semieje mayor
b = (alpha ** (1 / N) + alpha ** (-1 / N)) / 2
xc = np.linspace(-a, a, 50)
plt.plot(xc, b * np.sqrt(1 - (xc / a) ** 2), 'k--', lw=1)
plt.plot(xc, -b * np.sqrt(1 - (xc / a) ** 2), 'k--', lw=1)
# ceros
zeros_marker_style = dict(marker='o', linestyle='', markersize=8, markerfacecolor='w', markeredgecolor='k',
                          markeredgewidth=1.5)
# polos
polos_marker_style = dict(marker='x', linestyle='', markersize=8, markeredgecolor='k', markeredgewidth=1.5)
plt.plot(sk, wk, **polos_marker_style)

plt.annotate("", xytext=(0, 0), xycoords='data', xy=(a, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=3, headlength=8, facecolor='black', shrink=0.002))
plt.text(a / 2, 0.08, r'$a\Omega_c$', fontsize=fontsize, ha='center', va='baseline')
plt.annotate("", xytext=(0, 0), xycoords='data', xy=(0, b), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=3, headlength=8, facecolor='black', shrink=0.002))
plt.text(-0.04, b / 2 - 0.05, r'$b\Omega_c$', fontsize=fontsize, ha='right', va='center')

plt.plot([1, 1], [0, xtl], 'k-', lw=1)
plt.text(1, xtm, r'$\Omega_c$', fontsize=fontsize, ha='center', va='baseline')

# etiquetas de los ejes
plt.text(xmax_ax, xtm, r'$\textrm{Re}(s)$', fontsize=fontsize, ha='center', va='baseline')
plt.text(-ytm, ymax_ax, r'$\textrm{Im}(s)$', fontsize=fontsize, ha='left', va='center')
plt.text(xmin_ax, ymax_ax, r'$N={}$'.format(N), fontsize=fontsize, ha='left', va='center')


plt.axis('off')

# save as pdf image
plt.savefig('continuous_filters_chebyshev_poles.pdf', bbox_inches='tight')

plt.show()




