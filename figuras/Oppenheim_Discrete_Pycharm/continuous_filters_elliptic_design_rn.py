import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import polynomial as P
from scipy import special


__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


def compute_Rn(N, xL):
    """Cálculo de la función racional de Chebyshev Rn(x, L)

    Parámetros
    ----------
    N : entero
        Orden de  Rn
    xL : real mayor a 1
        Primer valor de x que se cumple que Rn(x, L) = L

    Retorna
    -------
    Rn_num :
        vector con los coeficientes de numerador de Rn
    Rn_den :
        vector con los coeficientes de denominador de Rn
    L :
        real tal que Rn(x, L) > L en x > XL
    """
    # Cálculo de la integral elíptica completa con módulo 1/xL, K(1/xL)
    K = special.ellipk((1 / xL) ** 2)
    # Cálculo de los ceros y polos de Rn
    Ni = 1 if N % 2 == 0 else 0
    Nh = N // 2
    xzi = np.zeros((Nh, ))
    xi = np.zeros((Nh, ))
    for i in np.arange(Nh):
        # Cálculo del seno elíptico de Jacobi sn(u, k)
        u = (2 * (i + 1) - Ni) * K / N
        sn, _, _, _ = special.ellipj(u, (1 / xL) ** 2)
        xzi[i] = sn
        xi[i] = xL / sn
    # Cálculo de la constante multiplicativa C
    C = np.prod((1 - xi ** 2)) / np.prod((1 - xzi ** 2))
    # Cálculo de los coeficientes de los polinomios del numerador y denominador de Rn
    Rn_num = C if N % 2 == 0 else [0, C]
    Rn_den = 1
    for i in np.arange(Nh):
        Rn_num = P.polymul(Rn_num, [-(xzi[i] ** 2), 0, 1])
        Rn_den = P.polymul(Rn_den, [-(xi[i] ** 2), 0, 1])
    # Cálculo de L, que es función de xL y N
    L_inv = (1 / xL) ** N
    for i in np.arange(Nh):
        u = (1 + 2 * i) * K / N
        sn, _, _, _ = special.ellipj(u, (1 / xL) ** 2)
        L_inv *= (sn ** 4)
    L = 1 / L_inv
    return Rn_num, Rn_den, L


# Cálculo para graficas de Rn
xL = 1.5
x = np.linspace(0, 7, 500)
Nmax = 6
Ns = np.arange(2, Nmax + 1)
Rns = np.zeros((x.shape[0], Ns.shape[0]))
Ls = np.zeros(Ns.shape)

for i, N in enumerate(Ns):
    Rn_num, Rn_den, L = compute_Rn(N, xL)
    Rn = P.polyval(x, Rn_num) / P.polyval(x, Rn_den)
    Rns[:, i] = Rn
    Ls[i] = L

# Cálculo de Hc(jw) a partir de Rn
# Parámetros
delta_p = 0.1
wB = 4

wH = xL * wB
Amax = -20 * np.log10(1 - delta_p)
# epsilon
ee = 10 ** (Amax / 10) - 1
w = np.linspace(0, 20, 500)
Hs = np.zeros((w.shape[0], Ns.shape[0]))
Amins = np.zeros(Ns.shape)

for i, N in enumerate(Ns):
    Rn_num, Rn_den, L = compute_Rn(N, xL)
    Rn = P.polyval(w / wB, Rn_num) / P.polyval(w / wB, Rn_den)
    Hs[:, i] = 1 / np.sqrt((1 + ee * (Rn ** 2)))
    Amins[i] = 10 * np.log10(1 + ee * (L ** 2))


gg = 0.4
grey = [gg, gg, gg]
# Gráfica de Rn
leg = [r'$N={}\;(L\approx{:.2f})$'.format(Ns[i], Ls[i]) for i in np.arange(Ns.shape[0])]
ymin = 10e-4

fig = plt.figure(0, figsize=(8, 4), frameon=False)
plt.semilogy(x, np.absolute(Rns))
plt.gca().set_prop_cycle(None)
plt.semilogy([xL, x[-1]], [Ls, Ls], ls='--', lw=1)
plt.plot([1, 1], [10e-4, 1], color=grey, ls='--', lw=1)
plt.plot([xL, xL], [10e-4, Ls[-1]], color=grey, ls='--', lw=1)
# se agregan ticks extra
plt.xticks(list(plt.xticks()[0]) + [xL])
plt.yticks(list(plt.yticks()[0]) + [1])
plt.xlim(0, x[-1])
plt.ylim(ymin, 10e6)

plt.legend(leg, loc='lower right', frameon=False, framealpha=1)
plt.xlabel(r'$x$')
plt.ylabel(r'$|R_N(x,\,L)|$')

# save as pdf image
plt.savefig('continuous_filters_elliptic_design_Rn.pdf', bbox_inches='tight')


# Gráfica de H
delta_ss = 10 ** (-Amins / 20)
leg = [r'$N={}\;(A_\textrm{{min}}\approx{:.2f}\textrm{{ dB o }}\delta_s\approx{:.4f})$'.
           format(Ns[i], Amins[i], delta_ss[i]) for i in np.arange(Ns.shape[0])]
ymin = 0
ymax = 1
wmin = 0
wmax = w[-1]

fig = plt.figure(1, figsize=(8, 4), frameon=False)
plt.xlim(wmin, wmax)
plt.ylim(ymin, ymax)
plt.plot(w, Hs)
plt.gca().set_prop_cycle(None)
plt.plot([wH, w[-1]], [delta_ss, delta_ss], ls='--', lw=1)
plt.plot([wB, wB], [0, 1 - delta_p], color=grey, ls='--', lw=1)
plt.plot([wH, wH], [0, delta_ss[0]], color=grey, ls='--', lw=1)
plt.plot([0, wB], [1 - delta_p, 1 - delta_p], color=grey, ls='--', lw=1)
# se agregan ticks extra
plt.xticks(list(plt.xticks()[0]) + [wB, wH])
plt.yticks(list(plt.yticks()[0]) + [1 - delta_p])

plt.legend(leg, loc='upper right', frameon=False, framealpha=1)
plt.xlabel(r'$\Omega\,(\textrm{rad/s})$')
plt.ylabel(r'$H_c(j\Omega)$')

# save as pdf image
plt.savefig('continuous_filters_elliptic_design_H.pdf', bbox_inches='tight')


plt.show()
