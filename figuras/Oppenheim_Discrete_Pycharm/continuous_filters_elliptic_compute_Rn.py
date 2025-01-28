import numpy as np
from numpy.polynomial import polynomial as P
from scipy import special


def compute_Rn(N, xL):
    """Cálculo de la función racional de Chebyshev Rn(x, L)
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssx rrrrrrrrrrrrrrrrrrrrrrrre                                                                                                                                                            rxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss                     ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssseeeeeeeeeeeeeeeeeeeeeeeex c
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
