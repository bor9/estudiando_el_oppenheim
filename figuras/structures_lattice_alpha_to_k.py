import numpy as np

def alpha_to_k(alpha):
    """Cálculo de los parámetros k a partir de los coeficientes de A(z).

    Parámetros
    ----------
    alpha : vector numpy de dimension (M, )
        vector con los coeficientes de A(z)

    Retorna
    -------
    numpy vector de dimensión (M, )
        vector con los parámetros k
    """
    M = len(alpha)
    # inicialización del vector con los parámetros k
    k = np.zeros((M, ))
    k[M - 1] = alpha[M - 1]
    alpha_iprev = alpha
    for i in np.arange(M - 1, 0, -1):
        # se almacena \alpha^{(i)}, necesario para calcular \alpha^{(i-1)}
        # en la siguiente iteración.
        alpha_i = alpha_iprev
        # inicialización del vector \alpha^{(i-1)}
        alpha_iprev = np.zeros((i, ))
        for m in np.arange(i):
            # cálculo de los elementos del vector \alpha^{(i-1)}[m]
            # para m = 0, 1, ..., i - 1
            alpha_iprev[m] = (alpha_i[m] + k[i] * alpha_i[i - m - 1]) / (1 - k[i] ** 2)
        # se asigna el valor de k, que es el último valor del vector \alpha^{(i-1)}
        k[i - 1] = alpha_iprev[i - 1]
    return k
