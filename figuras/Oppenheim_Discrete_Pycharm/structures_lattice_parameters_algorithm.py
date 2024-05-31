import numpy as np

__author__ = 'ernesto'

def k_to_alpha(k):
    """Cálculo de los coeficientes de A(z) a partir de los parámetros k.

    Parámetros
    ----------
    k : vector numpy de dimension (M, )
        vector con los parámetros k

    Retorna
    -------
    numpy vector de dimensión (M, )
        vector con los coeficientes
    """
    M = len(k)
    for i in np.arange(M):
        # inicialización del vector \alpha^{(i)}
        alpha_i = np.zeros((i + 1, ))
        for m in np.arange(i):
            # cálculo de los elementos del vector \alpha^{(i)}[m]
            # para m = 0, 1, ..., i - 1
            alpha_i[m] = alpha_iprev[m] - k[i] * alpha_iprev[i - m - 1]
        # cálculo del último elemento del vector: \alpha^{(i)}[i]
        alpha_i[i] = k[i]
        # se almacena \alpha^{(i)}, necesario para calcular \alpha^{(i+1)}
        # en la siguiente iteracción.
        alpha_iprev = alpha_i
        print("alpha^({}) = {}".format(i + 1, alpha_i))
    return alpha_i


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
        # en la siguiente iteracción.
        alpha_i = alpha_iprev
        # inicialización del vector \alpha^{(i-1)}
        alpha_iprev = np.zeros((i, ))
        for m in np.arange(i):
            # cálculo de los elementos del vector \alpha^{(i-1)}[m]
            # para m = 0, 1, ..., i - 1
            alpha_iprev[m] = (alpha_i[m] + k[i] * alpha_i[i - m - 1]) / (1 - k[i] ** 2)
        # se asigna el valor de k, que es el último valor del vector \alpha^{(i-1)}
        k[i - 1] = alpha_iprev[i - 1]
        print("alpha^({}) = {}".format(i, alpha_iprev))
    return k


alpha_in = np.array([0.9, -0.64, 0.576])
k_in = np.array([0.673, -0.182, 0.576])

print("Algoritmo conversión k a alpha")
alpha = k_to_alpha(k_in)
print("\n")
print("Algoritmo conversión alpha a k")
k = alpha_to_k(alpha_in)

