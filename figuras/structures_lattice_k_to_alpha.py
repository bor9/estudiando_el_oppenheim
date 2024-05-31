import numpy as np

def k_to_alpha(k):
    """Cálculo de los coeficientes de A(z) a partir de los parámetros k.

    Parámetros
    ----------
    k : vector numpy de dimension (M, )
        vector con los parámetros k

    Retorna
    -------
    vector numpy de dimensión (M, )
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
        # en la siguiente iteración.
        alpha_iprev = alpha_i
    return alpha_i 
