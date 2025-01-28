import numpy as np
import numpy.polynomial.polynomial as poly

def filter_transformation(b, a, gnum, gden):
    """
    Transformación en frecuencia de un filtro pasabajos.
    Cálculo de los coefcientes del numerador y denominador del filtro transformado a partir 
    de los coeficientes del numerador y denominador del filtro pasabajos prototipo y la 
    transformación.

    Parámetros
    ----------
    b, a : vectores de dimensión (N + 1, )
        Coeficientes del numerador y denominador del filtro prototipo
    gnum, gden : vectores de una dimensión de tamaño arbitrario
        Coeficientes del numerador y denominador de la transformación

    Retorna
    -------
    bt, at : vectores de una dimensión
        Coeficientes del numerador y denominador del filtro transformado
    """
    # Orden del filtro a transformar
    N = len(b) - 1
    # Lista de polinomios que multiplican a los coeficientes b[k], a[k].
    # Son de la forma (gden ** (N - k)) * (gnum ** k)
    pfactors = [1] * (N + 1)
    k = N
    for i in np.arange(N + 1):
        for j in np.arange(k):
            pfactors[i] = poly.polymul(pfactors[i], gden)
        for j in np.arange(N - k):
            pfactors[i] = poly.polymul(pfactors[i], gnum)
        k -= 1
    # Cálculo de los coeficientes del filtro transformado
    bt = 0
    at = 0
    for i in np.arange(N + 1):
        bt = poly.polyadd(bt, b[i] * pfactors[i])
        at = poly.polyadd(at, a[i] * pfactors[i])
    return bt, at
