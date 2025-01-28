import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import numpy.polynomial.polynomial as poly

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

def filter_transformation(b, a, gnum, gden):
    """
    Transformación en frecuencia de un filtro pasabajos.
    Cálculo de los coefcientes del numerador y denominador del filtro
    transformado a partir de los coeficientes del numerador y denominador
    del filtro pasabajos prototipo y la transformación.

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


# Ganancia mínima en banda pasante (dB) y frecuencia (rad)
gp_db = -1
thetap = 0.2 * np.pi
# Ganancia máxima en banda suprimida (dB) y frecuencia (rad)
gs_db = -15
thetas = 0.3 * np.pi

# Conversión de decibeles a escala lineal
delta_p = 10 ** (gp_db / 20)
delta_s = 10 ** (gs_db / 20)

print('delta_p = {:.6f}'.format(delta_p))
print('delta_s = {:.6f}'.format(delta_s))

### Filtro Chebyshev tipo I
N, Wn = signal.cheb1ord(thetap, thetas, -gp_db, -gs_db, analog=False, fs=2 * np.pi)
# Filtro como cascada de secciones de segundo orden
sos = signal.cheby1(N, -gp_db, Wn, btype='low', analog=False, output='sos', fs=2*np.pi)
# Conversión a zpk
z, p, k = signal.sos2zpk(sos)
np.set_printoptions(precision=5, suppress=True, floatmode='fixed')
print('Ganancia de H_lp(Z): {:.6f}'.format(k))
print('Ceros de H_lp(Z): {}'.format(z))
print('Denominador de H_lp(Z): {}'.format(sos[:, 3:]))

# Transformación a filtro pasaaltos
# parámetros de la transformación
wp = 0.6 * np.pi
alpha = -np.cos((thetap + wp) / 2) / np.cos((thetap - wp) / 2)
# numerador y denominador de la tranformación
tnum = [-alpha, -1]
tden = [1, alpha]
# conversión de sos a b, a
b, a = signal.sos2tf(sos)
# transformación
bt, at = filter_transformation(b, a, tnum, tden)
# conversión del filtro transformado a z, p, k
zt, pt, kt = signal.tf2zpk(bt, at)
print('Ganancia de H(z): {:.6f}'.format(kt))
print('Ceros de H(z): {}'.format(zt))
print('Denominador de H(z): {}'.format(pt))




