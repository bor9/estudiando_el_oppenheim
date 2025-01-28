import numpy as np
from scipy import signal

## Especificaciones del filtro.
# Se definenen los siguientes valores:
# delta_p: la ganancia mínima en la banda pasante es 1 - delta_p
# delta_s: ganancia máxima en la banda suprimida
# wp:      límite de la banda pasante (radianes)
# ws:      límite de la banda suprimida (radianes)

## Cálculo de la máxima atenuación en la banda pasante en decibeles como número positivo
rp = -20 * np.log10(1 - delta_p)
## Cálculo de la mínima atenuación en la banda suprimida en decibeles como número positivo
rs = -20 * np.log10(delta_s)

## Diseño del filtro: cálculo del orden N y la frecuencia de corte Wn 
## para cumplir con las tolerancias.
# wp y ws:      límites de frecuencia de las bandas pasante y suprimida
#               (en radianes si el parámetro fs se especifica como 2*np.pi)
# rp, rs: atenuación máxima y mínima en la banda pasante y suprimida respectivamente 
#               en decibeles como números positivos. 
#               La ganacia máxima en la banda pasante es 0 dB.
N, Wn = signal.ellipord(wp, ws, rp, rs, analog=False, fs=2*np.pi)
## Cálculo de los coeficientes a y b del filtro IIR 
## elíptico pasabajos que cumple las tolerancias. 
b, a = signal.ellip(N, rp, rs, Wn, btype='low', analog=False, output='ba', fs=2*np.pi)
## Cálculo de la respuesta en frecuencia del filtro
w, H = signal.freqz(b, a, worN=512, fs=2*np.pi)
