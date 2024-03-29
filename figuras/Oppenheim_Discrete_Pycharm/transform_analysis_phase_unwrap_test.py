import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from matplotlib.collections import LineCollection

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Filtro butter pasabajos
# Orden y frecuencia de corte
ord = 5
rp = 0.1
rs = 30
wc = 0.3 # frecuencia de corte, 0 < wc < 1

b, a = signal.ellip(ord, rp, rs, wc)
delta = np.zeros(100)
delta[1] = 1
M = 30
win = np.hanning(2 * M + 1)
h = signal.lfilter(b, a, delta)
h = h[:M + 1] * win[M:]

# Respuesta en frecuencia y retardo de grupo del filtro
ntheta = 1024
# respuesta en frecuencia
w, H = signal.freqz(b, a, ntheta)
# magnitud de la respuesta en frecuencia
magH = np.abs(H)
# fase de la respuesta en frecuencia
phiH = np.unwrap(np.angle(H))

# respuesta en frecuencia
w, H = signal.freqz(h, 1, ntheta)
# magnitud de la respuesta en frecuencia
magH = np.abs(H)
# fase de la respuesta en frecuencia
phiH = np.unwrap(np.angle(H))


cols = np.linspace(0, 1, len(H))
points = np.array([H.real, H.imag]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig = plt.figure(0, figsize=(7, 6), frameon=False)

ax = plt.subplot2grid((3, 1), (0, 0), rowspan=1, colspan=1)
plt.plot(w, magH)

ax = plt.subplot2grid((3, 1), (1, 0), rowspan=1, colspan=1)
plt.plot(w, np.angle(H))

ax = plt.subplot2grid((3, 1), (2, 0), rowspan=1, colspan=1)
plt.plot(w, phiH)

fig = plt.figure(1, figsize=(7, 6), frameon=False)
ax = fig.subplots()
ax.set_aspect('equal', adjustable='box')
ax.set(xlim=(-1, 1), ylim=(-1, 1))
lc = LineCollection(segments, cmap='viridis')
#lc = LineCollection(segments, cmap='rainbow')
lc.set_array(cols)
lc.set_linewidth(1.5)
line = ax.add_collection(lc)
#plt.plot(H.real, H.imag)
#plt.plot(0, 0, 'k.')

plt.show()







