import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import polynomial as P
from scipy import special
import math

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


def compute_Rn(Amax, Amin, FB, FH, w):
    ee = 10 ** (Amax / 10) - 1
    L = np.sqrt((10 ** (Amin / 10) - 1) / ee)
    print(L)
    xL = FH / FB

    # Cálculo de las integrales elípticas completas
    kkl = special.ellipk((1 / L) ** 2)
    kk1l = special.ellipk((1 - (1 / L) ** 2) ** 2)
    kk1 = special.ellipk((1 - (1 / xL) ** 2) ** 2)
    kk = special.ellipk((1 / xL) ** 2)

    N = math.ceil(kk1l * kk / (kkl * kk1))
    print(N)

    Ni = 1 if N % 2 == 0 else 0
    Nh = N // 2

    # cálculo de los polos y ceros de Rn
    xxz = np.zeros((Nh,))
    xx = np.zeros((Nh,))
    for i in np.arange(Nh):
        u = (2 * (i + 1) - Ni) * kk / N
        # cálculo de sn(u, k)
        s, _, _, _ = special.ellipj(u, 1 / xL)
        xxz[i] = s
        xx[i] = xL / s

    # cálculo de Rn
    C = np.prod((1 - xx ** 2)) / np.prod((1 - xxz ** 2))

    Rn_num = 1 if N % 2 == 0 else [0, 1]
    Rn_den = 1
    for i in np.arange(Nh):
        Rn_num = P.polymul(Rn_num, [-(xxz[i] ** 2), 0, 1])
        Rn_den = P.polymul(Rn_den, [-(xx[i] ** 2), 0, 1])

    ww = w / FB
    Rn = C * P.polyval(ww, Rn_num) / P.polyval(ww, Rn_den)
    H = 1 / np.sqrt((1 + ee * (Rn ** 2)))
    return Rn, H


def compute_Rn_order(N, Amax, FB, FH, w):
    ee = 10 ** (Amax / 10) - 1
    xL = FH / FB

    # Cálculo de las integrales elípticas completas
    kk = special.ellipk((1 / xL) ** 2)

    Ni = 1 if N % 2 == 0 else 0
    Nh = N // 2

    # cálculo de los polos y ceros de Rn
    xxz = np.zeros((Nh,))
    xx = np.zeros((Nh,))
    for i in np.arange(Nh):
        u = (2 * (i + 1) - Ni) * kk / N
        # cálculo de sn(u, k)
        s, _, _, _ = special.ellipj(u, (1 / xL) ** 2)
        xxz[i] = s
        xx[i] = xL / s

    # cálculo de Rn
    C = np.prod((1 - xx ** 2)) / np.prod((1 - xxz ** 2))

    Rn_num = 1 if N % 2 == 0 else [0, 1]
    Rn_den = 1
    for i in np.arange(Nh):
        Rn_num = P.polymul(Rn_num, [-(xxz[i] ** 2), 0, 1])
        Rn_den = P.polymul(Rn_den, [-(xx[i] ** 2), 0, 1])

    Rn = C * P.polyval(w, Rn_num) / P.polyval(w, Rn_den)
    ww = w / FB
    Rn_H = C * P.polyval(ww, Rn_num) / P.polyval(ww, Rn_den)
    H = 1 / np.sqrt((1 + ee * (Rn_H ** 2)))

    # cálculo de L
    # L_inv = (1 / xL) ** N
    # for i in np.arange(N):
    #     u = (1 + 2 * i) * kk / N
    #     s, _, _, _ = special.ellipj(u, 1 / xL)
    #     L_inv *= (s ** 2)
    L_inv = (1 / xL) ** N
    for i in np.arange(N // 2):
        u = (1 + 2 * i) * kk / N
        s, _, _, _ = special.ellipj(u, (1 / xL) ** 2)
        L_inv *= (s ** 4)
    L = 1 / L_inv
    Amin = 10 * np.log10(1 + ee * (L ** 2))
    return Rn, H, L, Amin


delta_p = 0.1
delta_s = 0.1
FB = 4
FH = 6
Amax = -20 * np.log10(1 - delta_p)
Amin = -20 * np.log10(delta_s)

ee1 = 1 / (1 - delta_p) ** 2 - 1
print("ee1 = {}".format(ee1))

w = np.linspace(0, 20, 500)

Nmax = 6
Ns = np.arange(2, Nmax + 1)
Rns = np.zeros((w.shape[0], Ns.shape[0]))
Hs = np.zeros(Rns.shape)
Ls = np.zeros(Ns.shape)
Amins = np.zeros(Ns.shape)

for i, N in enumerate(Ns):
    Rn, H, L, Amin = compute_Rn_order(N, Amax, FB, FH, w)
    Rns[:, i] = Rn
    Hs[:, i] = H
    Ls[i] = L
    Amins[i] = Amin

print(Ls)
print(Amins)

fig = plt.figure(0, figsize=(8, 4), frameon=False)
#plt.xlim(0, w[-1])
#plt.ylim(-1, 5000)
plt.plot(w, Rns)

fig = plt.figure(1, figsize=(8, 4), frameon=False)
#plt.xlim(0, w[-1])
#plt.ylim(-1, 5000)
plt.plot(w, 1 / Rns)


fig = plt.figure(2, figsize=(8, 4), frameon=False)
#plt.xlim(0, w[-1])
#plt.ylim(-1, 5000)
plt.plot(w, 20 * np.log10(np.absolute(Rns)))
plt.gca().set_prop_cycle(None)
plt.plot([FH / FB, w[-1]], [20 * np.log10(Ls), 20 * np.log10(Ls)], ls='--', lw=1)

delta_ss = 10 ** (-Amins / 20)
fig = plt.figure(3, figsize=(8, 4), frameon=False)
plt.xlim(0, w[-1])
plt.ylim(0, 1.1)
plt.plot(w, Hs)
plt.gca().set_prop_cycle(None)
plt.plot([FH, w[-1]], [delta_ss, delta_ss], ls='--', lw=1)
plt.show()
