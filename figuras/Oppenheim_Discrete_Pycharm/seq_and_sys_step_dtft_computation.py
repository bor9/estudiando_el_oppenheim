import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

__author__ = 'ernesto'

# use latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preview'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# rango de n
nmax = 5
nmin = -nmax
n0 = nmax

# vector n
n = np.arange(nmin, nmax + 1)
# construcción de u[n]
u = np.zeros(n.shape)
u[n0: ] = 1
# construcción de u[-n]
u_inv = np.flip(u)
# construcción de u_e[-n] y u_o[n]
ue = (u + u_inv) / 2
uo = (u - u_inv) / 2

# rango de los ejes
ymax = 1
ymin = -0.5

delta_n = 1
nmin_ax = nmin - delta_n
nmax_ax = nmax + delta_n
delta_y = 0.5
ymax_ax = ymax + delta_y
ymin_ax = ymin - delta_y

baseline = -0.3
fontsize1 = 12
fontsize2 = 13
y_sep = 0.4


fig = plt.figure(1, figsize=(8, 4), frameon=False)
ax = plt.subplot2grid((4, 4), (0, 0), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, u, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)
plt.text(nmax_ax, baseline, '$n$', fontsize=fontsize1, ha='center', va='baseline')
i = 0
plt.text(i, baseline, '${}$'.format(i), fontsize=fontsize1, ha='center', va='baseline')
plt.text(i - y_sep, 1, '$1$', fontsize=fontsize1, ha='right', va='center')
plt.text(nmin_ax, ymax_ax, '$u[n]$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 0), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, u_inv, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)
plt.text(nmax_ax, baseline, '$n$', fontsize=fontsize1, ha='center', va='baseline')
i = 0
plt.text(i, baseline, '${}$'.format(i), fontsize=fontsize1, ha='center', va='baseline')
plt.text(i + y_sep, 1, '$1$', fontsize=fontsize1, ha='left', va='center')
plt.text(nmin_ax, ymax_ax, '$u[-n]$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (0, 2), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, ue, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)
plt.text(nmax_ax, baseline, '$n$', fontsize=fontsize1, ha='center', va='baseline')
i = 0
plt.text(i, baseline, '${}$'.format(i), fontsize=fontsize1, ha='center', va='baseline')
plt.text(i - y_sep, 1, '$1$', fontsize=fontsize1, ha='right', va='center')
i = nmin
plt.text(i - y_sep, 0.5, r'$\dfrac{1}{2}$', fontsize=fontsize1, ha='right', va='center')
plt.text(nmin_ax, ymax_ax, r'$u_e[n]=\dfrac{1}{2}(u[n]+u[-n])$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

ax = plt.subplot2grid((4, 4), (2, 2), rowspan=2, colspan=2)
plt.xlim(nmin_ax, nmax_ax)
plt.ylim(ymin_ax, ymax_ax)
plt.annotate("", xytext=(nmin_ax, 0), xycoords='data', xy=(nmax_ax, 0), textcoords='data',
             arrowprops=dict(width=0.1, headwidth=6, headlength=8, facecolor='black', shrink=0.002))
(markers, stemlines, bl) = plt.stem(n, uo, linefmt='k', markerfmt='sk', use_line_collection=True)
plt.setp(markers, markersize=4)
plt.setp(bl, visible=False)
plt.text(nmax_ax, baseline, '$n$', fontsize=fontsize1, ha='center', va='baseline')
i = 0
plt.text(i, baseline, '${}$'.format(i), fontsize=fontsize1, ha='center', va='baseline')
i = nmin
plt.text(i - y_sep, -0.5, r'$-\dfrac{1}{2}$', fontsize=fontsize1, ha='right', va='center')
i = nmax
plt.text(i + y_sep, 0.5, r'$\dfrac{1}{2}$', fontsize=fontsize1, ha='left', va='center')
plt.text(nmin_ax, ymax_ax, r'$u_o[n]=\dfrac{1}{2}(u[n]-u[-n])$', fontsize=fontsize2, ha='left', va='center')
plt.axis('off')

plt.savefig('seq_and_sys_step_dtft_computation.pdf', bbox_inches='tight')
plt.show()