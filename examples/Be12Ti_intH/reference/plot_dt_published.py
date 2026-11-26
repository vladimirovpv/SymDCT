#!/usr/bin/env python3
"""Arrhenius plot of D_a, D_c of Be12Ti-H from the published barriers (nu0 = 1e13 1/s)."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dc_numeric_published as m
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
T = np.linspace(300, 1500, 61)
Da = np.array([m.dc_numeric(m.published_rates(t), 0) for t in T]) * 1e-20
Dc = np.array([m.dc_numeric(m.published_rates(t), 2) for t in T]) * 1e-20
W = m.published_rates(300.0)
print('D_aa 300 K num/wl', m.dc_numeric(W, 0) / m.dc_wl('dc_symbolic_1_0_0_wl.m', W))
def Eeff(t1, t2, D):
    f = lambda t: np.interp(t, T, np.log(D))
    return -(f(t2) - f(t1)) / (1/(m.KB*t2) - 1/(m.KB*t1))
for lo, hi in ((300, 500), (500, 800), (800, 1500), (300, 1500)):
    print('Eeff', lo, hi, 'a %.3f c %.3f eV' % (Eeff(lo, hi, Da), Eeff(lo, hi, Dc)))
for t in (300, 500, 600, 800, 1000, 1500):
    i = int(np.argmin(abs(T - t)))
    print(t, 'Da %.3e Dc %.3e ratio %.2f' % (Da[i], Dc[i], Dc[i] / Da[i]))
fig, ax = plt.subplots(figsize=(4.2, 3.2))
ax.semilogy(1000 / T, Dc, color='#1f4e9c', lw=1.6, label=r'$D_{cc}$ ($[001]$)')
ax.semilogy(1000 / T, Da, color='#c2491d', lw=1.6, ls='--', label=r'$D_{aa}$ ($[100]$)')
ax.set_xlabel(r'$1000/T$ (K$^{-1}$)'); ax.set_ylabel(r'$D$ (m$^2$/s)')
ax.legend(frameon=False)
ax2 = ax.secondary_xaxis('top', functions=(lambda x: 1000 / x, lambda x: 1000 / x))
ax2.set_xlabel(r'$T$ (K)'); ax2.set_xticks([300, 400, 600, 1000])
fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'own_preprints', 'figures', 'be12ti_h_arrhenius.pdf'))
