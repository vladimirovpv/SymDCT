#!/usr/bin/env python3
"""Numerical diffusion coefficient of the Be12Ti-H jump network (51 sites per primitive cell).

The full master equation is solved by perturbation theory in the wave vector.
Functions
  dc_numeric(W, axis)         general directed rates W[id], id = 0..18 (no detailed balance needed)
  published_rates(T)          rates from the published sites and barriers, nu0 = 1e13 1/s
  dc_conductance_cc(c, p)     compact series-parallel form of D_cc
Run as a script, it verifies the numerics against the Mathematica closed form and prints D(T).
Usage (from the public folder)  python3 examples/Be12Ti_intH/reference/dc_numeric_published.py
"""
import json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'output')
KB = 8.617333262e-5  # eV/K

def load_network():
    d = json.load(open(os.path.join(OUT, 'jump_structure.node_link.json')))
    L = open(os.path.join(OUT, 'jump_sites.poscar')).read().split('\n')
    A = np.array([[float(x) for x in L[i].split()] for i in (2, 3, 4)]) * float(L[1])
    n = sum(int(x) for x in L[6].split())
    pos = np.array([[float(x) for x in L[8 + i].split()[:3]] for i in range(n)])
    lbl = [x['lbl'] for x in d['nodes']]
    src = np.array([e['source'] for e in d['edges']])
    tgt = np.array([e['target'] for e in d['edges']])
    jid = np.array([e['glob_jump_cnt'] for e in d['edges']])
    dl = np.array([(pos[e['target']] + np.array(e['image']) - pos[e['source']]) @ A for e in d['edges']])
    return lbl, src, tgt, jid, dl, A

LBL, SRC, TGT, JID, DL, LAT = load_network()
N = len(LBL)
CELL_C = abs(LAT[0, 2])   # lattice constant c (A)

def dc_numeric(W, axis, return_p=False):
    """D along Cartesian axis (0 = x = [100], 2 = z = [001]) for directed rates W[0..18], in A^2 * rate unit."""
    w = np.array([W[j] for j in JID], float)
    M0 = np.zeros((N, N)); b = np.zeros(N); m2 = np.zeros(N)
    np.add.at(M0, (SRC, TGT), w)
    out = np.zeros(N); np.add.at(out, SRC, w)
    M0 -= np.diag(out)
    x = DL[:, axis]
    np.add.at(b, SRC, w * x); np.add.at(m2, SRC, w * x * x)
    # stationary row vector p M0 = 0, sum p = 1
    Aeq = np.vstack([M0.T[:-1], np.ones(N)])
    rhs = np.zeros(N); rhs[-1] = 1
    p = np.linalg.solve(Aeq, rhs)
    v = p @ b
    # M0 y = -(b - v), p.y = 0
    Ay = np.vstack([M0[:-1], p]); ry = np.concatenate([-(b - v)[:-1], [0.0]])
    y = np.linalg.solve(Ay, ry)
    D = 0.5 * p @ m2 + p @ (np.zeros(N) + _m1y(w, y, axis))
    return (D, p) if return_p else D

def _m1y(w, y, axis):
    r = np.zeros(N); np.add.at(r, SRC, w * DL[:, axis] * y[TGT]); return r

# ---- published data (Bachurin, Stihl, Vladimirov, Comput. Mater. Sci. 187, 109921 (2021)) ----
ES = dict(A=0.445, B=0.654, C=0.763, D=0.823, E=0.856, F=0.910, G=1.184)
# (jump class, i, j, E_b(i->j), E_b(j->i)); equivalent-site jumps carry one barrier
PUB = [('A-C', 'A', 'C', 0.48, 0.16), ('B-C', 'C', 'B', 0.08, 0.19), ('B-B', 'B', 'B', 0.04, 0.04),
       ('B-D', 'B', 'D', 0.19, 0.02), ('D-D', 'D', 'D', 0.10, 0.10), ('D-G', 'D', 'G', 0.44, 0.07),
       ('B-G', 'G', 'B', 0.04, 0.57), ('E-G', 'G', 'E', 0.04, 0.37), ('E-E', 'E', 'E', 0.06, 0.06),
       ('E-F', 'E', 'F', 0.09, 0.03)]
# directed id of (class, first site -> second site), W ids of the network
IDS = {('A', 'C'): 0, ('C', 'A'): 16, ('E', 'G'): 1, ('G', 'E'): 6, ('E', 'F'): 2, ('F', 'E'): 18, ('E', 'E'): 3,
       ('B', 'G'): 11, ('G', 'B'): 7, ('D', 'G'): 15, ('G', 'D'): 5, ('B', 'B'): 8, ('B', 'D'): 9, ('D', 'B'): 13,
       ('B', 'C'): 10, ('C', 'B'): 17, ('D', 'D'): 14}

def consistent_barriers():
    """Saddle energy of every jump = mean of the two estimates; barriers E_sp - E_site."""
    res = {}
    for name, i, j, f, r in PUB:
        if i == j:
            res[(i, j)] = f; continue
        esp = 0.5 * ((ES[i] + f) + (ES[j] + r))
        res[(i, j)] = esp - ES[i]; res[(j, i)] = esp - ES[j]
    return res

def published_rates(T, nu0=1e13, consistent=True):
    bar = consistent_barriers()
    W = np.zeros(19)
    for k, eb in bar.items():
        W[IDS[k]] = nu0 * np.exp(-eb / (KB * T))
    return W

def site_p(T):
    n = {s: LBL.count(s) for s in 'ABCDEFG'}
    pe = {s: np.exp(-ES[s] / (KB * T)) for s in 'ABCDEFG'}
    Z = sum(n[s] * pe[s] for s in n)
    return pe, n, Z

def dc_conductance_cc(c, pe, n, a_c=CELL_C):
    """D_cc = c^2 (C_E + C_B)/Z.  c: dict of conductances cAC cBC cBD cBG cEF cEG cEE."""
    g = 1 / (1 / c['BC'] + 2 / c['AC'])
    CE = 1 / (1 / (c['EG'] + c['EE']) + 2 / c['EF'])
    CB = 4 / (1 / c['BD'] + 1 / (c['BG'] + g))
    Z = sum(n[s] * pe[s] for s in n)
    return a_c ** 2 * (CE + CB) / Z

def conductances(W, pe):
    cc = lambda i, j: pe[i] * W[IDS[(i, j)]]
    return dict(AC=cc('A', 'C'), BC=cc('B', 'C'), BD=cc('B', 'D'), BG=cc('B', 'G'), EF=cc('E', 'F'),
                EG=cc('E', 'G'), EE=cc('E', 'E'))

def load_wl(fn):
    s = open(os.path.join(OUT, fn)).read().replace('^', '**')
    s = re.sub(r'(\d)\s*\n', r'\1\n', s)
    return s

def dc_wl(fn, W, a=None):
    s = load_wl(fn)
    env = {('W%d' % i): W[i] for i in range(19)}
    env['c'] = CELL_C
    env['a'] = LAT[1, 0] * -1 * 2 if a is None else a
    return eval('(' + s + ')', {}, env)

def check():
    rng = np.random.default_rng(1)
    worst_cc = worst_cond = worst_aa = 0
    for t in range(6):
        W = np.zeros(19)
        for i in (0,1,2,3,5,6,7,8,9,10,11,13,14,15,16,17,18):
            W[i] = rng.integers(1, 30) / rng.integers(1, 7)
        W[4] = W[12] = 0.0
        dn = dc_numeric(W, 2); dw = dc_wl('dc_symbolic_0_0_1_wl.m', W)
        worst_cc = max(worst_cc, abs(dn / dw - 1))
        # aa
        dna = dc_numeric(W, 0); dwa = dc_wl('dc_symbolic_1_0_0_wl.m', W)
        worst_aa = max(worst_aa, abs(dna / dwa - 1))
    # detailed balance sets (random energies): compact form vs full network and vs Mathematica
    for t in range(6):
        global ES
        es0 = dict(ES)
        ES = {s: rng.uniform(0, 1) for s in 'ABCDEFG'}
        T = 700.0
        pe, n, Z = site_p(T)
        cs = {k: rng.uniform(0.1, 10) for k in ('AC','BC','BD','BG','EF','EG','EE','GD','BB','DD')}
        pair = {('A','C'):'AC',('B','C'):'BC',('B','D'):'BD',('B','G'):'BG',('E','F'):'EF',('E','G'):'EG',('E','E'):'EE',('D','G'):'GD',('B','B'):'BB',('D','D'):'DD'}
        W = np.zeros(19)
        for (i, j), nm in pair.items():
            W[IDS[(i, j)]] = cs[nm] / pe[i]
            if i != j: W[IDS[(j, i)]] = cs[nm] / pe[j]
        dn = dc_numeric(W, 2)
        dcomp = dc_conductance_cc(conductances(W, pe), pe, n)
        worst_cond = max(worst_cond, abs(dn / dcomp - 1))
        ES = es0
    return worst_cc, worst_aa, worst_cond

if __name__ == '__main__':
    print('check (D_cc num vs wl, D_aa num vs wl, compact form vs num):', check())
    print('max DB deviation of published barriers:')
    for name, i, j, f, r in PUB:
        if i != j: print(name, round(f - r - (ES[j] - ES[i]), 3))
    print('T, D_a, D_c, D_c/D_a  (A^2/s)')
    for T in (300, 600, 1000, 1500):
        W = published_rates(T); pe, n, Z = site_p(T)
        da, dc = dc_numeric(W, 0), dc_numeric(W, 2)
        comp = dc_conductance_cc(conductances(W, pe), pe, n)
        print(T, da, dc, dc / da, comp / dc - 1)
