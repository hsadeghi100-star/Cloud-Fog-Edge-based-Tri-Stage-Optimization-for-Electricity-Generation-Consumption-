"""
Cloud-Fog-Edge tri-stage optimization (manuscript JII-D-26-00822) - reproducibility code.

Implements Eqs. (1)-(20) and the benchmark of Section 4.6 with SciPy's HiGHS MILP solver.
  python tristage_benchmark.py case        -> Table 13 (2 regions x 3 buildings, alpha = 0.01, 0.5, 0.99)
  python tristage_benchmark.py scale 6 30  -> rows of Table 14 (alpha = 0.5, random window shifts, seed 1)
Requirements: Python >= 3.10, numpy, scipy >= 1.9 (HiGHS).
Stage 1 uses the lexicographic rule of Section 4.2: minimize C^ED of the building, then, among the
cost-optimal schedules, use the latest equally priced hours (weights -t). This makes the result
independent of the solver; the GAMS/CPLEX code returns identical values.
"""
import math, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix

T = 24
TT = range(1, T + 1)
PV = [0,0,0,0,1.5,3.7,4.1,8.2,10.9,11.7,12.8,13.9,14.3,15,10.9,10.5,4.9,4.4,3,1.2,0,0,0,0]
WT = [8.9,11.5,13.4,12.1,9.7,10.5,7,5.9,6.4,2.8,0.5,1.7,3.2,5.1,6.7,8.9,12.1,15,11.2,9.3,7.5,6.5,6.2,7.3]
C  = [10.5,10,10,10,10,10,11,12,18,20,20,20,20,20,20,20,20,20,20,16,12,11,10.5,10.5]
# name: EN1, TDmin, TDW, TDR, window(start clock h, end clock h), interruptible, PFcon, PFcom
APPL = {
 'WM': (0.5, 2, 6, 4, (9,17), True, 0.8, 0.2),
 'DW': (1.5, 2, 6, 4, (9,17), True, 0.85, 0.15),
 'SD': (2.5, 1, 3, 2, (9,17), True, 0.75, 0.25),
 'VC': (0.9, 1, 3, 2, (11,19), True, 0.7, 0.3),
 'RE': (0.2, 24, 24, 0, (0,24), False, 0.5, 0.5),
 'TV': (0.13, 3, 9, 6, (14,23), True, 0.2, 0.8),
 'CP': (0.45, 2, 7, 5, (14,23), True, 0.02, 0.98),
 'AC': (2.0, 3, 10, 7, (14,23), False, 0.05, 0.95),
 'CO': (2.1, 2, 5, 3, (13,20), False, 0.21, 0.79),
 'PEV':(3.5, 2, 5, 3, (16,24), True, 0.5, 0.5),
}
WINDOW_MODE = 'next'   # clock window [h1,h2) -> periods h1+1..h2

def window(h1, h2):
    if WINDOW_MODE == 'next':
        return list(range(h1 + 1, h2 + 1))
    return list(range(max(h1, 1), min(h2, T) + 1))

def td_target(a, alpha):
    en, tdmin, tdw, tdr, win, intr, pfc, pfm = APPL[a]
    if a == 'RE':
        return 24.0
    return alpha * tdw + (1 - alpha) * tdr          # Eq. (4)

class LP:
    """Tiny sparse model builder for scipy.optimize.milp."""
    def __init__(s):
        s.n = 0; s.lb = []; s.ub = []; s.integ = []; s.c = []
        s.rows = []; s.cols = []; s.vals = []; s.rlo = []; s.rhi = []; s.m = 0
    def var(s, lb=0.0, ub=np.inf, integer=False, cost=0.0):
        s.lb.append(lb); s.ub.append(ub); s.integ.append(1 if integer else 0); s.c.append(cost)
        s.n += 1; return s.n - 1
    def con(s, terms, lo, hi):
        for j, v in terms:
            s.rows.append(s.m); s.cols.append(j); s.vals.append(v)
        s.rlo.append(lo); s.rhi.append(hi); s.m += 1
    def solve(s, time_limit=None, gap=1e-4):
        A = coo_matrix((s.vals, (s.rows, s.cols)), shape=(s.m, s.n)).tocsr()
        opts = {'mip_rel_gap': gap, 'disp': False}
        if time_limit: opts['time_limit'] = time_limit
        t0 = time.time()
        r = milp(np.array(s.c), constraints=LinearConstraint(A, s.rlo, s.rhi),
                 integrality=np.array(s.integ), bounds=Bounds(s.lb, s.ub), options=opts)
        r.wall = time.time() - t0
        return r

def add_edge(m, apps, alpha, cost_vec=None):
    """Eqs. (2)-(6) for a set of appliances apps=[(key, appl_name, window_shift)].
    Returns dict key -> list of x variable indices (t=1..T)."""
    X = {}
    for key, a, sh in apps:
        en, tdmin, tdw, tdr, (h1, h2), intr, _, _ = APPL[a]
        win = window(h1, h2) if a == 'RE' else [t + sh for t in window(h1, h2) if 1 <= t + sh <= T]
        td = math.floor(td_target(a, alpha) + 1e-9)
        cv = cost_vec if cost_vec is not None else C
        bi = [m.var(0, 1 if t in win else 0, True) for t in TT]
        x = [m.var(0, np.inf, False, cv[t-1]) for t in TT]
        m.con([(bi[t-1], 1) for t in win], td, td)                       # Eq. (2)
        if not intr:                                                      # Eq. (3)
            for t in TT:
                terms = [(bi[t-1], td)]
                if t > 1: terms.append((bi[t-2], -td))
                terms += [(bi[k-1], 1) for k in range(t + td, T + 1)]
                m.con(terms, -np.inf, td)
        for t in TT:                                                      # Eq. (5)
            m.con([(x[t-1], 1), (bi[t-1], -en)], 0, np.inf)
        X[key] = x
    return X

def add_fog(m, demand, B, reg, par, xvars=None):
    """Eqs. (7)-(19). demand[b][t] fixed (tri-stage) or xvars[b] = list of lists of x var idx."""
    eta, delta = par['eta'], par['delta']
    p = [par['p_ratio'] * c for c in C]
    nb = len(B)
    y = [[m.var(0, np.inf, False, -C[t-1]) for t in TT] for b in B]
    s = [[m.var(par['l'], par['u']) for t in TT] for b in B]
    ch = [[m.var() for t in TT] for b in B]
    dch = [[m.var() for t in TT] for b in B]
    im = [m.var() for t in TT]
    ex = [m.var(0, np.inf, False, -p[t-1]) for t in TT]
    BI = [m.var(0, 1, True) for t in TT]
    const = 0.0
    xs = None
    if xvars is not None:   # xs = sum_a x  (Eq. 8) as variables with cost c^t in C^FL
        xs = [[m.var(0, np.inf, False, C[t-1]) for t in TT] for b in B]
        for i in range(nb):
            for t in TT:
                m.con([(xs[i][t-1], 1)] + [(xv[t-1], -1) for xv in xvars[i]], 0, 0)
    else:
        const = sum(C[t-1] * demand[i][t-1] for i in range(nb) for t in TT)
    res = par['res']   # res[t] total regional wt+pv summed over regions
    for t in TT:
        # Eq. (9)
        terms = [(y[i][t-1], 1) for i in range(nb)] + [(ex[t-1], -1)]
        rhs = 0.0
        if xs is not None: terms += [(xs[i][t-1], -1) for i in range(nb)]
        else: rhs = sum(demand[i][t-1] for i in range(nb))
        m.con(terms, rhs, rhs)
        m.con([(im[t-1], 1), (BI[t-1], -par['lim_im'])], -np.inf, 0)      # Eq. (10)
        m.con([(ex[t-1], 1), (BI[t-1], par['lim_ex'])], -np.inf, par['lim_ex'])  # Eq. (11)
        # Eq. (12)
        terms = [(y[i][t-1], 1) for i in range(nb)]
        terms += [(ch[i][t-1], delta * eta) for i in range(nb)] + [(dch[i][t-1], -delta / eta) for i in range(nb)]
        terms += [(im[t-1], -1), (ex[t-1], 1)]
        m.con(terms, res[t-1], res[t-1])
        for i in range(nb):                                               # Eqs. (13)-(14)
            terms = [(s[i][t-1], 1), (ch[i][t-1], -delta * eta), (dch[i][t-1], delta / eta)]
            if t > 1: terms.append((s[i][t-2], -1))
            m.con(terms, par['s0'] if t == 1 else 0, par['s0'] if t == 1 else 0)
        m.con([(dch[i][t-1], 1) for i in range(nb)], -np.inf, par['DCap'])  # Eq. (16)
        m.con([(ch[i][t-1], 1) for i in range(nb)], -np.inf, par['GCap'])   # Eq. (17)
    return dict(y=y, s=s, ch=ch, dch=dch, im=im, ex=ex, BI=BI, const=const, xs=xs)

def city(R, Bper, seed=None):
    rng = np.random.default_rng(seed)
    blds = [(r, b) for r in range(R) for b in range(Bper)]
    shifts = {}
    for rb in blds:
        for a in APPL:
            shifts[(rb, a)] = 0 if (seed is None or a == 'RE') else int(rng.integers(-1, 2))
    nB = len(blds)
    par = dict(eta=0.9, delta=1.0, l=0.0, u=8.0, s0=0.0, DCap=3.0 * nB, GCap=3.0 * nB,
               lim_im=10.0 / 6 * nB, lim_ex=20.0 / 6 * nB, p_ratio=0.5,
               res=[R * (PV[t] + WT[t]) * Bper / 3 for t in range(T)])
    return blds, shifts, par


# ---------------------------------------------------------------- scenarios

def edge_stage(blds, shifts, alpha, cost_vec=None):
    """Stage 1: one independent MILP per building (Edge node). Returns demand[b][t], C_ED, max node time."""
    dem, ced, tmax, xa = [], 0.0, 0.0, {}
    for rb in blds:
        m = LP()
        X = add_edge(m, [((rb, a), a, shifts[(rb, a)]) for a in APPL], alpha, cost_vec)
        r = m.solve()
        tmax = max(tmax, r.wall)
        d = np.zeros(T)
        for (k, a), x in X.items():
            v = r.x[x]; xa[(rb, a)] = v; d += v
        dem.append(d); ced += sum(C[t] * d[t] for t in range(T))
    return dem, ced, tmax, xa

def fog_stage(blds, par, dem):
    m = LP(); F = add_fog(m, dem, blds, None, par)
    r = m.solve()
    return r.fun + F['const'], r, F

def uncoordinated(blds, shifts, alpha):
    dem, ced = [], 0.0
    for rb in blds:
        d = np.zeros(T)
        for a in APPL:
            en, tdmin, tdw, tdr, (h1, h2), intr, _, _ = APPL[a]
            win = window(h1, h2) if a == 'RE' else [t + shifts[(rb, a)] for t in window(h1, h2) if 1 <= t + shifts[(rb, a)] <= T]
            td = math.floor(td_target(a, alpha) + 1e-9)
            for t in win[:td]: d[t-1] += en
        dem.append(d); ced += sum(C[t] * d[t] for t in range(T))
    return dem, ced

def centralized(blds, shifts, par, alpha, tl=None):
    m = LP()
    X = add_edge(m, [((rb, a), a, shifts[(rb, a)]) for rb in blds for a in APPL], alpha)
    xv = [[X[(rb, a)] for a in APPL] for rb in blds]
    F = add_fog(m, None, blds, None, par, xvars=xv)
    r = m.solve(time_limit=tl)
    ced = sum(C[t] * r.x[X[(rb, a)][t]] for rb in blds for a in APPL for t in range(T))
    dem = [sum(np.array(r.x[X[(rb, a)]]) for a in APPL) for rb in blds]
    return r.fun, ced, r, dem

def par_(d):
    tot = np.sum(d, axis=0); return tot.max() / tot.mean(), tot
def peak_share(d):
    tot = np.sum(d, axis=0); return 100 * sum(tot[t] for t in range(T) if C[t] == 20) / tot.sum()


# ---------------------------------------------------------------- Stage 1 with lexicographic tie-break
def edge_lex(blds, shifts, alpha, pi):
    """Stage 1 with a lexicographic tie-break: (i) min C^ED_b; (ii) among the cost-optimal schedules,
    min sum_t pi^t x^t. C^ED is unchanged by construction."""
    dem, ced, tmax = [], 0.0, 0.0
    for rb in blds:
        m = LP()
        X = add_edge(m, [((rb, a), a, shifts[(rb, a)]) for a in APPL], alpha)
        r1 = m.solve(gap=0)
        allx = [j for x in X.values() for j in x]
        m.con([(j, m.c[j]) for j in allx], -np.inf, r1.fun + 1e-6)
        for j in allx: m.c[j] = 0.0
        for k, x in X.items():
            for t in range(T): m.c[x[t]] = pi[t]
        r2 = m.solve(gap=0)
        tmax = max(tmax, r1.wall + r2.wall)
        d = np.zeros(T)
        for k, x in X.items(): d += r2.x[x]
        dem.append(d); ced += sum(C[t] * d[t] for t in range(T))
    return dem, ced, tmax


LATE = -np.arange(1, T + 1, dtype=float)

def stage3_penalty(alpha, n_buildings):
    """Eq. (20) aggregated over buildings: time term with the comfort-dependent duration of Eq. (4),
    energy term with the Stage 1 energy (integer hours) minus xs^min = TD^min * EN^1."""
    out = {}
    for a, (en, tdmin, tdw, tdr, win, intr, pfcon, pfcom) in APPL.items():
        tdc = td_target(a, alpha); td = math.floor(tdc + 1e-9)
        out[a] = n_buildings * ((tdc - tdmin) * en * pfcon + (td * en - tdmin * en) * pfcom)
    return out

def run_case():
    blds, shifts, par = city(2, 3)
    print('alpha  scenario   C_ED      C_FL       J        PAR   peak%   time(s)')
    for alpha in (0.01, 0.5, 0.99):
        t0 = time.time(); dem, ced, _ = edge_lex(blds, shifts, alpha, LATE); cfl = fog_stage(blds, par, dem)[0]; tt = time.time() - t0
        t0 = time.time(); J, cedc, rc, demc = centralized(blds, shifts, par, alpha); tc = time.time() - t0
        du, cedu = uncoordinated(blds, shifts, alpha); cflu = fog_stage(blds, par, du)[0]
        for name, (a, b, d, tm) in {'Edge-only': (ced, 0.0, dem, None), 'Fog-only': (cedu, cflu, du, None),
                                     'Tri-stage': (ced, cfl, dem, tt), 'Cloud-only': (cedc, J - cedc, demc, tc)}.items():
            print('%-5s  %-10s %9.1f %9.1f %9.1f  %.2f  %5.1f  %s' % (alpha, name, a, b, a + b, par_(d)[0], peak_share(d), '' if tm is None else '%.1f' % tm))
        pen = stage3_penalty(alpha, len(blds))
        print('       PF_CC (Eq. 20) = %.3f  [%s]' % (sum(pen.values()), ', '.join('%s %.3f' % kv for kv in pen.items())))

SIZES = {6: (2, 3), 30: (3, 10), 120: (6, 20), 300: (10, 30), 600: (12, 50), 1000: (20, 50)}

def run_scale(ns, alpha=0.5):
    for n in ns:
        R, Bp = SIZES[n]
        blds, shifts, par = city(R, Bp, seed=1)
        t0 = time.time(); dem, ced, tmax = edge_lex(blds, shifts, alpha, LATE); te = time.time() - t0
        t0 = time.time(); cfl = fog_stage(blds, par, dem)[0]; tf = time.time() - t0
        du, cedu = uncoordinated(blds, shifts, alpha)
        line = '%5d buildings: edge mean %.3f s, Fog %.1f s, tri total %.1f s, J_tri %.1f, C_ED %.1f, peak %.1f/%.1f%%, PAR %.2f/%.2f' % (
            n, te / n, tf, te + tf, ced + cfl, ced, peak_share(dem), peak_share(du), par_(dem)[0], par_(du)[0])
        if n <= 300:
            t0 = time.time(); J, cedc, rc, _ = centralized(blds, shifts, par, alpha, tl=150); tc = time.time() - t0
            line += ' | cloud %.1f s, status %d, gap %.1f%%, J_cloud %.1f, J gap %.1f%% of C_ED' % (
                tc, rc.status, 100 * (getattr(rc, 'mip_gap', 0) or 0), J, 100 * (ced + cfl - J) / ced)
        print(line, flush=True)

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'scale':
        run_scale([int(a) for a in sys.argv[2:]] or [6, 30, 120, 300, 600, 1000])
    else:
        run_case()
