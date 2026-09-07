"""
Multi-follicle tumour construction for experiment 8.

A lymphoma is not one compact spheroid. Follicular lymphoma grows as many small follicles inside a
richly vascularised lymph node, each too small to develop a hypoxic core. This module grows M
separate nodules on one lattice (M = 1 is the compact mass used in every earlier experiment) to the
same total size N0, so that geometry is the only thing that changes.

Nodules are seeded on a sqrt(M) x sqrt(M) grid centred in the box, spaced so that at N0 they do not
touch. L = 160 for every M (including M = 1) so the wall stays >= 35 sites from every nodule centre.
"""
import numpy as np
from tumor import Tumor
import pilot

L = 160


def centres(M):
    k = int(round(np.sqrt(M)))
    assert k * k == M
    if M == 1:
        return [(L // 2, L // 2)]
    spacing = {4: 48, 9: 40, 16: 30}[M]
    off = (k - 1) * spacing / 2
    return [(int(round(L / 2 - off + i * spacing)), int(round(L / 2 - off + j * spacing))) for i in range(k) for j in range(k)]


def make_follicular(seed, M, params, init_radius=3):
    p = dict(params, L=L, init_resist_n=0)
    T = Tumor(seed=seed, **p)
    T.alive[:] = False; T.prolif[:] = 0; T.resist[:] = 0; T.gen[:] = 0
    yy, xx = np.mgrid[:L, :L]
    for (cy, cx) in centres(M):
        disc = (yy - cy) ** 2 + (xx - cx) ** 2 <= init_radius ** 2
        T.alive[disc] = True; T.prolif[disc] = T.base_prolif
    T._relax_o2()
    T.history = []
    while T.n < pilot.N0:
        T.step(0.0)
        if T.t > 6000:
            raise RuntimeError('did not reach N0')
    return T


def follicle_stats(T, M):
    """Per-nodule radius (from area) and merge check via connected components."""
    from scipy import ndimage
    lab, ncomp = ndimage.label(T.alive, structure=np.ones((3, 3)))
    sizes = ndimage.sum(T.alive, lab, range(1, ncomp + 1))
    r = np.sqrt(np.array(sizes) / np.pi)
    hyp = float((T.alive & (T.o2 < T.hypoxia)).sum() / T.n)
    return dict(ncomp=int(ncomp), r_med=float(np.median(r)), r_min=float(r.min()), r_max=float(r.max()), hypoxic=hyp)


def implant_centre(T, M, K=4, seed=0):
    """Compact K-cell resistant clone at the centre of the nodule nearest the box centre."""
    cs = centres(M)
    cy, cx = min(cs, key=lambda c: (c[0] - L / 2) ** 2 + (c[1] - L / 2) ** 2)
    # actual centroid of that nodule (it may have drifted)
    from scipy import ndimage
    lab, _ = ndimage.label(T.alive, structure=np.ones((3, 3)))
    lb = lab[cy, cx]
    if lb == 0:
        ys0, xs0 = np.nonzero(T.alive); j = np.argmin((ys0 - cy) ** 2 + (xs0 - cx) ** 2); lb = lab[ys0[j], xs0[j]]
    ys0, xs0 = np.nonzero(lab == lb)
    cy, cx = ys0.mean(), xs0.mean()
    d = (ys0 - cy) ** 2 + (xs0 - cx) ** 2
    pick = np.argsort(d)[:K]
    T.resist[ys0[pick], xs0[pick]] = 0.9
    T._relax_o2()
    return ys0[pick], xs0[pick], float(np.sqrt(len(ys0) / np.pi))


if __name__ == '__main__':
    import time
    for M in (1, 4, 9, 16):
        for a in (0.0005, 0.01):
            t0 = time.time()
            T = make_follicular(0, M, dict(pilot.BASE, mut_rate=0.0, apoptosis=a))
            s = follicle_stats(T, M)
            ys, xs, rf = implant_centre(T, M)
            q = int((T.o2[ys, xs] < T.hypoxia).sum())
            print(f"M={M:2d} a={a}: t_pre={T.t:5d} N={T.n} comps={s['ncomp']:2d} r_med={s['r_med']:5.1f} [{s['r_min']:.1f},{s['r_max']:.1f}] hypoxic={s['hypoxic']:.3f} clone quiescent={q}/4 r_f(implant)={rf:.1f}  {time.time()-t0:.0f}s ({(time.time()-t0)/T.t*1000:.1f} ms/step)")
