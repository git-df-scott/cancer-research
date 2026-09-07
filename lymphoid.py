"""
Agent-based model of a lymphoid tissue compartment under T-cell-engager therapy.

WHY THIS IS NOT THE SOLID-TUMOUR MODEL (../tumor.py)
----------------------------------------------------
The parent project models a spheroid with an oxygen field, where a diffusing drug kills
cells that attempt division. Three of those assumptions are wrong for follicular lymphoma
(see UNDERSTANDING_NHL.md):
  - Oxygen: FMISO-PET shows lymphoma is markedly less hypoxic than glioblastoma
    (tumour-to-normal 1.80 vs 2.75), and germinal-centre hypoxia is physiological
    signalling in a ~100 um structure, not diffusion-limited necrosis. NO OXYGEN FIELD HERE.
  - Drug: CD20xCD3 bispecifics do not diffuse-and-poison. They bridge a T cell to a B cell.
    The killer is a discrete, motile, exhaustible AGENT. Killing is not cycle-specific.
  - Resistance: FL relapse arises from a pre-existing common-progenitor reservoir, not from
    a clone evolving under drug pressure inside the lesion. No evolving resistance trait.

THE ONE STRUCTURAL FACT THIS MODEL EXISTS TO EXPRESS
---------------------------------------------------
"Peri-follicular regions represented a barrier for immune infiltration into the follicles...
 FL-cells in follicles were ... separated spatially from the attack by CD8+ T cells than that
 in the peri-follicular regions."
     -- imaging mass cytometry of paired FL/POD24 biopsies, J Hematol Oncol 2022, PMC9396877

A T cell cannot walk through a packed follicle. It can only advance into a site that killing
has emptied. So the cytotoxic front is SELF-PACED: T cells eat inward at a rate set by their
own killing, and the deeper they get the longer they have been engaging, which is what
exhausts them. Architecture and exhaustion are therefore coupled, and that coupling is absent
from both the non-spatial QSP models (Betts 2019; Susilo 2025) and from the one existing ABM
of engager scheduling (Obertopp/Basanta bioRxiv 2025), which seeds T cells randomly among
targets at 50% occupancy -- a petri dish, i.e. a leukaemia, not a lymphoma.

UNITS
-----
One lattice site = one cell = 10 um. One step = 1 MINUTE.
These are real units, unlike the parent project's abstract "steps", because the calibration
data (T-cell speed in um/min, kills per T cell per day, exhaustion over 28 days) are all in
real time. Every rate below cites its source.

CALIBRATION
-----------
  T-cell speed        ~11 um/min in lymph node cortex (Miller/Cahalan, Science 2002)
                      -> ~1 site/min. Inside packed tissue movement is blocked, not slowed.
  CTL killing         2-16 targets per CTL per day in vivo (Halle et al., Immunity 2016,
                      PMID 26872694). T cells stay motile ("kinapses"), they do not form
                      stable clusters.
  Cooperativity       death probability rises when >2 CTLs contact one target (same paper).
                      Implemented as an optional multi-hit rule.
  Engagement cycle    engage-kill-detach within ~25 min in vivo (Cazaux, J Exp Med 2019).
  B-cell cycle        ~2 days mean human tumour cell-cycle time -> p_div = 1/2880 per min.
  Exhaustion          continuous CD19xCD3 exposure: specific lysis 88.4% (day 7) -> 8.6%
                      (day 28); a treatment-free interval restores day-14 lysis to 93.4% vs
                      34.9% continuous (Philipp et al., Blood 2022, PMID 35878001).
"""
import numpy as np
from scipy import ndimage

_K3 = np.ones((3, 3), np.float32); _K3[1, 1] = 0.0


def nbr_sum(a):
    """Moore-neighbourhood sum, one optimized pass instead of eight np.roll allocations."""
    return ndimage.convolve(np.asarray(a, np.float32), _K3, mode='constant', cval=0.0)


OFFS = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
_DY = np.array([o[0] for o in OFFS]); _DX = np.array([o[1] for o in OFFS])
EMPTY, BCELL, TCELL = 0, 1, 2
MIN_PER_DAY = 1440


def _shift(a, dy, dx):
    """Shift by (dy,dx) with zero fill. Slicing, not np.roll: no wraparound to undo."""
    out = np.zeros_like(a)
    n, m = a.shape
    ys, yd = slice(max(0, -dy), n - max(0, dy)), slice(max(0, dy), n - max(0, -dy))
    xs, xd = slice(max(0, -dx), m - max(0, dx)), slice(max(0, dx), m - max(0, -dx))
    out[yd, xd] = a[ys, xs]
    return out


class Lymphoid:
    def __init__(self, L=160, seed=0,
                 p_div=1.0 / 2880,        # B-cell division per min (2-day cycle)
                 p_death=1.0 / 14400,     # B-cell background apoptosis per min
                 t_motility=1.0,          # T-cell move attempt per min (11 um/min)
                 p_kill=0.0042,           # per-min kill hazard per T-B contact -> ~6 kills/T/day (Halle 2-16)
                 t_influx=0.002,          # T-cell entry per free border site per min
                 t_div=1.0 / 2880,        # activated T-cell division per min
                 t_death=1.0 / 20160,     # T-cell background death per min
                 exhaust_per_kill=0.0,    # exhaustion is dwell-time driven, not kill-count driven
                 recover_tau=10080.0,     # 7 d; TFI reinvigoration (Philipp 2022, Weber Science 2021)
                 exhaust_tonic=2.48e-5,   # per min in contact under engager; 28 d contact -> E=1
                 multi_hit=False,         # require cooperative hits (Halle 2016)
                 hits_needed=3, hit_decay=3000.0, dt=1.0,
                 swap_prob=0.0,
                 # ---- two-component exhaustion (see calib_philipp.py, PHASE2_PREREG.md) ----
                 exhaust_model='legacy',  # 'legacy' = single-state linear accrual (experiment L1);
                                          # 'twostate' = reversible + durable, lagged accrual,
                                          # calibrated to Philipp 2022 (5 measured points).
                 k_exh=8.0e-5,            # max accrual per min, acting on REMAINING function
                 c50_exh=11520.0,         # cumulative engaged-exposure minutes at half-max rate
                 n_exh=4.0,               # Hill exponent of the lag
                 frac_durable=0.35,       # rho: share of accrual that does NOT recover on a TFI
                 recover_tau_r=2880.0):   # reversible-component recovery time constant, minutes
        self.L = L
        self.rng = np.random.default_rng(seed)
        self.p_div, self.p_death = p_div, p_death
        self.t_motility, self.p_kill = t_motility, p_kill
        self.t_influx, self.t_div, self.t_death = t_influx, t_div, t_death
        self.exhaust_per_kill, self.recover_tau = exhaust_per_kill, recover_tau
        self.exhaust_tonic = exhaust_tonic
        self.multi_hit, self.hits_needed, self.hit_decay = multi_hit, hits_needed, hit_decay
        self.dt = float(dt)      # minutes per step; motility is sub-stepped so physics is unchanged
        # swap_prob: probability that a moving T cell exchanges places with a malignant B cell it
        # bumps into, instead of being blocked. 0 = absolute volume exclusion (T cells can only
        # enter vacancies, so the tissue must be eaten from the rim inward); 1 = occupancy does not
        # impede motility. Real lymphocytes migrate through densely packed lymphoid tissue, so 0 is
        # almost certainly too restrictive. See TRAFFICKING_PREREG.md.
        self.swap_prob = float(swap_prob)
        self.exhaust_model = exhaust_model
        self.k_exh, self.c50_exh, self.n_exh = float(k_exh), float(c50_exh), float(n_exh)
        self.frac_durable, self.recover_tau_r = float(frac_durable), float(recover_tau_r)

        self.B = np.zeros((L, L), bool)     # tumour B cells
        self.T = np.zeros((L, L), bool)     # T cells
        self.E = np.zeros((L, L))           # per-T-cell exhaustion in [0,1] (total, = Er + Ed)
        self.Er = np.zeros((L, L))          # reversible component, recovers during a TFI
        self.Ed = np.zeros((L, L))          # durable component, does not recover
        self.C = np.zeros((L, L))           # per-T-cell cumulative engaged-exposure, minutes
        self.hits = np.zeros((L, L))        # sublethal damage on B cells (multi-hit mode)
        self.drug = 0.0
        self.t = 0
        self.kills = 0
        self.cum_contact = 0
        self._n_engaged = 0
        self.history = []

    # ------------------------------------------------------------------ setup
    def seed_follicles(self, centres, radius):
        """Dense packed follicles: every site inside the radius is a B cell."""
        L = self.L
        yy, xx = np.mgrid[:L, :L]
        for (cy, cx) in centres:
            self.B[(yy - cy) ** 2 + (xx - cx) ** 2 <= radius ** 2] = True

    def seed_dispersed(self, n, occupancy=0.5, box=None):
        """Leukaemia-like control: B cells scattered at sub-unity occupancy, no architecture.
        This is the Obertopp/Basanta seeding regime and is the falsifying comparator."""
        L = self.L
        y0, y1, x0, x1 = box if box else (0, L, 0, L)
        region = np.zeros((L, L), bool); region[y0:y1, x0:x1] = True
        idx = np.flatnonzero(region)
        pick = self.rng.choice(idx, min(n, int(len(idx) * occupancy)), replace=False)
        flat = np.zeros(L * L, bool); flat[pick] = True
        self.B |= flat.reshape(L, L)

    def seed_tcells(self, n):
        """T cells start in the interfollicular space (any empty site)."""
        idx = np.flatnonzero(~self.B & ~self.T)
        pick = self.rng.choice(idx, min(n, len(idx)), replace=False)
        flat = np.zeros(self.L * self.L, bool); flat[pick] = True
        self.T |= flat.reshape(self.L, self.L)

    # ------------------------------------------------------------------- step
    def step(self, drug=None):
        if drug is not None:
            self.drug = drug
        rng = self.rng
        L = self.L
        killed_now = 0
        self._n_engaged = int(((nbr_sum(self.B) > 0) & self.T).sum()) if self.T.any() else 0

        # ---- 1. T-cell killing. A T cell adjacent to a B cell, with engager present,
        #         kills it. Kill hazard scales with drug and with remaining (1 - exhaustion).
        if self.drug > 0 and self.T.any() and self.B.any():
            # effective per-contact hazard carried by each T cell
            pot = np.where(self.T, self.p_kill * self.dt * self.drug * (1.0 - self.E), 0.0)
            # accumulate, for each B site, the hazard from all adjacent T cells
            haz = nbr_sum(pot) * self.B
            ncontact = nbr_sum(self.T)
            ncontact_B = nbr_sum(self.B)      # targets adjacent to each T cell
            if self.multi_hit:
                # cooperative killing: sublethal hits accumulate and decay (Halle/Weigelin)
                self.hits *= np.exp(-self.dt / self.hit_decay)
                newhit = self.B & (rng.random((L, L)) < 1 - np.exp(-haz))
                self.hits[newhit] += 1
                dead = self.B & (self.hits >= self.hits_needed)
            else:
                dead = self.B & (rng.random((L, L)) < 1 - np.exp(-haz))
            killed_now = int(dead.sum())
            if killed_now:
                self.B[dead] = False
                self.hits[dead] = 0
                # exhaust the T cells responsible, sharing the cost among contacts
                share = dead.astype(float) / np.maximum(ncontact, 1.0)
                credit = nbr_sum(share)
                self.E += self.exhaust_per_kill * credit * self.T
            self.kills += killed_now
            # chronic-stimulation exhaustion: any T cell in contact with a target while the
            # engager is present accrues exhaustion per unit time, whether or not it kills.
            # Calibrated so ~28 d of continuous contact drives E->1 (Philipp 2022: specific
            # lysis 88.4% at day 7 -> 8.6% at day 28 under continuous exposure).
            engaged = self.T & (ncontact_B > 0)
            if self.exhaust_model == 'twostate':
                # Accrual acts on REMAINING function and is gated by a Hill function of the cell's
                # own cumulative engaged-exposure C. The lag is required by the measured curve:
                # specific lysis is still 88.4% at day 7 (function intact) while >60% of T cells
                # already coexpress PD-1/Tim-3/LAG-3, then collapses to 34.9% by day 14.
                if self.c50_exh > 0:
                    cn = self.C ** self.n_exh
                    h = cn / (cn + self.c50_exh ** self.n_exh)
                else:
                    h = 1.0          # no lag; guard against 0/0 when C == 0 and C50 == 0
                rate = self.k_exh * h * (1.0 - self.E) * self.dt * self.drug
                self.Er += (1.0 - self.frac_durable) * rate * engaged
                self.Ed += self.frac_durable * rate * engaged
                self.C += self.dt * engaged
                np.clip(self.Ed, 0.0, 1.0, out=self.Ed)
                np.clip(self.Er, 0.0, 1.0, out=self.Er)
                self.E = np.clip(self.Er + self.Ed, 0.0, 1.0)
            elif self.exhaust_tonic:
                self.E += self.exhaust_tonic * self.dt * engaged
            np.clip(self.E, 0.0, 1.0, out=self.E)

        # ---- 2. Exhaustion recovery when the engager is absent (Philipp 2022 TFI effect)
        if self.drug == 0:
            if self.exhaust_model == 'twostate':
                # only the reversible component recovers. The durable component is required by
                # Philipp's day-28 TFI point: after a full 7-day rest the TFI arm returns to 58.7%
                # specific lysis, not to the ~93% it reached after the first rest. Recovery is
                # partial, and the shortfall grows with cumulative exposure.
                if self.recover_tau_r > 0:
                    self.Er *= np.exp(-self.dt / self.recover_tau_r)
                self.E = np.clip(self.Er + self.Ed, 0.0, 1.0)
            elif self.recover_tau > 0:
                self.E *= np.exp(-self.dt / self.recover_tau)

        # ---- 3. B-cell death and division
        occupied = self.B | self.T
        bdie = self.B & (rng.random((L, L)) < self.p_death * self.dt)
        self.B[bdie] = False
        occupied = self.B | self.T
        empty = ~occupied
        wants = self.B & (rng.random((L, L)) < self.p_div * self.dt)
        self._place(wants, empty, what='B')

        # ---- 4. T-cell division (antigen-driven expansion) and death
        if self.T.any():
            adjB = nbr_sum(self.B) > 0
            tdie = self.T & (rng.random((L, L)) < self.t_death * self.dt)
            self.T[tdie] = False
            self.E[tdie] = self.Er[tdie] = self.Ed[tdie] = self.C[tdie] = 0.0
            occupied = self.B | self.T
            empty = ~occupied
            twants = self.T & adjB & (self.drug > 0) & (rng.random((L, L)) < self.t_div * self.dt * (1.0 - self.E))
            self._place(twants, empty, what='T')

        # ---- 5. T-cell motility. This is the crux: a T cell can ONLY move into an EMPTY
        #         site, so a packed follicle is impenetrable until killing opens it.
        if self.t_motility > 0 and self.T.any():
            nsub = max(1, int(round(self.dt)))          # ~1 site per minute of simulated time
            for _ in range(nsub):
                empty = ~(self.B | self.T)
                movers = self.T if self.t_motility >= 1.0 else (self.T & (rng.random((L, L)) < self.t_motility))
                self._move(movers, empty)

        # ---- 6. T-cell influx from outside the tissue block (recruitment)
        if self.t_influx > 0 and self.drug >= 0:
            border = np.zeros((L, L), bool)
            border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
            free = border & ~self.B & ~self.T
            newT = free & (rng.random((L, L)) < self.t_influx * self.dt)
            self.T[newT] = True
            self.E[newT] = self.Er[newT] = self.Ed[newT] = self.C[newT] = 0.0

        self.t += self.dt
        self._record(killed_now)

    # --------------------------------------------------------------- helpers
    def _place(self, wants, empty, what):
        """Index-based: each dividing cell picks one random neighbour direction; collisions on
        the same target site are resolved by keeping one parent at random."""
        ys, xs = np.nonzero(wants)
        if not len(ys):
            return
        rng = self.rng
        L = self.L
        k = rng.integers(0, 8, size=len(ys))
        ty = ys + _DY[k]; tx = xs + _DX[k]
        ok = (ty >= 0) & (ty < L) & (tx >= 0) & (tx < L)
        ys, xs, ty, tx = ys[ok], xs[ok], ty[ok], tx[ok]
        if not len(ys):
            return
        ok = empty[ty, tx]
        ys, xs, ty, tx = ys[ok], xs[ok], ty[ok], tx[ok]
        if not len(ys):
            return
        flat = ty * L + tx                      # one winner per contested target site
        order = rng.permutation(len(flat))
        _, first = np.unique(flat[order], return_index=True)
        w = order[first]
        ys, xs, ty, tx = ys[w], xs[w], ty[w], tx[w]
        if what == 'B':
            self.B[ty, tx] = True
        else:
            self.T[ty, tx] = True
            for a in (self.E, self.Er, self.Ed, self.C):
                a[ty, tx] = a[ys, xs]

    def _move(self, movers, empty):
        """T-cell random walk. A T cell can only step into an EMPTY site, so a packed follicle
        is impenetrable until killing opens it. This is the model's central rule.
        Index-based over the (sparse) T-cell population."""
        ys, xs = np.nonzero(movers)
        if not len(ys):
            return
        rng = self.rng
        L = self.L
        k = rng.integers(0, 8, size=len(ys))
        ty = ys + _DY[k]; tx = xs + _DX[k]
        ok = (ty >= 0) & (ty < L) & (tx >= 0) & (tx < L)
        ys, xs, ty, tx = ys[ok], xs[ok], ty[ok], tx[ok]
        if not len(ys):
            return
        free = empty[ty, tx]
        if self.swap_prob > 0:
            # bump into a malignant B cell and squeeze past it, exchanging places
            swap = self.B[ty, tx] & (rng.random(len(ty)) < self.swap_prob)
        else:
            swap = np.zeros(len(ty), bool)
        ok = free | swap
        ys, xs, ty, tx, swap = ys[ok], xs[ok], ty[ok], tx[ok], swap[ok]
        if not len(ys):
            return
        flat = ty * L + tx                      # one winner per contested target site
        order = rng.permutation(len(flat))
        _, first = np.unique(flat[order], return_index=True)
        w = order[first]
        ys, xs, ty, tx, swap = ys[w], xs[w], ty[w], tx[w], swap[w]
        # a swapping T cell pushes a B cell back into its old site, so that site must not also be
        # somebody else's destination; drop those swaps rather than double-occupy
        if swap.any():
            tgt_set = set((ty * L + tx).tolist())
            clash = np.array([(swap[i] and (ys[i] * L + xs[i]) in tgt_set) for i in range(len(ys))])
            keep = ~clash
            ys, xs, ty, tx, swap = ys[keep], xs[keep], ty[keep], tx[keep], swap[keep]
        if not len(ys):
            return
        carried = [(a, a[ys, xs].copy()) for a in (self.E, self.Er, self.Ed, self.C)]
        sy, sx = ys[swap], xs[swap]             # B cells displaced backwards
        by, bx = ty[swap], tx[swap]
        self.T[ys, xs] = False
        for a, _ in carried:
            a[ys, xs] = 0.0
        self.B[by, bx] = False
        self.T[ty, tx] = True
        for a, v in carried:
            a[ty, tx] = v
        self.B[sy, sx] = True

    def _record(self, killed):
        nT = int(self.T.sum())
        nEng = int(self._n_engaged) if hasattr(self, '_n_engaged') else 0
        rec = dict(t=self.t, nB=int(self.B.sum()), nT=nT, drug=self.drug, killed=killed,
                   meanE=float(self.E[self.T].mean()) if nT else 0.0,
                   engaged=nEng, engFrac=(nEng / nT) if nT else 0.0)
        self.history.append(rec)
        self.cum_contact += nEng * self.dt   # total T-cell-minutes spent in contact

    # ------------------------------------------------------------- measurement
    def radial_profile(self, centre, rmax, nbins=12):
        """T-cell occupancy and mean exhaustion vs distance from a follicle centre.
        This is the quantity the imaging-mass-cytometry papers measure."""
        L = self.L
        yy, xx = np.mgrid[:L, :L]
        r = np.hypot(yy - centre[0], xx - centre[1])
        edges = np.linspace(0, rmax, nbins + 1)
        out = []
        for i in range(nbins):
            m = (r >= edges[i]) & (r < edges[i + 1])
            nb, nt = int((self.B & m).sum()), int((self.T & m).sum())
            tot = nb + nt
            out.append(dict(r_lo=float(edges[i]), r_hi=float(edges[i + 1]),
                            nB=nb, nT=nt,
                            t_frac=(nt / tot) if tot else np.nan,
                            meanE=float(self.E[self.T & m].mean()) if nt else np.nan))
        return out

    @property
    def nB(self):
        return int(self.B.sum())

    def run(self, steps, schedule=None, record_every=60, stop_when_clear=True):
        for _ in range(steps):
            d = schedule(self) if schedule else None
            self.step(d)
            if stop_when_clear and self.nB == 0:
                break
        return self
