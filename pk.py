"""
Two-compartment pharmacokinetics for a half-life-extended CD3 bispecific.

WHY THIS EXISTS
---------------
`lymphoid.py` treats the engager as a binary: `self.drug` is 0 or 1, and every schedule is an
on/off square wave. That is a defensible approximation for blinatumomab, whose terminal half-life
is ~2 hours and which is given by continuous infusion -- stopping the pump really does take
concentration to zero within hours.

It is wrong for tarlatamab, and the wrongness is the entire point of the lung question.

    Terminal half-life  5.8 days   (DeLLphi-300 mean, PMID 39589690)
                       11.2 days   (population PK median, 420 patients, PMID 40261494)
    Clearance          0.649 L/day
    Central volume     3.44 L      (typical 73 kg subject)
    Clinical schedule  1 mg C1D1 step-up, 10 mg C1D8, 10 mg C1D15, then 10 mg every 2 weeks

An agent with an 11-day half-life dosed every 14 days does not produce a treatment-free interval.
Concentration never approaches zero; the profile is a sawtooth riding on a high trough. If T-cell
exhaustion is driven by cumulative engaged dwell time, then the half-life extension that made
tarlatamab clinically dosable may have simultaneously removed the disengagement window that the
entire treatment-free-interval literature assumes exists.

That is a falsifiable claim about a real dosing regimen, and it cannot be asked of a model whose
drug term is a switch. Hence this module.

WHAT IS PUBLISHED AND WHAT IS NOT
---------------------------------
Published and used directly:  CL = 0.649 L/day, V1 = 3.44 L, terminal half-life 5.8-11.2 days.
NOT published in what could be retrieved:  intercompartmental clearance Q, peripheral volume V2.

Two unknowns, one constraint (the terminal half-life), so the system is underdetermined. Rather
than invent values, `solve_two_compartment` fixes the peripheral-to-central volume ratio as an
explicit free parameter and solves Q analytically for the target terminal half-life. Experiments
sweep both the ratio and the half-life across the published 5.8-11.2 day range, and report the
span. Nothing here is fitted to any outcome.

UNITS
-----
Time in days externally; `occupancy_at` accepts model-minutes for use as a `Lymphoid` schedule.
Amounts in mg, volumes in L, so concentrations are mg/L.
"""
import numpy as np

LN2 = np.log(2.0)

# Published population-PK point estimates (PMID 40261494), typical 73 kg subject.
CL_DEFAULT = 0.649      # L/day
V1_DEFAULT = 3.44       # L
THALF_POPPK = 11.2      # days, population-PK median terminal half-life
THALF_D300 = 5.8        # days, DeLLphi-300 mean terminal half-life (PMID 39589690)

# Clinical regimen: step-up priming then Q2W maintenance.
TARLATAMAB_DOSES = [(0.0, 1.0), (7.0, 10.0), (14.0, 10.0)]   # (day, mg); Q2W continues after



def feasible_v_ratio_max(t_half_terminal, CL=CL_DEFAULT, V1=V1_DEFAULT):
    """Largest physically admissible V2/V1 for a given terminal half-life.

    A two-compartment model can only produce a terminal half-life slower than its own elimination
    rate, and the slower the terminal phase, the more peripheral volume it can support. Concretely
    the solution requires k10/beta > 1 + v_ratio, so

        v_ratio < CL/V1 * t_half_terminal / ln2 - 1

    This is a real constraint, not a numerical one, and it bites: at the DeLLphi-300 half-life of
    5.8 days the peripheral compartment must be under ~0.58x central, whereas at the population-PK
    median of 11.2 days it may be up to ~2.05x. Experiments must sweep inside the feasible range
    for each half-life rather than holding v_ratio fixed across them.
    """
    return (CL / V1) * t_half_terminal / LN2 - 1.0


def default_v_ratio(t_half_terminal, CL=CL_DEFAULT, V1=V1_DEFAULT, frac=0.5):
    """Midpoint of the feasible range. A stated convention, not a published value."""
    return frac * feasible_v_ratio_max(t_half_terminal, CL, V1)


def solve_two_compartment(t_half_terminal, CL=CL_DEFAULT, V1=V1_DEFAULT, v_ratio=None):
    """Return (k10, k12, k21) reproducing a target terminal half-life.

    v_ratio = V2/V1 is an explicit assumption, not a published value. It is swept, not chosen.

    Derivation. For a two-compartment IV model,
        alpha + beta = k10 + k12 + k21          alpha * beta = k10 * k21
    With k12 = (V2/V1) * k21 = v_ratio * k21 and beta = ln2 / t_half_terminal, eliminating alpha
    gives a closed form for k21:
        k21 = (k10 - beta) / (k10/beta - 1 - v_ratio)
    """
    if v_ratio is None:
        v_ratio = default_v_ratio(t_half_terminal, CL, V1)
    k10 = CL / V1
    beta = LN2 / t_half_terminal
    denom = k10 / beta - 1.0 - v_ratio
    if beta >= k10 or denom <= 0:
        raise ValueError(
            f'no physical solution: terminal half-life {t_half_terminal} d with CL/V1={k10:.4f}/day '
            f'and v_ratio={v_ratio}. Need beta < k10 and k10/beta > 1 + v_ratio.')
    k21 = (k10 - beta) / denom
    k12 = v_ratio * k21
    return k10, k12, k21


class TwoCompartmentPK:
    """IV-bolus two-compartment model, integrated analytically per inter-dose interval.

    Doses superpose linearly (the published model is linear elimination), so the concentration
    after a dose history is the sum of single-dose responses. That is exact here, not an
    approximation, and it avoids any ODE-solver tolerance question.
    """

    def __init__(self, t_half_terminal=THALF_POPPK, CL=CL_DEFAULT, V1=V1_DEFAULT, v_ratio=None):
        self.CL, self.V1 = CL, V1
        self.t_half_terminal = t_half_terminal
        self.v_ratio = default_v_ratio(t_half_terminal, CL, V1) if v_ratio is None else v_ratio
        self.k10, self.k12, self.k21 = solve_two_compartment(t_half_terminal, CL, V1, self.v_ratio)
        s = self.k10 + self.k12 + self.k21
        disc = np.sqrt(max(s * s - 4.0 * self.k10 * self.k21, 0.0))
        self.alpha = 0.5 * (s + disc)
        self.beta = 0.5 * (s - disc)
        self.doses = []

    # ------------------------------------------------------------------ dosing
    def add_dose(self, day, mg):
        self.doses.append((float(day), float(mg)))
        return self

    def add_regimen(self, doses=TARLATAMAB_DOSES, q_days=14.0, until_day=84.0, maint_mg=10.0):
        """Step-up priming then maintenance every q_days until `until_day`."""
        for d, mg in doses:
            self.add_dose(d, mg)
        t = doses[-1][0] + q_days
        while t < until_day:
            self.add_dose(t, maint_mg)
            t += q_days
        return self

    def scale_doses(self, factor):
        """Multiply every dose by `factor`. Elimination is linear, so AUC scales identically."""
        self.doses = [(t, mg * factor) for t, mg in self.doses]
        return self

    def auc(self, t0, t1, n=200001):
        """Trapezoidal AUC of the central compartment over a FINITE window.

        Exposure comparisons must be made on the window actually simulated. Equal eventual dose
        does not imply equal AUC over a finite horizon: a bolus late in the window contributes
        almost nothing to it, and a terminal tail extending past the horizon is never realised.
        """
        t = np.linspace(t0, t1, n)
        return float(np.trapezoid(self.conc(t), t))

    def add_infusion_auc_matched(self, target_auc, start=0.0, end=84.0, n_steps=2000):
        """Infusion delivering the same AUC over [start, end] as `target_auc`.

        Matching TOTAL DOSE is the wrong comparator, and measurably so. Codex's review integrated
        the two arms over days 0-84 at the 11.2-day half-life and found Q2W gave AUC 87.278 against
        the dose-matched infusion's 95.522, a 9% exposure advantage to the infusion -- because the
        Q2W history contained a bolus at exactly day 84 that cannot influence the simulated
        trajectory, yet its milligrams were spread across the whole infusion. Any timing conclusion
        drawn against that comparator could be a dose effect.

        Elimination is linear, so AUC is proportional to infused dose: scale a unit infusion.
        """
        probe = TwoCompartmentPK(self.t_half_terminal, self.CL, self.V1, self.v_ratio)
        probe.add_infusion_matched(1.0, start, end, n_steps)
        per_mg = probe.auc(start, end)
        self.add_infusion_matched(target_auc / per_mg, start, end, n_steps)
        return self

    def add_infusion_matched(self, total_mg, start=0.0, end=84.0, n_steps=2000):
        """Continuous infusion delivering the same total dose, as many small boluses.

        This is the matched-AUC comparator. It is a discretisation of an infusion, and n_steps is
        large enough that the sawtooth is negligible relative to the Q2W profile being compared.
        """
        per = total_mg / n_steps
        for t in np.linspace(start, end, n_steps, endpoint=False):
            self.add_dose(t, per)
        return self

    # ----------------------------------------------------------- concentration
    def conc(self, day):
        """Central-compartment concentration (mg/L) at `day`, summed over the dose history."""
        day = np.asarray(day, float)
        out = np.zeros_like(day)
        A = (self.alpha - self.k21) / (self.alpha - self.beta)
        B = (self.k21 - self.beta) / (self.alpha - self.beta)
        for t0, mg in self.doses:
            dt = day - t0
            m = dt >= 0
            if not np.any(m):
                continue
            c0 = mg / self.V1
            out = out + np.where(
                m, c0 * (A * np.exp(-self.alpha * np.where(m, dt, 0.0))
                         + B * np.exp(-self.beta * np.where(m, dt, 0.0))), 0.0)
        return out

    def total_mg(self):
        return sum(mg for _, mg in self.doses)

    # -------------------------------------------------------------- occupancy
    def occupancy(self, day, ec50):
        """Emax target occupancy in [0,1]. This is what replaces Lymphoid's binary `drug`."""
        c = self.conc(day)
        return c / (c + ec50)

    def as_schedule(self, ec50, minutes_per_day=1440.0):
        """Adapter for `Lymphoid.run(schedule=...)`, which passes the model and reads `s.t` in min."""
        def sched(s):
            return float(self.occupancy(s.t / minutes_per_day, ec50))
        return sched


# --------------------------------------------------------------------- checks
def _selftest():
    """Verify the PK in isolation BEFORE it is wired into the ABM.

    A silent PK bug would invalidate every downstream conclusion, and it would not announce
    itself: the ABM would simply receive a plausible-looking drug trace. So these run against the
    published quantities, not against internal consistency alone.
    """
    ok = True

    def check(name, got, want, tol):
        nonlocal ok
        good = abs(got - want) <= tol
        ok &= good
        print(f'  [{"PASS" if good else "FAIL"}] {name}: got {got:.4f}, want {want:.4f} +-{tol}')

    for th in (THALF_D300, THALF_POPPK):
        pk = TwoCompartmentPK(t_half_terminal=th).add_dose(0.0, 10.0)
        # 1. Terminal slope must reproduce the published half-life. Measured late, where the
        #    beta phase dominates, by regressing log-concentration on time.
        t = np.linspace(4 * th, 8 * th, 200)
        slope = np.polyfit(t, np.log(pk.conc(t)), 1)[0]
        check(f't_half terminal ({th} d target)', LN2 / (-slope), th, 0.05 * th)

        # 2. Extrapolated AUC of a single IV dose must equal dose/CL. This is the strongest
        #    available check that CL is entering the model correctly.
        tt = np.linspace(0, 200 * th, 400000)
        auc = np.trapezoid(pk.conc(tt), tt)
        check(f'AUC = dose/CL ({th} d)', auc, 10.0 / CL_DEFAULT, 0.02 * 10.0 / CL_DEFAULT)

        # 3. Initial concentration of an IV bolus is dose/V1 by definition.
        check(f'C(0) = dose/V1 ({th} d)', float(pk.conc(0.0)), 10.0 / V1_DEFAULT, 1e-6)

    # 4. Feasible-range helper must agree with where the solver actually fails.
    for th in (THALF_D300, THALF_POPPK):
        vmax = feasible_v_ratio_max(th)
        solve_two_compartment(th, v_ratio=vmax * 0.999)          # must succeed
        try:
            solve_two_compartment(th, v_ratio=vmax * 1.001)      # must fail
            print(f'  [FAIL] feasible_v_ratio_max({th}) did not bound the solver')
            ok = False
        except ValueError:
            print(f'  [PASS] feasible_v_ratio_max({th} d) = {vmax:.3f} bounds the solver')

    # 5. Linearity: doubling the dose doubles the concentration (published as dose-proportional).
    a = TwoCompartmentPK().add_dose(0.0, 10.0)
    b = TwoCompartmentPK().add_dose(0.0, 20.0)
    check('dose proportionality', float(b.conc(5.0) / a.conc(5.0)), 2.0, 1e-9)

    # 6. Superposition sanity: two 5 mg doses at the same instant equal one 10 mg dose.
    c = TwoCompartmentPK().add_dose(0.0, 5.0).add_dose(0.0, 5.0)
    check('superposition', float(c.conc(3.0)), float(a.conc(3.0)), 1e-12)

    print('\nPK self-test:', 'PASS' if ok else 'FAIL')
    return ok


if __name__ == '__main__':
    import sys
    if not _selftest():
        sys.exit(1)
