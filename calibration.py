"""
Versioned calibration artifact and the hard gate on downstream experiments.

WHY A GATE
----------
`exp_replicate.py` and `exp_tarlatamab.py` construct `Lymphoid` with shipped defaults. Running the
handoff commands would therefore have used UNCALIBRATED parameters while appearing to run the
calibrated experiment. Codex's review raised this as finding 3 and it has been open since.

The gate makes that failure impossible rather than merely documented: R and T load a calibration
artifact and refuse to execute unless its status is VALIDATED.

WHY "FITTED" IS NOT "VALIDATED"
-------------------------------
The status is not a description of fit quality. It records whether the model has earned predictive
authority, which requires a provenance-clean external prediction. As `provenance.py` asserts, every
Philipp observation is contaminated through `recover_tau` and `exhaust_tonic`, so no currently
available observation can raise the status above FIT_NOT_VALIDATED however good the fit becomes.

That is the intended behaviour. A model that fits its own calibration data is not licensed to
predict treatment schedules.
"""
import json, hashlib, os, subprocess
from dataclasses import dataclass, asdict, field

STATUS = ('VALIDATED', 'FIT_NOT_VALIDATED', 'UNDERIDENTIFIED', 'MODEL_CLASS_REJECTED', 'UNCALIBRATED')
ARTIFACT = 'results/calibration.json'


def _git_rev():
    try:
        return subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    except Exception:
        return 'unknown'


def _hash_files(paths):
    h = hashlib.sha256()
    for p in sorted(paths):
        if os.path.exists(p):
            h.update(p.encode())
            h.update(open(p, 'rb').read())
    return h.hexdigest()[:16]


@dataclass
class Calibration:
    status: str
    model: str
    params: dict
    residuals: dict
    fit_targets: list
    held_out: dict
    seeds: list
    notes: str
    boundary: list = field(default_factory=list)
    code_hash: str = ''
    git_rev: str = ''

    def save(self, path=ARTIFACT):
        assert self.status in STATUS, f'bad status {self.status}'
        self.code_hash = _hash_files(['lymphoid.py', 'philipp_assay.py', 'exhaustion.py',
                                      'fit_calibration.py', 'provenance.py'])
        self.git_rev = _git_rev()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            json.dump(asdict(self), f, indent=1)
        return self

    @staticmethod
    def load(path=ARTIFACT):
        if not os.path.exists(path):
            return Calibration(status='UNCALIBRATED', model='', params={}, residuals={},
                               fit_targets=[], held_out={}, seeds=[],
                               notes='no calibration artifact exists')
        return Calibration(**json.load(open(path)))


def require_validated(experiment_name, path=ARTIFACT):
    """Hard gate. Call at the top of any downstream experiment.

    Refuses on anything but VALIDATED, including a good fit. Deliberately not overridable by a
    keyword argument: an override would be used.
    """
    c = Calibration.load(path)
    if c.status != 'VALIDATED':
        raise RuntimeError(
            f"{experiment_name} refused: calibration status is {c.status}, not VALIDATED.\n"
            f"  model:  {c.model or '(none)'}\n"
            f"  notes:  {c.notes}\n"
            f"  A fitted model is not a validated one. Downstream experiments require a "
            f"provenance-clean external prediction, which provenance.py currently shows is "
            f"unavailable: every Philipp observation is contaminated via recover_tau and "
            f"exhaust_tonic.")
    return c


def _gate_tests():
    """Verify the gate actually blocks. A gate that has never refused anything is not a gate."""
    import tempfile
    ok = True

    def check(name, cond):
        nonlocal ok
        ok &= bool(cond)
        print(f'  [{"PASS" if cond else "FAIL"}] {name}')

    with tempfile.TemporaryDirectory() as d:
        base = dict(model='M0', params={}, residuals={}, fit_targets=[], held_out={},
                    seeds=[], notes='test')
        for st in ('FIT_NOT_VALIDATED', 'UNDERIDENTIFIED', 'MODEL_CLASS_REJECTED', 'UNCALIBRATED'):
            p = os.path.join(d, f'{st}.json')
            Calibration(status=st, **base).save(p)
            try:
                require_validated('experiment R', p)
                check(f'{st} refuses R', False)
            except RuntimeError:
                check(f'{st} refuses R', True)

        p = os.path.join(d, 'missing.json')
        try:
            require_validated('experiment T', p)
            check('absent artifact refuses T', False)
        except RuntimeError:
            check('absent artifact refuses T', True)

        p = os.path.join(d, 'ok.json')
        Calibration(status='VALIDATED', **base).save(p)
        try:
            require_validated('experiment R', p)
            check('VALIDATED permits R', True)
        except RuntimeError:
            check('VALIDATED permits R', False)

    print('\nGate tests:', 'PASS' if ok else 'FAIL')
    return ok


if __name__ == '__main__':
    import sys
    if not _gate_tests():
        sys.exit(1)
    print('\nCurrent artifact status:', Calibration.load().status)
