"""
Source-to-parameter provenance graph, and an automated contamination check.

WHY THIS EXISTS
---------------
A held-out observation was designated in `fit_calibration.py` and turned out to be contaminated:
`lymphoid.py` sets `recover_tau = 10080` with the comment "7 d; TFI reinvigoration (Philipp 2022,
Weber Science 2021)". The recovery constant was chosen by looking at the same paper's
treatment-free-interval result that the held-out point is drawn from. The fit left the parameter
alone, which is not the same as the parameter being independent of the observation.

Nothing in the pipeline checked for that, which is why it got through. This module makes the
failure mechanically detectable rather than dependent on someone remembering.

THE RULE
--------
A validation observation is clean only if its SOURCE informed none of:
  - any fitted parameter,
  - any fixed default relevant to the prediction,
  - the structure of the observation model, where that structure was chosen from this result,
  - the acceptance threshold, where that was chosen after seeing the outcome.

A parameter does not become independent of an observation because a later optimiser did not
modify it.
"""

# --------------------------------------------------------------------------- observations
OBSERVATIONS = {
    'philipp_d7_cont':  dict(source='philipp2022', assay='72h specific lysis, E:T 1:1, vs control construct',
                             value=88.4, figure='Fig 2E', n=6,
                             used_for_model_design=False, used_for_parameter_default=True,
                             used_for_fitting=True, used_for_validation=False),
    'philipp_d14_cont': dict(source='philipp2022', assay='72h specific lysis', value=34.9, figure='Fig 3E', n=6,
                             used_for_model_design=False, used_for_parameter_default=True,
                             used_for_fitting=True, used_for_validation=False),
    'philipp_d28_cont': dict(source='philipp2022', assay='72h specific lysis', value=8.6, figure='Fig 2E', n=6,
                             used_for_model_design=False, used_for_parameter_default=True,
                             used_for_fitting=True, used_for_validation=False),
    'philipp_d14_tfi':  dict(source='philipp2022', assay='72h specific lysis after 7d drug-free', value=93.4,
                             figure='Fig 3E', n=6,
                             used_for_model_design=False, used_for_parameter_default=True,
                             used_for_fitting=False, used_for_validation=True),
    # A fifth value exists and was not previously recorded. Same source, so equally contaminated:
    # it cannot rescue the validation problem.
    'philipp_d28_tfi':  dict(source='philipp2022', assay='72h specific lysis after TFI', value=58.7,
                             figure='Fig 3E', n=6,
                             used_for_model_design=False, used_for_parameter_default=True,
                             used_for_fitting=False, used_for_validation=False),
}

# --------------------------------------------------------------------------- parameters
PARAMETERS = {
    'p_kill':          dict(meaning='per-contact kill hazard per minute', status='fitted',
                            informed_by=['halle2016'], consumers=['chronic', 'readout', 'R', 'T']),
    'exhaust_tonic':   dict(meaning='exhaustion accrued per contact-minute at full occupancy',
                            status='fitted', informed_by=['philipp2022'],
                            consumers=['chronic', 'R', 'T']),
    'recover_tau':     dict(meaning='exhaustion recovery time constant', status='FIXED default',
                            informed_by=['philipp2022', 'weber2021'],
                            consumers=['chronic', 'R', 'T']),
    't_motility':      dict(meaning='T-cell move attempts per minute', status='fixed',
                            informed_by=['miller2002'], consumers=['chronic', 'readout', 'R', 'T']),
    'exhaust_per_kill': dict(meaning='exhaustion accrued per kill event', status='fixed at 0',
                             informed_by=['obertopp2025'], consumers=['chronic', 'R', 'T']),
}

SOURCES = {
    'philipp2022': 'Philipp et al., Blood 2022, PMID 35878001 (AMG 562 chronic stimulation + TFI)',
    'weber2021':   'Weber et al., Science 2021 (rest restores CAR-T function)',
    'halle2016':   'Halle et al., Immunity 2016, PMID 26872694 (CTL killing capacity in vivo)',
    'miller2002':  'Miller/Cahalan, Science 2002 (T-cell motility in lymph node)',
    'obertopp2025': 'Obertopp/Basanta, PMC12667981 (kill-driven PD-1 accrual)',
}


def contaminating_parameters(obs_id):
    """Parameters whose value was informed by the source of `obs_id`. Non-empty means contaminated."""
    src = OBSERVATIONS[obs_id]['source']
    return {p: d for p, d in PARAMETERS.items() if src in d['informed_by']}


def validation_is_clean(obs_id):
    return len(contaminating_parameters(obs_id)) == 0


def audit():
    ok = True
    print('Provenance audit of every observation marked for validation:\n')
    any_val = False
    for oid, o in OBSERVATIONS.items():
        if not o['used_for_validation']:
            continue
        any_val = True
        bad = contaminating_parameters(oid)
        clean = not bad
        ok &= clean
        print(f"  {oid}: source={o['source']}  -> {'CLEAN' if clean else 'CONTAMINATED'}")
        for p, d in bad.items():
            print(f"      contaminated via {p} ({d['status']}), informed by {d['informed_by']}")
    if not any_val:
        print('  (no observation is currently marked for validation)')

    # This assertion is the permanent regression test for the leak that occurred.
    assert not validation_is_clean('philipp_d14_tfi'), \
        'philipp_d14_tfi must never be reported as a clean validation target'
    print("\n  [PASS] philipp_d14_tfi is permanently flagged contaminated (recover_tau cites its source)")

    print('\nProvenance audit:', 'CLEAN' if ok else 'CONTAMINATED VALIDATION PRESENT')
    return ok


if __name__ == '__main__':
    audit()
