"""
Candidate screen: which cancer should a compute-only effort attack?

CRITERIA FIXED BEFORE SCORING. The point is to KILL candidates fast, the same discipline that
has worked all session at the hypothesis level, applied at the disease level.

We have: public data, arithmetic, literature access, and a validated instinct for checking our own
work. We do NOT have: a wet lab, clinical data access, patients, or money. The screen must select
for what that hand can actually play.

SEVEN CRITERIA
--------------
1  CURE_EXISTS   Is anyone cured of this today? If yes, the mechanism is proven to exist and the
                 question becomes "why not everyone", which is far more tractable than "is cure
                 possible at all".
2  DEPENDENCY    Is there a defined molecular dependency? Curable cancers overwhelmingly have one
                 (PML-RARA, BCR-ABL, platinum-sensitivity). No dependency, no lever.
3  GAP_SIZE      How many patients fail current therapy? Small gap = little to win. Huge gap with
                 no mechanism = not our problem to solve.
4  PUBLIC_DATA   Is there downloadable quantitative data? DepMap, TCGA, SEER, clinicaltrials.gov,
                 published PK. Without this we cannot do anything at all.
5  ACTIONABLE    Does a drug or trial exist that a finding could redirect? A finding nobody can act
                 on is a paper, not a contribution.
6  UNCROWDED     Can we add something 500 better-resourced labs have not? Inverse of field size.
7  COMPUTE_ONLY  Can the specific gap be attacked by analysis rather than experiment?

Scored 0-3 each. Any candidate scoring 0 on CURE_EXISTS, PUBLIC_DATA or ACTIONABLE is KILLED
outright regardless of total -- those are necessary, not merely desirable.
"""

C = ['cure_exists', 'dependency', 'gap_size', 'public_data', 'actionable', 'uncrowded', 'compute_only']
NECESSARY = ['cure_exists', 'public_data', 'actionable']

CANDIDATES = [
    # name, scores, note
    ('Testicular germ cell', [3,3,0,2,1,2,1],
     'Cured >95%. Gap is tiny. Nothing to win.'),
    ('Hodgkin lymphoma',     [3,2,1,2,2,1,1],
     'Cured >90%. Remaining gap is late toxicity, not cure.'),
    ('Papillary thyroid',    [3,2,0,2,1,2,1],
     '~98% survival. Overtreatment is the issue, not cure.'),
    ('APL',                  [3,3,2,3,2,2,1],
     'Cured >90% in trials but 17-40% early death in population registries. '
     'DELIVERY gap, already well documented, ECOG-ACRIN already showed a fix.'),
    ('Childhood ALL',        [3,3,2,3,3,0,2],
     '92% cured. The 8% who fail is a real gap, but the most heavily worked '
     'paediatric cancer on earth.'),
    ('DLBCL',                [3,2,3,3,3,0,2],
     '~60% cured, 40% fail. Real gap but CAR-T and bispecifics are swarming it.'),
    ('CML',                  [2,3,2,3,3,1,3],
     'Not cured, controlled. TFR question is a dynamics problem with superb serial '
     'data - but modelling groups already own it.'),
    ('Neuroblastoma HR',     [2,3,3,3,3,2,2],
     '~50% survival high-risk. MYCN/ecDNA dependency. Anti-GD2 proves immune approach '
     'works. Genuinely underworked relative to its lethality.'),
    ('Ewing sarcoma',        [2,3,3,3,2,2,3],
     'EWSR1-FLI1 fusion in ~85% - a SINGLE defining lesion, the property curable '
     'cancers share. Ewing lines are the most PARP-sensitive in DepMap, yet clinical '
     'PARP trials failed. That unexplained gap is a data question.'),
    ('Osteosarcoma',         [2,1,3,3,2,3,2],
     'Flat survival 40 years. Chromothripsis >70%, no recurrent driver. CIN itself may '
     'be the dependency (KIF18A). Rare, so commercially orphaned.'),
    ('Rhabdomyosarcoma',     [2,3,3,2,2,2,2],
     'PAX3-FOXO1 fusion driver in alveolar subtype. Similar shape to Ewing, less data.'),
    ('AML',                  [2,2,3,3,3,0,2],
     'Heterogeneous. Some subsets curable. Extremely crowded.'),
    ('Multiple myeloma',     [1,2,3,3,3,0,2],
     'Not curable but deeply treatable. Very crowded.'),
    ('Glioblastoma',         [0,1,3,3,2,1,1],
     'KILLED: nobody is cured. No proven mechanism to extend.'),
    ('Pancreatic',           [0,2,3,3,2,1,1],
     'KILLED: ~12% 5-year. Cure essentially absent outside early resection.'),
    ('SCLC',                 [0,2,3,3,3,1,2],
     'KILLED on cure_exists. Tarlatamab tail is promising but not cure.'),
    ('NSCLC metastatic',     [1,2,3,3,3,0,2],
     'Immunotherapy tail looks like functional cure in ~15-20%. Enormous field.'),
]


def score(s):
    return sum(s)


def killed(s):
    return [C[i] for i, v in enumerate(s) if v == 0 and C[i] in NECESSARY]


def report():
    rows = []
    for name, s, note in CANDIDATES:
        k = killed(s)
        rows.append((name, s, score(s), k, note))
    rows.sort(key=lambda r: (-len(r[3]) == 0, r[2]), reverse=True)

    print(f"{'candidate':>22s} {'tot':>4s}  " + ' '.join(f'{c[:4]:>4s}' for c in C) + '   status')
    print('-' * 96)
    for name, s, tot, k, note in rows:
        st = 'KILLED: ' + ','.join(k) if k else ''
        print(f"{name:>22s} {tot:4d}  " + ' '.join(f'{v:4d}' for v in s) + f'   {st}')
    print('\nSurvivors ranked:\n')
    for name, s, tot, k, note in rows:
        if k:
            continue
        print(f"  {tot:2d}  {name}")
        print(f"      {note}\n")


if __name__ == '__main__':
    report()
