"""Re-run only the L1 cells lost to the 2026-09-07 out-of-memory stall at 125/150.
Writes each result to its own file immediately, so a further stall loses at most one run.
Identical code path and seeds as exp_schedule.run - this is a resumption, not a new experiment."""
import json, os, sys, glob
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exp_schedule as S

OUT = S.OUT
DONE = json.load(open(f'{OUT}/expL1_salvage125.json'))
have = {(r['arch'], r['tfi'], r['seed']) for r in DONE}
need = [(a, k, s) for a in S.ARCHS for k in S.TFIS for s in S.SEEDS if (a, k, s) not in have]
os.makedirs(f'{OUT}/L1_rerun', exist_ok=True)


def one(job):
    a, k, s = job
    f = f'{OUT}/L1_rerun/{a}_{k}_{s}.json'
    if os.path.exists(f):
        return f
    r = S.run((a, k, s))
    json.dump(r, open(f, 'w'))
    return f


if __name__ == '__main__':
    print(f'missing {len(need)} cells: {sorted(set((a,k) for a,k,_ in need))}', flush=True)
    with Pool(3) as p:                      # 3 not 6: the stall was an out-of-memory event
        for i, f in enumerate(p.imap_unordered(one, need), 1):
            print(f'  {i}/{len(need)} {os.path.basename(f)}', flush=True)
    merged = list(DONE)
    for f in glob.glob(f'{OUT}/L1_rerun/*.json'):
        merged.append(json.load(open(f)))
    json.dump(merged, open(f'{OUT}/expL1.json', 'w'))
    open(f'{OUT}/expL1_log.txt', 'w').write(
        f'expL1: {len(merged)} runs (125 from first launch, {len(need)} re-run after an '
        f'out-of-memory stall at 125/150 on 2026-09-07)\n')
    print(f'merged {len(merged)} runs -> expL1.json')
