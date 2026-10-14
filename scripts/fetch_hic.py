"""python -m scripts.fetch_hic cis|trans [res]   (remote reads; caches sparse/dense derived matrices in data/hic/)"""
import sys
import time
from multiprocessing import Pool
from pathlib import Path

from src.hic.genes import HIC_CHROMS
from src.hic.io import chrom_lengths, fetch_cis, fetch_pair, open_hic

CIS_RES, MAX_DIST, TRANS_RES = 25_000, 10_000_000, 500_000

def _cis(c): return fetch_cis(c, CIS_RES, MAX_DIST, "data/hic/cis_25kb")

if __name__ == "__main__":
    what = sys.argv[1]
    t0 = time.time()
    if what == "cis":
        Path("data/hic/cis_25kb").mkdir(parents=True, exist_ok=True)
        L = chrom_lengths(open_hic())
        order = sorted(HIC_CHROMS, key=lambda c: -L[c])
        with Pool(4) as p:
            for path, n in p.imap_unordered(_cis, order):
                print(f"{path} records={n} t={time.time()-t0:.0f}s", flush=True)
    else:   # `trans` (500 kb, default) or `trans 250000`: ALL chromosome pairs incl. the full cis range, raw counts
        res = int(sys.argv[2]) if len(sys.argv) > 2 else TRANS_RES
        out = f"data/hic/pairs_{res // 1000}kb"
        Path(out).mkdir(parents=True, exist_ok=True)
        jobs = [(a, b, res, out) for i, a in enumerate(HIC_CHROMS) for b in HIC_CHROMS[i:]]
        with Pool(4) as p:
            for path, n in p.imap_unordered(fetch_pair, jobs):
                print(f"{path} records={n} t={time.time()-t0:.0f}s", flush=True)
