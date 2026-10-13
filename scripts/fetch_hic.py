"""python -m scripts.fetch_hic cis|trans   (remote reads; caches sparse/dense derived matrices in data/hic/)"""
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
    else:
        Path("data/hic/pairs_500kb").mkdir(parents=True, exist_ok=True)
        jobs = [(a, b, TRANS_RES, "data/hic/pairs_500kb") for i, a in enumerate(HIC_CHROMS) for b in HIC_CHROMS[i:]]
        with Pool(4) as p:
            for path, n in p.imap_unordered(fetch_pair, jobs):
                print(f"{path} records={n} t={time.time()-t0:.0f}s", flush=True)
