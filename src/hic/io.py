"""Remote (HTTP range-request) reading of the K562 Hi-C map with hicstraw; only derived sparse matrices are cached locally.

Source: ENCODE ENCFF616PUW.hic, experiment ENCSR545YBD (K562 in situ Hi-C, Rao et al. 2014 data), GRCh38, 4 replicates merged,
'mapping quality thresholded contact matrix', 22.64 GB, ENCODE-listed md5 d4b999471d4c17a74f1e66370c0012e3 (not verifiable without
downloading the file; the file is NOT downloaded). Raw counts are fetched ('NONE'); balancing (ICE) is done in graphs.py because this
file has no KR vectors at the resolutions used."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .genes import HIC_CHROMS

URL = "https://encode-public.s3.amazonaws.com/2022/10/10/0332f11a-aa2a-4532-9615-5b42b55ea5c3/ENCFF616PUW.hic"
ENCODE_ID, ENCODE_MD5 = "ENCFF616PUW", "d4b999471d4c17a74f1e66370c0012e3"


def open_hic():
    import hicstraw
    return hicstraw.HiCFile(URL)


def chrom_lengths(h) -> dict:
    return {c.name: c.length for c in h.getChromosomes()}


def _recs(rec) -> np.ndarray:
    return np.array([(r.binX, r.binY, r.counts) for r in rec], dtype=np.float64).reshape(-1, 3)


def fetch_cis(chrom: str, res: int, max_dist: int, out_dir: str | Path, win: int = 30_000_000):
    """Raw intra-chromosomal counts within `max_dist` bp, as (i, j, count) with i <= j in bin units. Fetched in diagonal windows
    so that only a band (not the full triangle) is ever held in memory."""
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out = Path(out_dir) / f"{chrom}.npz"
    if out.exists():
        return str(out), 0
    h = open_hic()
    L = chrom_lengths(h)[chrom]
    z = h.getMatrixZoomData(chrom, chrom, "observed", "NONE", "BP", res)
    parts = []
    for x0 in range(0, L, win):
        r = _recs(z.getRecords(x0, min(x0 + win, L), x0, min(x0 + win + max_dist, L)))
        if len(r):
            keep = (r[:, 0] >= x0) & (r[:, 0] < x0 + win) & (r[:, 1] - r[:, 0] <= max_dist)
            parts.append(r[keep])
    a = np.concatenate(parts) if parts else np.zeros((0, 3))
    np.savez_compressed(out, i=(a[:, 0] // res).astype(np.int32), j=(a[:, 1] // res).astype(np.int32), c=a[:, 2].astype(np.float32),
                        n_bins=-(-L // res), res=res, chrom=chrom)
    return str(out), len(a)


def fetch_pair(args):
    """Raw counts for one chromosome pair (cis or trans) at coarse resolution -> dense float32 block (bins_a x bins_b)."""
    a, b, res, out_dir = args
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out = Path(out_dir) / f"{a}__{b}.npy"
    if out.exists():
        return str(out), 0
    h = open_hic()
    L = chrom_lengths(h)
    z = h.getMatrixZoomData(a, b, "observed", "NONE", "BP", res)
    r = _recs(z.getRecords(0, L[a], 0, L[b]))
    M = np.zeros((-(-L[a] // res), -(-L[b] // res)), dtype=np.float32)
    if len(r):
        np.add.at(M, ((r[:, 0] // res).astype(int), (r[:, 1] // res).astype(int)), r[:, 2].astype(np.float32))
    np.save(out, M)
    return str(out), len(r)
