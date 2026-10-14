"""STUDY 2 / Phase 1: RPE1 Hi-C QC gate ONLY (no Perturb-seq data of any kind is read here).
python -m scripts.study2_rpe1_hic_qc {header|fetch|qc}
  header : remote header only (genome, chromosomes, resolutions) + traffic
  fetch  : read ALL chromosome-pair blocks of ENCODE ENCFF876OWE remotely at 500 kb and 250 kb (raw counts; the .hic is never downloaded) -> data/hic_rpe1/pairs_<res>kb
  qc     : total contacts, ICE valid-bin fraction (same filter as Study 1), node-level all-zero-row fractions of the Hi-C graph and of one matched random-trans control at K=10/50
Reuses the Study 1 functions unchanged (src/hic/bins.py, graphs.valid_bins, spectral); the K562 reference numbers come from the cached Study 1 pair blocks."""
import json
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

URL = "https://encode-public.s3.amazonaws.com/2021/12/10/499fc588-5601-463d-9da3-fdb675e9224f/ENCFF876OWE.hic"   # ENCODE ENCFF876OWE, GRCh38, 2,240,397,369 B, md5 b6828cb642d1fce2b64c5383f22ec4f9
CH = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
OUT = Path("data/hic_rpe1")


def net_in() -> int:
    return int(subprocess.check_output("netstat -ib | awk '$1==\"en0\" && $3 ~ /Link/ {print $7; exit}'", shell=True, text=True))


def fetch_pair(args):
    a, b, res = args
    f = OUT / f"pairs_{res // 1000}kb" / f"{a}__{b}.npy"
    if f.exists():
        return 0
    import hicstraw
    h = hicstraw.HiCFile(URL); L = {c.name: c.length for c in h.getChromosomes()}
    z = h.getMatrixZoomData(a, b, "observed", "NONE", "BP", res)
    rec = z.getRecords(0, L[a], 0, L[b])
    M = np.zeros((-(-L[a] // res), -(-L[b] // res)), dtype=np.float32)
    if len(rec):
        r = np.array([(x.binX, x.binY, x.counts) for x in rec], dtype=np.float64)
        np.add.at(M, ((r[:, 0] // res).astype(int), (r[:, 1] // res).astype(int)), r[:, 2].astype(np.float32))
    f.parent.mkdir(parents=True, exist_ok=True); np.save(f, M)
    return len(rec)


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "header":
        import hicstraw
        n0 = net_in(); h = hicstraw.HiCFile(URL)
        print("genome:", h.getGenomeID()); print("chromosomes:", [(c.name, c.length) for c in h.getChromosomes()][:6], "...", len(h.getChromosomes()))
        print("resolutions:", h.getResolutions()); print(f"traffic for the header read: {(net_in()-n0)/1e6:.1f} MB (all en0 traffic)")
    elif what == "fetch":
        n0, t0 = net_in(), time.time()
        for res in (500_000, 250_000):
            jobs = [(a, b, res) for i, a in enumerate(CH) for b in CH[i:]]
            with Pool(3) as p:
                tot = sum(p.imap_unordered(fetch_pair, jobs))
            print(f"res {res}: {len(jobs)} blocks, {tot} records, elapsed {time.time()-t0:.0f}s, traffic so far {(net_in()-n0)/1e6:.0f} MB", flush=True)
    elif what == "qc":
        from src.hic.bins import BIN_CFG, Genome, build_cells, ice, load_raw, oe_matrix
        from src.hic.spectral import spectral_embedding
        from src.hic.bins import rewire_trans
        rep = {}
        for name, root in (("RPE1", OUT), ("K562 (Study 1 reference)", Path("data/hic"))):
            for res in (500_000, 250_000):
                d = root / f"pairs_{res // 1000}kb"
                G = Genome(res)
                tot = 0.0
                for ia, a in enumerate(CH):
                    for b in CH[ia:]:
                        M = np.load(d / f"{a}__{b}.npy"); tot += float(np.triu(M).sum() if a == b else M.sum())
                R = load_raw(G, d); cov = R.sum(1, dtype=np.float64)
                x, valid = ice(R, BIN_CFG["ice_iters"], BIN_CFG["ice_tol"])
                rep[f"{name} {res//1000}kb"] = {"total_contacts_in_pair_blocks": tot, "bins": G.n, "valid_bins": int(valid.sum()), "valid_fraction": float(valid.mean()),
                                                "bins_zero_coverage": int((cov == 0).sum()), "median_bin_coverage": float(np.median(cov[valid])), "cv_bin_coverage": float(cov[valid].std() / cov[valid].mean())}
                if name == "RPE1" and res == 500_000:                                                         # graph-level node zero-row gate (no perturbation targets needed)
                    OE, mt, _ = oe_matrix(R, x, valid, G, BIN_CFG["min_pairs"]); cells = build_cells(OE, valid, G); T = cells.pop("_trans_only")
                    Tr = rewire_trans(T, G.chr_idx, 0); ctrl = (cells["lin_cis"] + Tr).tocsr()
                    zr = {}
                    for gname, W in (("hic_trans", cells["hic_trans"]), ("lin_randtrans_s0", ctrl)):
                        r_ = spectral_embedding(W, k_max=50)
                        E, ok = r_["emb"], r_["embedded"].astype(bool)
                        zr[gname] = {"components": r_["n_components"], "isolated_nodes": r_["n_isolated"], **{f"node_rows_allzero_K{K}": float((np.abs(E[ok][:, :K]) < 1e-8).all(1).mean()) for K in (10, 50)}}
                    rep["RPE1 500kb graph gate (node level)"] = zr
        Path("results/study2").mkdir(parents=True, exist_ok=True)
        json.dump(rep, open("results/study2/rpe1_hic_qc.json", "w"), indent=1)
        print(json.dumps(rep, indent=1))
