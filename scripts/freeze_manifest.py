"""python -m scripts.freeze_manifest -> results/hic_v2/MANIFEST.json : checksums of every artifact the pre-registered analysis depends on."""
import hashlib, json, subprocess
from pathlib import Path

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
G = Path("data/hic/graphs")
files = sorted(G.glob("b500*")) + sorted(G.glob("b250*")) + sorted(G.glob("hic_local25k*")) + sorted(Path("configs/splits").glob("*.json")) + \
        sorted(Path("configs/prereg").glob("*.yaml")) + sorted(Path("data/k562_regions_raw").glob("*")) + [Path("data/liftover/hg19ToHg38.over.chain.gz")] + \
        [Path(p) for p in ("results/hic/gene_table.csv", "results/hic/symbol_map.csv", "data/processed/k562_gears.npz", "configs/k562_rearranged_regions_hg38.bed",
                           "PREREGISTRATION.md", "results/hic_v2/ab_500k.csv", "results/hic_v2/ab_250k.csv", "results/hic_v2/zero_rows.csv")]
files = [f for f in files if "__rewire" not in f.name]       # rewired graphs are built lazily from fixed seeds (not frozen as files)
m = {"git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
     "git_dirty": bool(subprocess.check_output(["git", "status", "--porcelain", "--", ".", ":!results", ":!data"], text=True).strip()),
     "sources": {"hic": "ENCODE ENCFF616PUW (md5 d4b999471d4c17a74f1e66370c0012e3, remote read)", "annotation": "GENCODE v44 md5 fc4731226d2d315e6b2415fd1690403d",
                 "gears_h5ad_zip": "md5 84eed779531de88f83d3eb16a773a261", "k562_cnv_bed": "ENCFF486MJU md5 b4a78fa77789291cb5e3d528b32aa55a",
                 "k562_sv_vcf": "ENCFF863MPP md5 04f1c1f3aad184c0e05c863da7a7d3e0", "chain": "UCSC hg19ToHg38 md5 35887f73fe5e2231656504d1f6430900"},
     "sha256": {str(f): sha(f) for f in files if f.exists()}}
Path("results/hic_v2/MANIFEST.json").write_text(json.dumps(m, indent=1))
print(len(m["sha256"]), "files hashed; commit", m["git_commit"][:7], "dirty" if m["git_dirty"] else "clean")
assert not m["git_dirty"], "commit all code/doc changes first: the manifest must be generated on a clean tree"
