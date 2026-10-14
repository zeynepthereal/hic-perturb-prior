"""python -m scripts.freeze_manifest -> results/hic_v2/MANIFEST.json : checksums of every artifact the pre-registered analysis depends on."""
import hashlib, json, subprocess
from pathlib import Path

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
files = sorted(Path("data/hic/graphs").glob("b500_*")) + sorted(Path("data/hic/graphs").glob("b250_*")) + sorted(Path("configs/splits").glob("*.json")) + \
        [Path(p) for p in ("results/hic/gene_table.csv", "results/hic/symbol_map.csv", "data/processed/k562_gears.npz", "configs/k562_gears_hic_primary.yaml")]
m = {"git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
     "git_dirty": bool(subprocess.check_output(["git", "status", "--porcelain", "--", ".", ":!results", ":!data"], text=True).strip()),
     "sources": {"hic": "ENCODE ENCFF616PUW (md5 d4b999471d4c17a74f1e66370c0012e3, remote read)", "annotation": "GENCODE v44 md5 fc4731226d2d315e6b2415fd1690403d",
                 "gears_h5ad_zip": "md5 84eed779531de88f83d3eb16a773a261"},
     "sha256": {str(f): sha(f) for f in files if f.exists()}}
Path("results/hic_v2/MANIFEST.json").write_text(json.dumps(m, indent=1))
print(len(m["sha256"]), "files hashed; commit", m["git_commit"][:7], "dirty" if m["git_dirty"] else "clean")
