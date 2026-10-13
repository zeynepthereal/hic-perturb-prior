"""Build results/hic/gene_table.csv (graph nodes) + results/hic/symbol_map.csv + results/hic/mapping_report.txt.

Nodes = all GENCODE v44 protein-coding genes on chr1-22,X  UNION  every GEARS dataset gene (read-outs and perturbation targets)
that can be placed on chr1-22,X. Targets/read-outs are resolved by Ensembl ID first, then GENCODE symbol, then (read-outs only)
the coordinates stored in the h5ad."""
import json
from collections import Counter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from src.hic.genes import GTF, HIC_CHROMS, read_gencode_genes

OUT = Path("results/hic"); OUT.mkdir(parents=True, exist_ok=True)
gen = read_gencode_genes(GTF)
by_id, by_sym = gen.set_index("gene_id"), gen.drop_duplicates("symbol").set_index("symbol")
print("GENCODE v44 genes:", len(gen), "| protein-coding on chr1-22,X:",
      int(((gen.gene_type == "protein_coding") & gen.chr.isin(HIC_CHROMS)).sum()))

# --- GEARS dataset genes
g = ad.read_h5ad("data/gears/replogle_k562_essential/perturb_processed.h5ad", backed="r")
ro = g.var.reset_index()[["gene_name", "gene_id", "chr", "start", "end", "strand"]].copy()
ro["gene_name"] = ro.gene_name.astype(str)
cond = [c for c in g.obs["condition"].astype(str).unique() if c != "ctrl"]
targets = sorted({c.replace("+ctrl", "") for c in cond})
# symbol -> Ensembl ID for targets, from the scPerturb obs (Replogle's own target IDs)
sp = ad.read_h5ad("data/raw/ReplogleWeissman2022_K562_essential.h5ad", backed="r").obs[["perturbation", "gene_id"]]
sp = sp.astype(str).drop_duplicates().set_index("perturbation").gene_id.to_dict()
ro_sym2id = dict(zip(ro.gene_name, ro.gene_id.astype(str)))


def resolve(sym, gid, fallback_coords=None):
    """returns (gene_id, symbol, chr, tss, strand, method, fail_reason)"""
    if gid and gid in by_id.index:
        r = by_id.loc[gid]; return gid, r.symbol, r.chr, int(r.tss), r.strand, "id", None
    if sym in by_sym.index:
        r = by_sym.loc[sym]; return r.gene_id, sym, r.chr, int(r.tss), r.strand, "symbol", None
    if fallback_coords is not None:
        c = fallback_coords; tss = int(c.start) + 1 if c.strand == "+" else int(c.end)
        return gid, sym, c.chr, tss, c.strand, "h5ad_coords", None
    return gid, sym, None, None, None, None, "id_and_symbol_not_in_gencode_v44" if gid else "symbol_not_found"


recs = []
for _, r in ro.iterrows():
    out = resolve(r.gene_name, str(r.gene_id), r)
    recs.append(("readout", r.gene_name, *out))
for t in targets:
    gid = ro_sym2id.get(t) or sp.get(t)
    out = resolve(t, gid)
    recs.append(("target", t, *out))
m = pd.DataFrame(recs, columns=["role", "query", "gene_id", "symbol", "chr", "tss", "strand", "method", "fail_reason"])
m.loc[m.chr.notna() & ~m.chr.isin(HIC_CHROMS) & m.fail_reason.isna(), "fail_reason"] = \
    m.chr.map(lambda c: "chrY_no_hic_signal(K562 female)" if c == "chrY" else ("chrM_no_hic" if c == "chrM" else "non_standard_contig"))
ok = m.fail_reason.isna()

lines = []
for role in ("readout", "target"):
    s = m[m.role == role]
    so = s[s.fail_reason.isna()]
    lines.append(f"{role}: {len(s)} queried | mapped {len(so)} "
                 f"(by id {int((so.method=='id').sum())}, by symbol {int((so.method=='symbol').sum())}, "
                 f"by h5ad coords {int((so.method=='h5ad_coords').sum())}) | failed {int(s.fail_reason.notna().sum())}")
    for k, v in Counter(s.fail_reason.dropna()).most_common():
        lines.append(f"    fail: {k}: {v}  e.g. {', '.join(s[s.fail_reason==k]['query'].head(8))}")
tgt_in_ro = m[(m.role == "target")]["query"].isin(set(ro.gene_name))
lines.append(f"targets also present as read-outs: {int(tgt_in_ro.sum())}; targets NOT read-outs: {int((~tgt_in_ro).sum())} "
             f"(mapped {int((m.role.eq('target') & ~m['query'].isin(set(ro.gene_name)) & ok).sum())})")
lines.append("id_changes: symbol of Ensembl ID differs between dataset and GENCODE v44: " +
             str(int((ok & (m.method == 'id') & (m['query'] != m.symbol)).sum())))

# --- node table
ds_nodes = m[ok & m.gene_id.notna()].drop_duplicates("gene_id")[["gene_id", "symbol", "chr", "tss", "strand"]]
pc = gen[(gen.gene_type == "protein_coding") & gen.chr.isin(HIC_CHROMS)][["gene_id", "symbol", "chr", "tss", "strand"]]
nodes = pd.concat([ds_nodes, pc]).drop_duplicates("gene_id").reset_index(drop=True)
nodes["is_readout"] = nodes.gene_id.isin(set(m[(m.role == "readout") & ok].gene_id))
nodes["is_target"] = nodes.gene_id.isin(set(m[(m.role == "target") & ok].gene_id))
nodes["bin25k"] = ((nodes.tss - 1) // 25_000).astype(int)
nodes["bin500k"] = ((nodes.tss - 1) // 500_000).astype(int)
nodes = nodes.sort_values(["chr", "tss"], key=lambda s: s.map({c: i for i, c in enumerate(HIC_CHROMS)}) if s.name == "chr" else s).reset_index(drop=True)
nodes.to_csv(OUT / "gene_table.csv", index=False)
# symbol -> gene_id (dataset names + current GENCODE names), for embeddings
sm = {}
for _, r in m[ok].iterrows():
    sm[r["query"]] = r.gene_id
for _, r in nodes.iterrows():
    sm.setdefault(r.symbol, r.gene_id)
# scPerturb targets too (secondary pipeline)
for k, v in sp.items():
    if v in set(nodes.gene_id):
        sm.setdefault(k, v)
pd.Series(sm, name="gene_id").rename_axis("symbol").reset_index().to_csv(OUT / "symbol_map.csv", index=False)
m.to_csv(OUT / "gears_gene_mapping.csv", index=False)
lines.append(f"nodes: {len(nodes)} (protein-coding {len(pc)}, extra dataset genes {len(nodes)-len(pc)}), "
             f"readouts on graph {int(nodes.is_readout.sum())}, targets on graph {int(nodes.is_target.sum())}")
lines.append(f"genes sharing a 25 kb TSS bin with another node: {int(nodes.duplicated(['chr','bin25k'], keep=False).sum())}; "
             f"sharing a 500 kb bin: {int(nodes.duplicated(['chr','bin500k'], keep=False).sum())}")
(OUT / "mapping_report.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
