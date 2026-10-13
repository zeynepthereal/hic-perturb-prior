"""Gene annotation (GENCODE v44, GRCh38) and gene -> TSS-bin mapping for the Hi-C graphs.

Genome build consistency: Hi-C = ENCODE ENCFF616PUW (GRCh38/hg38), annotation = GENCODE v44 primary assembly (GRCh38.p14),
GEARS/scPerturb read-out coordinates (also GRCh38) are only a fallback for genes whose Ensembl ID is absent from v44.
Chromosomes: chr1-22 and chrX only (K562 is female -> no chrY signal; chrM has no Hi-C).
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

HIC_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
GTF = "data/annotation/gencode.v44.primary_assembly.basic.annotation.gtf.gz"


def read_gencode_genes(gtf_path: str | Path) -> pd.DataFrame:
    rows = []
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t")
            if c[2] != "gene":
                continue
            a = c[8]
            rows.append((re.search(r'gene_id "([^"]+)"', a).group(1).split(".")[0],
                         re.search(r'gene_name "([^"]+)"', a).group(1),
                         re.search(r'gene_type "([^"]+)"', a).group(1), c[0], int(c[3]), int(c[4]), c[6]))
    df = pd.DataFrame(rows, columns=["gene_id", "symbol", "gene_type", "chr", "start", "end", "strand"])
    df["tss"] = np.where(df.strand == "+", df.start, df.end)          # 1-based 5' end of the gene
    return df.drop_duplicates("gene_id")   # PAR genes appear on chrX and chrY with the same ID: keep first (chrX)


def tss_bin(tss: np.ndarray, res: int) -> np.ndarray:
    return (np.asarray(tss) - 1) // res


# GRCh38 chromosome lengths as reported by the Hi-C file header (ENCFF616PUW)
CHROM_LEN = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
             "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
             "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
             "chr19": 58617616, "chr20": 64444167, "chr21": 46709983, "chr22": 50818468, "chrX": 156040895}
