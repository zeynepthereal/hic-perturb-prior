"""python -m scripts.make_k562_regions -> configs/k562_rearranged_regions_hg38.bed  (secondary analysis S10 exclusion list)
Inputs (hg19, ENCODE; downloaded to data/k562_regions_raw/): ENCFF486MJU.bed.gz, ENCFF863MPP.vcf.gz; chain: data/liftover/hg19ToHg38.over.chain.gz (UCSC)."""
import gzip, hashlib, re
from pyliftover import LiftOver

CHAIN = "data/liftover/hg19ToHg38.over.chain.gz"
lo = LiftOver(CHAIN)
def lift(c, p):
    r = lo.convert_coordinate(c, p)
    return (r[0][0], r[0][1]) if r else None
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
amp, fail = [], 0
for l in gzip.open("data/k562_regions_raw/ENCFF486MJU.bed.gz", "rt").read().splitlines():
    c = l.split("\t")
    if c[3] != "amp":
        continue
    a, b = lift(c[0], int(c[1])), lift(c[0], int(c[2]) - 1)
    if a and b and a[0] == b[0] and b[1] >= a[1]:
        amp.append((a[0], a[1], b[1] + 1, "amp", f"logR={c[9]}"))
    else:
        fail += 1
W, bp, fail2, seen = 1_000_000, [], 0, set()
for l in gzip.open("data/k562_regions_raw/ENCFF863MPP.vcf.gz", "rt"):
    if l[0] == "#":
        continue
    c = l.split("\t"); m = re.search(r"(chr\w+):", c[4])
    if not m or m.group(1) == c[0]:
        continue                                   # intra-chromosomal breakends are not translocations
    p = lift(c[0], int(c[1]))
    if not p:
        fail2 += 1; continue
    if (p[0], p[1] // 1000) in seen:
        continue
    seen.add((p[0], p[1] // 1000))
    bp.append((p[0], max(0, p[1] - W), p[1] + W, "tra_bp", f"hg19 {c[0]}:{c[1]} -> {m.group(1)}"))
reg = sorted(amp + bp, key=lambda x: (x[0], x[1]))
hdr = ['# K562 rearranged/amplified regions, GRCh38 -- secondary-analysis exclusion list (PREREGISTRATION.md S10)',
       f'# (i)  all segments called "amp" in ENCODE ENCFF486MJU (ENCSR000BBO, HAIB Illumina genotyping array, hg19, md5 {md5("data/k562_regions_raw/ENCFF486MJU.bed.gz")})',
       f'# (ii) +-1 Mb around every inter-chromosomal breakpoint in ENCODE ENCFF863MPP (ENCSR053AXS, 10X linked-read SV calls grocsvs, hg19, md5 {md5("data/k562_regions_raw/ENCFF863MPP.vcf.gz")}), all PASS',
       f'# hg19->hg38: UCSC hg19ToHg38.over.chain.gz (md5 {md5(CHAIN)}) via pyliftover; segment ends lifted separately, kept if same chromosome and ordered',
       '# columns: chrom start end type info']
open("configs/k562_rearranged_regions_hg38.bed", "w").write("\n".join(hdr) + "\n" + "\n".join("\t".join(map(str, r)) for r in reg) + "\n")
print(f"amp segments lifted {len(amp)} (failed {fail}), total {sum(e-s for _,s,e,_,_ in amp)/1e6:.1f} Mb | breakpoint windows {len(bp)} (failed {fail2}) | regions {len(reg)}")
