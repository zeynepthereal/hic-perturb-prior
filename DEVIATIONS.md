# DEVIATIONS — log of everything that deviates from, or post-dates, the frozen pre-registration

`PREREGISTRATION.md` is **frozen at git tag `prereg-v1`** and is **never edited again**. Everything that changes after the freeze — tooling, corrections, departures from the registered plan —
is recorded here, append-only, with date, reason and the commits involved. An entry never retroactively alters what was registered.

**OSF registration link:** (to be added; registration postponed, see D4)

---

## D1 — Tag `prereg-v1` moved once, for a header-only correction (2026-10-14)

- **What happened.** The tag was created on the freeze commit and, minutes later (commit timestamps: 05:10:52 → 05:14:18 local time), moved to a new commit.
  The first version of the frozen file still carried a draft header ("DRAFT v0.5 … Not frozen") and a future-tense line 4; the owner asked for exactly those two lines to be corrected.
- **Before any use of the tag.** No analysis, no Hi-C perturbation model and no OSF upload had been run or made; nothing referenced the tag.
- **Old tag:** object `c2e530d5fa1dc5ceef902c61b44d0b8ad9df4e70` → commit `43fa01c32d1f00beb61204c42d01f37d2f193ded` (deleted locally and on `origin`).
- **New tag:** object `bdda47a99d9201226eba99d5f63e9cbd99e77be5` → commit `7bb6c399992a06a236b8698c28e51eea2934a8be` (annotated, pushed; `git ls-remote --tags origin` shows the peeled commit equal to the then-HEAD and to `origin/main`).
- **Intermediate commit** (header correction only): `a7db2dd965928319fd7a0ebb3cc30ba03cd6da9e`.
- **Exact difference between the two tag targets** (`git diff 43fa01c 7bb6c39`): 2 files.
  - `PREREGISTRATION.md`, lines 3–4 only (status line → "Status: FROZEN v1.0 (2026-10-14), tag prereg-v1. No Hi-C perturbation model had been run at freeze."; line 4 → "The manifest (`results/hic_v2/MANIFEST.json`) was regenerated on a clean tree at freeze, and the commit was tagged `prereg-v1`.").
  - `results/hic_v2/MANIFEST.json`: the recorded code commit (`534b22c…` → `a7db2dd…`) and the SHA-256 of `PREREGISTRATION.md` (`2f3588311eaee432…` → `d389a4655dbfd07e…`).
- **Analytical content:** identical to v0.5 (the body of the document, all configs, graphs, splits and code are byte-identical between the two tag targets).

## D2 — Tooling added after the freeze: `scripts/run_primary_analysis.py` (2026-10-14)

- A thin command-line wrapper around the frozen `src/stats.py::primary_analysis` (and, optionally, the `nulls.json` written by `src/nulls.py` via `h1_support`). Arguments: run directory, arm names
  (`--hic`, `--control` / `--control-template` + `--n-control`), model kind (`--kind paper|fixed`, `--K`, `--lam`), and the frozen defaults for eval set, target handling, gene subset, metric, permutations (100,000) and seed (0).
- **No analytical change.** It computes nothing itself: every number is returned by the frozen functions. The optional `--nulls` / `--tier2-from` flags only transcribe the Tier-1 / Tier-2 rule of §3 from results of the frozen code;
  **if the wrapper and the prose of §3 ever disagree, the prose governs.**
- **Test:** `src/test_cli.py` runs the one-liner that was reported at the freeze as a subprocess (default 100,000 permutations) and the wrapper on the same synthetic metrics file and asserts identical output JSON; a second test covers fixed-K / residualized model names and the tier logic.
- The code and configs registered at `prereg-v1` are untouched; this file, the wrapper and its test (`src/test_cli.py`) are the only additions of this commit.

## D3 — Open items

- OSF registration link: (to be added when the owner has uploaded the frozen `PREREGISTRATION.md`).

## D4 — OSF registration postponed; pre-analysis timestamp by OpenTimestamps (2026-10-14)

OSF registration is **postponed**. The pre-analysis timestamp is instead provided by **OpenTimestamps** (`opentimestamps-client` 0.7.2, `ots stamp`), made at **2026-10-14T20:39:09Z**, before any pre-registered analysis was run.

| file (byte-identical to the file at tag `prereg-v1` = commit `7bb6c399992a06a236b8698c28e51eea2934a8be`) | SHA-256 | proof |
|---|---|---|
| `PREREGISTRATION.md` | `d389a4655dbfd07e866e9ee7300b001bb8faacf945f5ac30114aa80bcd09e0ee` | `PREREGISTRATION.md.ots` |
| `results/hic_v2/MANIFEST.json` | `82421ab7098f6a5e8a1f55c6fc2e522ae9703ac19fd3ff6bd0feab87df3bba09` | `results/hic_v2/MANIFEST.json.ots` |

- Only the SHA-256 digests were submitted, to the public calendars `a.pool.opentimestamps.org`, `b.pool.opentimestamps.org`, `a.pool.eternitywall.com` and `ots.btc.catallaxy.com`.
- **Status: pending.** Bitcoin confirmation takes a few hours; run `ots upgrade <file>.ots` later (then commit the upgraded proofs) and verify with `ots verify PREREGISTRATION.md.ots`.
- **OSF link: to be added later.**

## D5 — Notes on the execution of the pre-registered run (2026-10-14)

None of these changes an analytical choice; they are recorded for completeness.

- **Sanity-check configs.** §7 checks 1 and 5 were run through the registered code path on subsets of the registered arms with two additional, non-registered configs, `configs/sanity/sanity1_pca_only.yaml` (PCA + mean baselines) and `configs/sanity/sanity5_lean_vs_standard.yaml` (two primary arms). They were committed before running (clean tree logged). Check 5 was then confirmed on the full primary run (22 arm×K combinations). No Hi-C vs control contrast was evaluated from the sanity runs.
- **Order of operations.** Checks 1, 3, 4 and 5 (on the two-arm sanity run) were evaluated and passed before the primary run was started; check 2 (reported, not a gate) is reported after the null run because it needs the null arms.
- **Storage of large outputs.** `per_perturbation_metrics.csv` of the primary (47.9 MB) and S3 (100.9 MB) runs are committed as `.csv.gz` (round-trip SHA-256 verified); the raw CSVs are git-ignored (`.gitignore`) to respect the 50 MB file limit. `scripts/run_primary_analysis.py` reads the CSV: `gunzip -k` the file first.
- **Rewired graphs.** The 20 degree-preserving-rewire graphs (and their embeddings) of the S2 null are built lazily from fixed seeds (0–19) at analysis time and are therefore not in the frozen manifest; the SHA-256 of all `*__rewire*` artifacts is recorded in `results/<null run>/rewire_artifacts_sha256.json`.
- **Results record.** `results/PREREG_ANALYSIS_SUMMARY.md` collects the numbers and the pre-registered labels; per-run details are in each run directory (config, git.json, per-case metrics).

## D6 — Pre-registered S4 / S9 executed with a new script (2026-10-14)

- **What:** the pre-registered secondary S4 (each arm vs the mean baseline on the primary test cases) and the S4/S9 PCA reference were run with `scripts/run_s4_vs_mean.py` (committed before the run, code commit `049388b`; run directory `results/*_prereg_s4_vs_mean`). The script only calls frozen functions (`paired_table`, `signflip_test`, `NullEvaluator`); D is defined as in the primary with the mean baseline as reference (positive = arm better).
- **Settings not fixed by the registration:** the 1,000 random-Gaussian arms use 10,000 permutations and 1,000 bootstrap draws per arm (the named arms use the frozen 100,000 / 10,000); the PCA reference uses the `common` set of the PCA-only sanity run (PCA-predictable subset, 517 cases), which is bit-identical to the earlier `k562_gears` run (sanity check 1). K = 50 uses the fixed K = 50, λ = 0.1 models of the primary run.

## D7 — EXPLORATORY, POST HOC, NOT PRE-REGISTERED: heterogeneity of the primary effect (2026-10-14)

- Done **after** the pre-registered results were known, at the owner's request. Script `scripts/EXPLORATORY_posthoc_heterogeneity.py`; output directory `results/*_EXPLORATORY_posthoc_primary_heterogeneity` (README states the same). Contents: distribution of the per-perturbation contributions a_i (quantiles, histogram), concentration of D in the largest contributors, per-split concentration, and the ten largest contributors to D in split 1 (target gene, chromosome, membership of the target bin in an S10 rearranged region).
- **Descriptive only: no tests, no p-values.** It changes no registered label or conclusion and is not part of the confirmatory analysis.

## D8 — OpenTimestamps upgrade attempt (2026-10-14T21:07Z)

`ots upgrade` on both proofs (about 28 minutes after stamping) did **not** complete: the calendars `alice.btc.calendar.opentimestamps.org` (tx `c781f9850e206433f0e42ba36c5c73fb67b1d868756b3225a713b227dcf91380`) and `bob.btc.calendar.opentimestamps.org` (tx `88c3c62fc5510473b91926b9fc5badff136f8cd505416167b12769fafaf6e2f6`) report the digests are included in Bitcoin transactions and are waiting for 6 confirmations; `btc.calendar.catallaxy.com` is pending; `finney.calendar.eternitywall.com` timed out. The `.ots` files are unchanged (identical SHA-256 before and after), so nothing was committed. To be retried.
