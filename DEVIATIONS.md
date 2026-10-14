# DEVIATIONS — log of everything that deviates from, or post-dates, the frozen pre-registration

`PREREGISTRATION.md` is **frozen at git tag `prereg-v1`** and is **never edited again**. Everything that changes after the freeze — tooling, corrections, departures from the registered plan —
is recorded here, append-only, with date, reason and the commits involved. An entry never retroactively alters what was registered.

**OSF registration link:** (to be added)

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
