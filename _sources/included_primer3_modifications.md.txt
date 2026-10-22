# Included Primer3 modifications

The C sources under `primer3/src/libprimer3/` are the included copy of the
open-source Primer3 v2.6.1 library, modified to make them amenable to the
Python C API / Cython bindings (see also [](miscellany.md#derivative-work)). In
addition to the pervasive build-glue and API-surface changes, we occasionally
patch this included C to fix bugs that are reachable through the bindings.

**The included source is Primer3 v2.6.1 plus the bug fixes listed on this
page.** Primer3 has not tagged a release since v2.6.1, but its `main` branch
has accumulated a number of fixes; the ones that affect results or memory
safety through the bindings have been backported here (see
[Upstream fixes backported onto v2.6.1](#upstream-fixes-backported-onto-v261)),
so calculations can differ from a vanilla v2.6.1 build in the cases described
below.

This page records those **intentional, behavior-affecting divergences** from
upstream Primer3 so they can be reviewed and re-applied if the included tree is
ever updated to a newer upstream release. It is a living document: add an entry
here whenever you change the included C in a way that alters behavior (pure
mechanical/build changes and comments do not need to be listed).

Each entry notes the file and function, the category of change, what changed
and why, and whether it is reachable from the Python API.

```{note}
The extension modules are built from the `*flex` variants of the library
(`libprimer3flex.c` / `thalflex.c`, derivatives of `libprimer3.c` / `thal.c`),
plus `oligotm.c`, `read_boulder.c` and `thal_parameters.c` directly. Where a
fix touches `libprimer3` or `thal` code it is applied to both the compiled
`*flex` file and its non-flex counterpart so the two do not drift apart.
```

---

## Upstream fixes backported onto v2.6.1

Each of these corresponds to a commit on the Primer3 `main` branch
([primer3-org/primer3](https://github.com/primer3-org/primer3)) after the
v2.6.1 tag. Where an upstream hash is given, it is the commit the fix was taken
from.

### `oligotm.c` — length term of the Owczarzy 2008 Mg²⁺ correction

- **Category:** correctness fix (integer division)
- **Upstream:** `c5b092b`; primer3-py issue
  [#178](https://github.com/libnano/primer3-py/issues/178)
- **Reachable from Python:** yes — `calc_tm` / `bindings.calc_tm` with
  `salt_corrections_method='owczarzy'` and `dv_conc > 0`, in the regimes where
  Eq. 16 of Owczarzy et al. (2008) applies (no monovalent cations, or
  √[Mg²⁺]/[Mon⁺] ≥ 0.22).

The length-dependent term `1/(2(N−1))·[e + f·ln[Mg²⁺] + g·ln²[Mg²⁺]]` was
written `1/(2 * (len - 1))` with `len` an `int`, so it was integer division and
always `0`. The term is now computed in floating point. Tm changes by roughly
3.5–5.5 °C for 10–12-mers in Mg²⁺ buffers and by about 0.3–0.6 °C for typical
20–25-mer primers in 50 mM K⁺ / 1.5–3 mM Mg²⁺; against the measured duplexes in
Owczarzy et al. 2008 Table S1 the mean absolute error drops from 2.59 °C to
1.55 °C. Other salt-correction methods are unaffected.

### `thal_parameters.c` / `primer3_config/stackmm.dh` — single-mismatch stack enthalpies

- **Category:** correctness fix (thermodynamic parameters)
- **Upstream:** `d2844b5` (only the three numeric changes were taken; upstream
  also reformatted the parameter files and added citations)
- **Reachable from Python:** yes — `calc_hairpin`, `calc_homodimer`,
  `calc_heterodimer` (and the self-/pair-complementarity checks inside
  `design_primers`). `calc_tm` is not affected.

Three entries of the mismatch-stack enthalpy table (used for internal 1×1
loops, i.e. a single mismatch flanked by matched base pairs) were wrong:

| Stack (5'→3' / 3'→5') | Was | Now |
| --- | --- | --- |
| `AC_TC` | −10 cal/mol | 0 |
| `CT_CA` | −10 cal/mol | 0 |
| `TA_AA` | `inf` (forbidden) | 4700 cal/mol |

The defaults are the `stackmm_dh` string literal in `thal_parameters.c`; the
reference copy `primer3_config/stackmm.dh` was updated to match. The `AC_TC` /
`CT_CA` change moves ΔH by 10 cal/mol per such stack (negligible). The `TA_AA`
change makes a previously impossible mismatch context available to the
dynamic programming, so a small fraction of sequences (a T·A pair stacked on an
A·A mismatch) can now fold into a structure that was previously invisible;
for those, hairpin/dimer ΔG, ΔH and Tm can change substantially (e.g. a
hairpin that reported no structure under v2.6.1 now reports Tm ≈ 27 °C).
Results agree with upstream's `ntthal` built from `d2844b5`.

### `thal.c` / `thalflex.c` — `readTLoop` one-byte overwrite

- **Category:** memory-safety fix (out-of-bounds write)
- **Upstream:** `6461c0e`
- **Reachable from Python:** yes, on every thermodynamic calculation (the
  default triloop/tetraloop tables are parsed when the parameters are loaded).

`readTLoop` NUL-terminated its output with `s[5] = '\0'` / `s[6] = '\0'`, but
`struct triloop::loop` / `struct tetraloop::loop` are exactly 5 / 6 bytes, so
each call wrote one byte past the array. The terminator was never used (the
bytes are immediately converted in place to numeric codes), so it was removed;
results are unchanged.

### `libprimer3flex.c` / `libprimer3.c` — negative interval start after `first_index` correction

- **Category:** correctness / robustness fix
- **Upstream:** `a70770e`
- **Reachable from Python:** yes — `design_primers` with a non-zero
  `PRIMER_FIRST_BASE_INDEX` and an interval tag (`SEQUENCE_TARGET`,
  `SEQUENCE_EXCLUDED_REGION`, `SEQUENCE_INTERNAL_EXCLUDED_REGION`,
  `SEQUENCE_PRIMER_PAIR_OK_REGION_LIST`).

In `_check_and_adjust_1_interval`, `intervals[i][0] -= first_index` could make a
start position negative when the user's start was below `first_index`. The
start is now clamped to `0`. For example, with `PRIMER_FIRST_BASE_INDEX=1` an
excluded region starting at `0` previously raised a spurious `EXCLUDED_REGION
outside of INCLUDED_REGION` warning; it is now treated as starting at the first
base.

### `libprimer3flex.c` / `libprimer3.c` — ok-region optimization could discard valid product sizes

- **Category:** correctness fix
- **Upstream:** `e1232c2` (with one deliberate difference, below)
- **Reachable from Python:** yes — `design_primers` with
  `SEQUENCE_PRIMER_PAIR_OK_REGION_LIST`.

`_optimize_ok_regions_list` narrows each left/right ok-region pair using the
requested product-size range and primer length. The narrowing could produce
ranges that no longer admit the requested products. It now only narrows when the
result remains consistent, and a region pair that cannot yield any product in
the requested range is collapsed to a 1 bp region and reported as a
`PRIMER_WARNING`:
`PAIR_OK_REGIONS[l_start,l_len,r_start,r_len] does not allow products of min-maxbp`.
The function gained a `warning` argument to report this.

**Difference from upstream:** upstream's final feasibility check is
`pmin < re - ls + 1`, which wrongly rejects a region pair whose largest possible
product is *exactly* the smallest requested size (for example a fixed
`200-200` product range with regions spanning exactly 200 bp). primer3-py uses
`pmin <= re - ls + 1`.

### Reviewed but not ported

The upstream `main` branch also contains the following changes that were
evaluated and deliberately **not** taken. They are listed here so a future
re-vendoring can revisit them:

- A large, behavior-preserving refactor of `thal.c` (split into several
  translation units, removal of macros and global state, 2-D dynamic
  programming tables, default parameters compiled in as header tables).
- `oligotm.c` rewritten without macros, and the SantaLucia 2004 parameter
  table (`c5b092b`, `23fee1d`); only the Owczarzy fix was taken.
- Hardening of `malloc` failure paths (`f8b6706`), which mostly touches code
  that is not exposed through the bindings (`masker`, `amplicontm`).
- Downgrading "left and right primer are identical" from an error to a warning
  for user-supplied primers (`8cc57fe`, `027eb29`), a policy change rather than
  a bug fix.
- New settings such as `PRIMER_INTERNAL_REV_COMP_G` /
  `PRIMER_INTERNAL_NO_5_PRIME_G` (`010c4e8`).
- Fix for a `thal` dynamic-programming-table leak (`8964d5e`); not applicable
  because the included code allocates those tables as flat arrays freed with
  `free()`.

---

## Safety / correctness fixes

### `read_boulder.c` — `parse_intron_list` bounds check

- **Category:** memory-safety fix (out-of-bounds write)
- **Reachable from Python:** yes — `run_design` /
  `bindings.design_primers` via the `SEQUENCE_OVERLAP_JUNCTION_LIST` and
  `SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST` tags.

`parse_intron_list` wrote `list[*count]` and incremented `*count` with no
check against `PR_MAX_INTERVAL_ARRAY` (200). The destination arrays
(`sa->primer_overlap_junctions`, `sa->intl_overlap_junctions`) are fixed
`int[PR_MAX_INTERVAL_ARRAY]` fields of `seq_args`, so a junction list with more
than 200 integers overran the array and corrupted adjacent struct fields. The
sibling interval parsers already bounds-check; this one was the outlier.

The parser now resets `*count` and returns `0` (its existing error signal, on
which the caller appends an `"Error in SEQUENCE_..._JUNCTION_LIST"` message)
when the array is full.

### `oligotm.c` — reject sequences shorter than 2 nt in `oligotm()`

- **Category:** memory-safety + correctness fix (out-of-bounds read / nonsense
  result)
- **Reachable from Python:** yes — `calc_tm` / `bindings.calc_tm` via `seqtm`.

`oligotm()` computed `len = strlen(s) - 1` and then applied the
nearest-neighbor terminal penalty at `s[len]`. For an empty string `len`
underflowed to `-1`, so the code read `s[-1]` (out of bounds); a 1-nt sequence
produced a nonsense Tm instead of the `OLIGOTM_ERROR` sentinel that the
function's own `ERROR:` label documents for "length < 2".

`oligotm()` now returns `OLIGOTM_ERROR` up front when `strlen(s) < 2`. Empty
input still returns the same sentinel (now without the OOB read), and 1-nt
input returns the sentinel instead of a meaningless value.

---

## Notes

- Line numbers are intentionally omitted; search for the function name and the
  `primer3-py:` comment that marks each patch in the source.
- When updating the included tree from a newer upstream Primer3, check each
  entry above against the incoming source and re-apply any patch that upstream
  has not itself adopted.
