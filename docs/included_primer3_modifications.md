# Included Primer3 modifications

The C sources under `primer3/src/libprimer3/` are the included copy of the
open-source Primer3 v2.6.1 library, modified to make them amenable to the
Python C API / Cython bindings (see also [](miscellany.md#derivative-work)). In
addition to the pervasive build-glue and API-surface changes, we occasionally
patch this included C to fix bugs that are reachable through the bindings.

This page records those **intentional, behavior-affecting divergences** from
upstream Primer3 so they can be reviewed and re-applied if the included tree is
ever updated to a newer upstream release. It is a living document: add an entry
here whenever you change the included C in a way that alters behavior (pure
mechanical/build changes and comments do not need to be listed).

Each entry notes the file and function, the category of change, what changed
and why, and whether it is reachable from the Python API.

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
