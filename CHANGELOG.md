# CHANGELOG

This file contains the canonical "Instructions for Future Agents" entry referenced
by `public/.clinerules`. It is intentionally minimal and designed to be read
before any non-trivial refactor or architectural change.

---

## 2026-04-15  - Initial placeholder entry

**Context:** The `public/.clinerules` file expects a `CHANGELOG.md` entry named
"Instructions for Future Agents" to be consulted before making architectural
changes. This is a repository-level placeholder to satisfy that requirement and
to surface a small set of high-priority instructions for automated agents.

**Instructions for Future Agents**

- BEFORE starting any refactor, architectural change, or addition that affects
  the public pipeline, read this file (the latest entry) and then `public/.clinerules`.
- If the current prompt contains an explicit override token `!!!`, the explicit
  directive in the prompt takes highest priority.
- Priority for resolving conflicting instructions (highest → lowest):
  1. Explicit directives in current prompt marked with `!!!`
  2. Latest entry in this `CHANGELOG.md` (section "Instructions for Future Agents")
  3. `public/.clinerules`
  4. Files under `public/.memory-bank/`
  5. Subpackage `.clinerules` files (for local constraints)

- DO NOT introduce runtime imports from any package located outside
  `public/`. If additional functionality is required, implement it
  inside `public/src/<subpkg>/` and document it in `public/.memory-bank/`.
- Follow the project rule: Python modules must live in `public/src/<subpkg>/` —
  avoid placing loose `.py` files directly under `public/src/`.

**How to add a new changelog entry**

When completing a non-trivial task, append a new entry here following the
template below. Ensure the section title includes a timestamp so agents can
select the latest entry easily.

### [YYYY-MM-DD HH:MM] - Brief description

**Context:** Why the change was needed
**Solution:** What was implemented
**Result:** Success / Partial / Rolled Back

**Key Files:**
- `path/to/file` — brief note

**Instructions for Future Agents:**
- Short bullet-list of agent-facing instructions (if any)

---

### [2026-10-28 15:00] - Lumpability check, merge policy, axis-orbit equivalence

**Context:** Merging equivalent sites along a transport axis was exact only if
the partition is lumpable (`R A(k) = A'(k) R` along `k = hkl * t`). The old
criterion (same Wyckoff letter / orbit and same plane) could merge sites related
only by operations that reverse the transport direction (e.g. inversion), and
merged orientation pseudo-species through the bare Wyckoff letter.
**Solution:**
- `src/matrix_builder/lumpability.py`: symbolic lumpability check,
  coarsest lumpable refinement, suggested rate relations.
- `matrix_builder.py`: YAML keys `merge_classes: check|force` and
  `site_equivalence: plane_orbit|orbit|wyckoff`.
- `src/symmetry_context/axis_equivalence.py`: new default `axis_orbit`
  (operations with `hkl W = hkl`, union-find over sites of equal species;
  screw axes and glides along xi may carry a site into another plane, as in
  Ishioka and Koiwa 1985, Sec. 3(b)). `plane_orbit` additionally requires
  integer `hkl . w`. Classes are numbered by the smallest `s = (hkl . frac) mod 1`.
- `examples/hcp_int_OT` and `test_hcp_interstitial_ot.py`: O/T interstitial in
  hcp, reproduces eqs. (29) and (31) of the paper; 3x3 along c.
- `conftest.py`: spglib error-handling flag set via `SPGLIB_OLD_ERROR_HANDLING`.
- Tests: `src/matrix_builder/tests/test_lumpability.py`,
  `src/symmetry_context/tests/test_axis_equivalence.py` (incl. P-1 counterexample).
- `scripts/check_lumpability_examples.py` prints the check for an example YAML.
**Result:** Success (tests pass; all examples lumpable, D unchanged between `plane_orbit` and `axis_orbit`).

**Key Files:**
- `src/matrix_builder/lumpability.py`, `src/matrix_builder/matrix_builder.py`
- `src/symmetry_context/axis_equivalence.py`
- `src/matrix_builder/README.md` — section on `site_equivalence` and `merge_classes`

**Instructions for Future Agents:**
- Do not merge sites by operations that reverse the transport direction (`hkl W = hkl` is required).
- Screw axes and glides along the transport direction are allowed; the matrix is in the position gauge (phases `exp(i k . h)` from real displacements).
- Keep the lumpability check as a safeguard when changing the merge rules.

---

### [2026-11-10 09:30] - Class labels use the smallest projection of a class

**Context:** Reduced-matrix labels are the site label plus a number that orders
classes of the same label along the transport axis. The number was taken from
the projection of the first site of a class. A class spanning several planes
(screw axis, glide plane) is placed at its smallest projection in
`axis_equivalence.py`, so the two orderings could differ.
**Solution:** `prepare_label_data_for_reordered_matrix` uses the smallest
projection over all sites of a class.
**Result:** Success (new test in `test_matrix_reduction_helpers.py`; label
numbering of the existing examples is unchanged).


**Key Files:**
- `src/matrix_builder/matrix_reduction_helpers.py`
- `src/matrix_builder/tests/test_matrix_reduction_helpers.py`

**Instructions for Future Agents:**
- Keep the class position used for labels consistent with `find_equiv_classes_for_axis` (smallest projection of the class).


---

### [2026-11-10 10:10] - Report lumpable merges, fix docstring of `to_structure_hkl`

**Context:** The lumpability check printed nothing when a merge of equivalent
sites passed, so a passed check was indistinguishable from a skipped one.
The docstring of `JumpNetworkSymmetryContext.to_structure_hkl` described the
conversion in the wrong direction (primitive to conventional); the code
converts conventional Miller indices to the structure basis.
**Solution:** `_apply_merge_policy` prints one line per axis when classes
were merged and the merge is lumpable. The docstring is corrected.
**Result:** Success (two new tests in `test_lumpability.py`; 78 tests pass with
the hcp vacancy tests below).

**Additional test:** `test_hcp_vacancy.py` runs the hcp vacancy example with one
species for both sites. Along (2,-1,0) and (0,0,1) the two sites are merged
(1x1, `axis_orbit`). Along (1,0,0), the <10-10> direction, no operation that
keeps the direction relates them, so they stay separate (2x2). D is the same
along (2,-1,0) and (1,0,0). `plane_orbit` merges only along (2,-1,0).

**Key Files:**
- `src/matrix_builder/matrix_builder.py`
- `src/matrix_builder/tests/test_lumpability.py`
- `src/matrix_builder/tests/test_hcp_vacancy.py`
- `src/symmetry_context/jump_symmetry_context.py`

---

_If this placeholder is insufficient for a particular task, add a clear
entry to this file explaining the intended behaviour before proceeding._

