# Host-host diplon (fcc), test of the multi-orientation-state technique

This example tests whether SymDCT's site/jump model can represent a
"diplon" defect (two atoms sharing one substitutional site, in one of
several discrete orientations) as discussed in Ishioka & Koiwa (1985),
Section 4.3. A diplon isn't a simple site occupancy problem -- each site
carries an internal orientation degree of freedom, and diffusion involves
two physically distinct processes: **rotation** (in-place reorientation,
rate `w_R`) and **translation** (hop to a neighbouring site, rate `w_T`).

## Modelling technique

SymDCT has no built-in notion of "site with internal state." Instead,
each orientation is represented as a distinct pseudo-species occupying
the *same* coincident position:

- `Qx`, `Qy`, `Qz` -- three coincident dummy species per diplon site,
  one per cube-axis orientation. (Pymatgen's dummy-species mechanism is
  used here; note it rejects any label whose prefix matches a real
  element or isotope symbol -- e.g. `Dx` is rejected because `D` is
  deuterium -- which is why `Qx`/`Qy`/`Qz` were chosen instead.)
- **Rotation edges**: a small positive cutoff (0.05 Å) between different
  orientation-species at the same site, which the code's normal
  nearest-neighbour search finds as zero-distance neighbour pairs.
- **Translation edges**: the ordinary nearest-neighbour cutoff (2.0 Å,
  the simple-cubic cell edge) between *same*-orientation species at
  different sites.

## Simplifying assumption (differs from the literal fcc geometry)

The real fcc host-host diplon (Ishioka & Koiwa, Fig. 6-7) has
*directional* translation correlated with orientation -- an x-oriented
diplon does not hop isotropically. Reproducing that exactly would need
direction-specific cutoffs, which the geometry stage doesn't support
(cutoffs are purely radial); it would require pruning the
auto-discovered jump list by hand via `tweak_transition_rates`, the way
the paper's own two-frequency correlation problems are handled elsewhere
in this package.

This example instead uses a **decoupled, isotropic-translation** model
placed on a simple-cubic lattice (lattice parameter 2.0 Å): every
orientation hops to all 6 nearest neighbours at the same rate, and
reorients independently at each site. This is a different (simpler)
physical model than the literal fcc interstitialcy mechanism, so the
numeric prefactor in the result below does not match Ishioka & Koiwa's
published `D = (2/3) w_T a^2` -- but it tests exactly the mechanism in
question: coincident multi-orientation sites, rotation edges, and
translation edges, all discovered and reduced correctly by the existing
pipeline.

## Result

Running the full pipeline (`jump_network_builder` -> `matrix_builder` ->
`transport_coefficient_solver` -> `transport_coefficient_evaluator`)
along the [111] direction gives:

- All three orientations (`Qx`, `Qy`, `Qz`) are correctly grouped into a
  single equivalence class along [111], exactly as Ishioka & Koiwa state
  for the real fcc case.
- Interplanar spacing along [111] comes out as `a**2/3`, matching
  `(a/sqrt(3))**2` -- independently confirming the `a/sqrt(3)` period
  Ishioka & Koiwa give in their Fig. 7 for the true fcc lattice.
- Symbolic result: `D = a**2 * (W_{QxQx} + W_{QyQy} + W_{QzQz}) / 3`,
  i.e. `D = w_T * a**2` once the three (symmetry-equivalent) translation
  rates are set equal. **The rotation rates do not appear in the result
  at all** -- reproducing Ishioka & Koiwa's stated finding that the
  diffusion coefficient is independent of `w_R`.
- Numeric evaluation (`energetics.yaml`: translation barrier 0.5 eV,
  rotation barrier 0.2 eV) gives a perfectly linear Arrhenius plot
  (R² = 1.0000) with an apparent migration energy of exactly 0.500 eV --
  the translation barrier alone, with zero contribution from the 0.2 eV
  rotation barrier. This is an independent numerical confirmation of the
  symbolic result above.

## Files

- `input/host-host_diplon.poscar` -- simple-cubic structure, one diplon
  site per cell, three coincident `Qx`/`Qy`/`Qz` pseudo-species.
- `input/host-host_diplon.yaml` -- run parameters (axis [1,1,1]).
- `input/cutoffs.yaml` -- rotation and translation cutoffs.
- `input/energetics.yaml` -- example numeric energetics (filled in by
  hand for this test, not derived from any real material).
- `output/` -- generated jump network, symbolic/numeric results, and the
  Arrhenius plot.
