<p align="center">
  <img src="../images/SymDCT_logo.png" alt="SymDCT logo" width="100%"/>
</p>

<p align="center">
  A symbolic toolkit for analytical diffusion coefficient tensor calculations
  in crystals from atomic-scale kinetic data.
</p>

---
# Current Status
> **Status:** Pre-release / active development. First public version (v0.1.0)
> targeted for **September 1, 2026**. See [ROADMAP.md](symdct/ROADMAP.md) for details.

# What is SymDCT?
**SymDCT** is

- **Graph-based** — represents diffusion networks as their underlying connectivity, independent of a particular lattice realization
- **Symmetry-based** — exploits crystallographic point-group symmetry to reduce the transport problem to its irreducible components
- **Exact** — yields closed-form (analytical) diffusion coefficients, not numerical approximations
- **General** — applicable to complex, multi-site crystal structures beyond simple Bravais lattices
- **Interoperable** — compatible with both first-principles (DFT) and ML interatomic potential-based energetics/kinetics as input

# Main Features

## Theoretical foundation
- **Exact, closed-form solutions** — diffusion coefficients derived analytically from the transport matrix, not fitted or extracted from simulation trajectories
- **Symmetry-adapted formulation** — crystallographic point-group symmetry reduces the full transport problem to its irreducible sublattice/site-type components
- **Graph-based network representation** — diffusion pathways are encoded as connectivity graphs, decoupling the topology of hops from the specific geometry of any one lattice

## Computational machinery
Present version uses **sympy** for symbolic calculations.

## Generality and scope
- **Applicable beyond simple Bravais lattices** — handles complex, multi-site crystal structures (e.g., tetragonal Be₁₂Ti with distinct interstitial site types)
- **Topological transferability** — results carry over across topologically equivalent diffusion networks, even when the underlying lattice differs
- **Separation of geometry and energetics** — the same symbolic/graph machinery applies regardless of the energy source

## Interoperability
- **First-principles or ML-potential input** — migration barriers, binding energies, and attempt frequencies can come from DFT or machine-learned interatomic potentials interchangeably
- **Symbolic + numerical backends** — implemented with SymPy for symbolic derivation, with Mathematica as an alternative backend for cross-validation

# Scientific Background
Follows and significantly extends the framework developed by Ishioka and Koiwa in their 1985 paper [1].


# Roadmap

# Installation

# Examples

# Documentation

# Citation

# References

[1] S. Ishioka and M. Koiwa, "Diffusion coefficient in crystals with
  multiple jump frequencies," *Philosophical Magazine A* **52**,
  267–277 (1985). [DOI: 10.1080/01418618508237623](https://doi.org/10.1080/01418618508237623)


# Contact
Email: pavel.vladimirov@kit.edu; pavel.v.vladimirov1962@gmail.com
