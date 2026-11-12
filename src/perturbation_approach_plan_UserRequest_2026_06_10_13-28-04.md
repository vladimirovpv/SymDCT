# Task Specification: Perturbative Construction of Diffusion Coefficients from a Transport Matrix

## Background

A previous Python module named `matrix_builder` has already been implemented.

`matrix_builder`:

1. Reads a transport graph stored on disk.
2. Constructs the corresponding transport matrix.
3. Exports the matrix to disk in Mathematica format.

The next step is to develop a new Python script that implements a perturbative hierarchy for calculating diffusion coefficients from this transport matrix.

This requirements file refines the initial plan saved in:

`perturbation_approach_plan_UserRequest_2026_06_10_13-06-43.md`

---

## Input Requirements

The script must read from disk:

1. A square transport matrix.
2. Row labels corresponding to the matrix indices.
3. Metadata associated with the matrix export.

The matrix is square, and row and column labels are the same. Therefore the file only needs to store one ordered label list, unless a future export format requires row and column labels separately for compatibility.

---

## Recommended File Format

The preferred file format for the next implementation is JSON.

The file must contain:

- matrix values,
- row labels,
- metadata.

JSON is preferred at this stage because the matrices are currently not large and because the format is easy to inspect while the algorithm is still being refined.

---

## Required Metadata

The JSON export produced by `matrix_builder` should be slightly extended so that the perturbative diffusion script can reconstruct the physical context of the matrix.

The metadata should include:

- Miller indices of the selected transport orientation,
- interplanar distance in the selected direction,
- `min_arg`, if `min_arg` was substituted during matrix construction,
- optional site class information if this is not already recoverable from the row labels.

The earlier phrase "node type information" should be interpreted carefully. A generic node type may not be necessary. What is required for the perturbative hierarchy is the ability to group sites into ordered energy/site classes such as `A`, `B`, `C`, etc. If this class can be reliably parsed from the row labels, no separate node type field is required. If it cannot be parsed reliably, then the JSON metadata or per-site records should include an explicit site class for each matrix label.

---

## Assumptions About Node Ordering

The transport matrix is expected to be ordered according to increasing solution energy.

Site classes follow:

A < B < C < D < E < ...

where:

- A represents the lowest solution-energy site class,
- B the next-lowest site class,
- etc.

Before any calculations begin, the script must verify that the matrix ordering is consistent with this convention.

If the ordering is invalid, the script should terminate with a clear error message.

---

## Main Algorithm

### Objective

Construct a sequence of increasingly accurate approximations to the diffusion coefficient by progressively adding higher-energy site classes.

The procedure should not stop after the first connected or percolating partial network is found.

Instead, it should continue until all site classes have been included.

---

## Iterative Procedure

### Step 1

Identify all sites belonging to class A.

Restrict the full transport matrix to the rows and columns corresponding only to A-type sites.

The full transport matrix read from disk should not be modified in place. The restricted matrix can be copied for each perturbation level to avoid accidental changes to the original matrix.

### Step 2

Evaluate whether transport is possible in the selected direction.

The matrix being analyzed has already been projected or constructed for the desired transport direction, so no further directional processing is required.

For the first implementation, it is acceptable to calculate the determinant of the restricted transport matrix itself, because the matrices are currently not large.

As an initial criterion:

- det(M) != 0 -> transport/percolation exists.
- det(M) = 0 -> transport/percolation does not exist.

### Step 3

If transport exists:

- Compute the diffusion coefficient associated with the current restricted matrix.
- Store the result.

The diffusion-coefficient calculation should be delegated to the implementation in `transport_coefficient_solver/`. The specific solver function is not yet implemented and should be added there rather than embedded directly in the perturbation script.

### Step 4

Add the next site class.

Example:

A -> A+B -> A+B+C -> A+B+C+D -> ...

For each enlargement:

1. Restrict the full transport matrix to the rows and columns for the currently included site classes.
2. Evaluate connectivity/percolation using the encapsulated percolation criterion.
3. If transport exists:
   - Compute the diffusion coefficient.
   - Store the result.
4. Continue to the next level regardless of whether transport already existed at a previous level.

---

## Connectivity and Percolation Criterion

For a partial graph, where only a subset of the full graph sites is included, connectivity in the selected transport direction can be tested in several ways.

Possible criteria include:

- determinant of the restricted transport matrix,
- determinant of the graph Laplacian,
- a numerically stable rank or null-space criterion,
- a graph-theory path search checking whether paths span the crystal cell from one edge to another in the selected direction.

For the time being, use the determinant of the restricted transport matrix itself.

This check must be encapsulated in a dedicated function or class named `PercolationCriterion`, so that the determinant-based implementation can later be replaced without changing the perturbative loop.

---

## Termination Condition

The algorithm stops only after all site classes have been included.

At that point the final diffusion coefficient corresponds to the complete transport network in the selected orientation.

---

## Expected Output

For every perturbation level k:

- included site classes,
- matrix size,
- determinant or criterion value,
- percolation/connectivity status,
- calculated diffusion coefficient,
- metadata copied from the source matrix export, including orientation and interplanar distance.

Example:

| Level | Included Classes | Matrix Size | Criterion Value | Percolates | Diffusion Coefficient |
|---------|---------|---------|---------|---------|---------|
| 1 | A | ... | ... | No | N/A |
| 2 | A+B | ... | ... | Yes | D2 |
| 3 | A+B+C | ... | ... | Yes | D3 |
| 4 | A+B+C+D | ... | ... | Yes | D4 |

The complete history should be saved to disk.

---

## Software Design Requirements

Implement the code in a modular way.

### Matrix Loader

Responsible for:

- reading the JSON matrix export,
- reading matrix labels,
- reading metadata,
- validating that the matrix is square,
- validating that row and column labels are identical if both are present,
- validating ordering of site classes.

### Matrix Restrictor

Responsible for:

- selecting included site classes for each perturbation level,
- constructing restricted copies of the full matrix for those sites.

This replaces the earlier "Matrix Reducer" concept. The matrix has already been reduced or projected to a particular transport direction before it is read from disk. The perturbative script should only restrict the already-loaded full matrix to the necessary subset of sites.

### PercolationCriterion

Responsible for:

- determinant evaluation for the first implementation,
- percolation/connectivity detection,
- providing an interface that can later support graph-theory, Laplacian, rank, or null-space criteria.

### Diffusion Calculator

Responsible for:

- computing the diffusion coefficient for a given restricted matrix.

The actual diffusion-coefficient formula should be implemented in `transport_coefficient_solver/` and called from the perturbative script.

### Result Manager

Responsible for:

- storing intermediate results,
- exporting tables,
- generating convergence plots.

---

## Future Extensions

The determinant criterion is only the first implementation.

The code should be designed so that later it can be replaced by:

- graph Laplacian determinant checks,
- rank-based criteria,
- null-space analysis,
- graph-theoretical spanning-path checks,
- Green-function methods,
- perturbative corrections based on block-matrix theory.

Therefore the percolation test should remain encapsulated in `PercolationCriterion`.
