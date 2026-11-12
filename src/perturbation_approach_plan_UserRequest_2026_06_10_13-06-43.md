# Task Specification: Perturbative Construction of Diffusion Coefficients from a Transport Matrix

## Background

A previous Python script named Matrix Builder has already been implemented.

Matrix Builder:

1. Reads a transport graph stored on disk.
2. Constructs the corresponding transport matrix.
3. Exports the matrix to disk in Mathematica format.

The next step is to develop a new Python script that implements a perturbative hierarchy for calculating diffusion coefficients from this transport matrix.

---

## Input Requirements

The script must read from disk:

1. The transport matrix.
2. Row and column labels corresponding to the matrix indices.
3. Metadata identifying the node type associated with each row/column.

### Recommended File Format

In addition to Mathematica export, the matrix should also be exported in a Python-friendly format.

Preferred options:

- HDF5 (.h5) - recommended for large sparse matrices.
- NumPy compressed archive (.npz) - simple and efficient.
- JSON (only for small matrices and debugging).

The file must contain:

- matrix values,
- row labels,
- column labels,
- node type information.

---

## Assumptions About Node Ordering

The transport matrix is expected to be lexicographically ordered.

Node labels correspond to increasing solution energy:

A < B < C < D < E < ...

where:

- A represents the lowest solution-energy site,
- B the next-lowest,
- etc.

Before any calculations begin, the script must verify that the matrix ordering is consistent with this convention.

If the ordering is invalid, the script should terminate with a clear error message.

---

## Main Algorithm

### Objective

Construct a sequence of increasingly accurate approximations to the diffusion coefficient by progressively adding higher-energy node classes.

The procedure should not stop after the first percolating network is found.

Instead, it should continue until all node classes have been included.

---

## Iterative Procedure

### Step 1

Identify all nodes belonging to class A.

Construct the reduced transport matrix containing only A-type nodes.

### Step 2

Evaluate whether transport is possible in the selected direction.

The matrix being analyzed is already projected onto the desired transport direction, so no further directional processing is required.

As a first implementation, use the determinant as the criterion:

- det(M) != 0 -> transport/percolation exists.
- det(M) = 0 -> transport/percolation does not exist.

### Step 3

If transport exists:

- Compute the diffusion coefficient associated with the current reduced matrix.
- Store the result.

### Step 4

Add the next node class.

Example:

A -> A+B -> A+B+C -> A+B+C+D -> ...

For each enlargement:

1. Construct the corresponding reduced matrix.
2. Evaluate its determinant.
3. If transport exists:
   - Compute the diffusion coefficient.
   - Store the result.
4. Continue to the next level regardless of whether transport already existed at a previous level.

---

## Termination Condition

The algorithm stops only after all node classes have been included.

At that point the final diffusion coefficient corresponds to the complete transport network.

---

## Expected Output

For every perturbation level k:

- included node classes,
- matrix size,
- determinant,
- percolation status,
- calculated diffusion coefficient.

Example:

| Level | Included Classes | Matrix Size | Determinant | Percolates | Diffusion Coefficient |
|---------|---------|---------|---------|---------|---------|
| 1 | A | ... | ... | No | N/A |
| 2 | A+B | ... | ... | Yes | D2 |
| 3 | A+B+C | ... | ... | Yes | D3 |
| 4 | A+B+C+D | ... | ... | Yes | D4 |

The complete history should be saved to disk.

---

## Software Design Requirements

Implement the code in a modular way:

### Matrix Loader

Responsible for:

- reading matrices,
- reading labels,
- validating ordering.

### Matrix Reducer

Responsible for:

- selecting node classes,
- constructing reduced matrices.

### Percolation Analyzer

Responsible for:

- determinant evaluation,
- percolation detection.

### Diffusion Calculator

Responsible for:

- computing the diffusion coefficient for a given reduced matrix.

The actual diffusion-coefficient formula should be implemented as a separate function so that alternative transport models can be inserted later.

### Result Manager

Responsible for:

- storing intermediate results,
- exporting tables,
- generating convergence plots.

---

## Future Extensions

The determinant criterion is only the first implementation.

The code should be designed so that later it can be replaced by:

- rank-based criteria,
- null-space analysis,
- graph-theoretical percolation checks,
- Green-function methods,
- perturbative corrections based on block-matrix theory.

Therefore the percolation test should be encapsulated in a dedicated class or function.
