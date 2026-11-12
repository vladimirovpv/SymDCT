# Third-party licenses

This file lists the third-party software that SymDCT needs, with version and license.
SymDCT itself is released under the MIT license (see `LICENSE`).

SymDCT does not copy or redistribute any third-party Python package. Users install the
packages themselves from PyPI, for example with `uv sync`. The licenses below therefore
apply to the packages as installed. They do not change the license of the SymDCT source code.

Date of the audit: 2026-11-10.

A machine-readable software bill of materials (SBOM) is provided in `sbom.cdx.json`
(CycloneDX 1.5). It lists the same packages with name and version, as recorded in `uv.lock`,
and contains no license fields. The licenses are given in the tables below.
It is regenerated after every change of the dependencies with

```bash
uv export --format cyclonedx1.5 --output-file sbom.cdx.json
```

## Summary

- Total: 10 direct runtime packages, 30 transitive runtime packages, 4 test packages
  (pytest and its three dependencies).
- No package under GPL, LGPL-only or AGPL terms is required.
- No package conflicts with distributing SymDCT under MIT.
- Four items need a remark (flags F1 to F4 below). None blocks the release.
- Two further packages appear in `uv.lock` only for Windows (colorama, tzdata). They were not
  installed in the audit environment. Their licenses are listed in the section "Platform-specific packages".

## Method

1. Direct dependencies: `[project].dependencies` and `[dependency-groups]` in `pyproject.toml`.
2. Real imports: all `import` and `from ... import` lines in `src/`, `main.py`, `scripts/`,
   `tests/`, `examples/` and `conftest.py` were collected. Third-party top-level modules found:
   matplotlib, networkx, numpy, pint, pymatgen, ruamel.yaml, scipy, spglib, sympy, tqdm, pytest.
   Before this audit numpy, scipy, sympy and tqdm were imported but not declared in
   `pyproject.toml`. They are now declared.
3. Transitive dependencies: the dependency graph was resolved from the installed package metadata
   (`Requires-Dist`, markers evaluated for Linux) in a Python 3.13 virtual environment built from `uv.lock`.
   The result was compared with `uv.lock` (47 entries). The only differences are `colorama` and `tzdata`
   (Windows-only markers) and the project itself.
4. Licenses: read from the installed metadata of each package (`License`, `License-Expression`,
   trove classifiers) and from the license file shipped in the package. The metadata of
   colorama, tzdata, bibtexparser, orjson, tqdm, certifi, pillow, pymatgen and pymatgen-core was
   confirmed on the PyPI JSON API (pypi.org). All other entries rely on the installed metadata.
5. Compatibility with MIT: permissive licenses (MIT, BSD, Apache-2.0, PSF, HPND) are compatible.
   Weak copyleft (MPL-2.0) and dual licenses are flagged and explained.

This list is a technical audit and not legal advice.

## Direct runtime dependencies

| Package | Version | License | Source of license information | MIT compatible | Role | Note |
|---|---|---|---|---|---|---|
| matplotlib | 3.10.8 | PSF-based matplotlib license | LICENSE file; classifier | OK | direct | Permissive. Notice must stay with copies of matplotlib. |
| networkx | 3.6.1 | BSD-3-Clause | METADATA (License-Expression) | OK | direct |  |
| numpy | 2.4.4 | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 | METADATA (License-Expression) | OK | direct | Wheels bundle GCC runtime libraries under GPL-3.0 with GCC exception. See flag F4. |
| Pint | 0.25.3 | BSD-3-Clause | METADATA "BSD"; LICENSE file | OK | direct |  |
| pymatgen | 2026.3.23 | MIT | METADATA; PyPI | OK | direct |  |
| ruamel.yaml | 0.19.1 | MIT | METADATA | OK | direct |  |
| scipy | 1.17.1 | BSD-3-Clause | METADATA text; LICENSE.txt | OK | direct | Wheels bundle GCC runtime libraries under GPL-3.0 with GCC exception. See flag F4. |
| spglib | 2.7.0 | BSD-3-Clause | METADATA (License-Expression) | OK | direct |  |
| sympy | 1.14.0 | BSD-3-Clause | METADATA "BSD"; LICENSE file | OK | direct |  |
| tqdm | 4.67.3 | MPL-2.0 AND MIT | METADATA; LICENCE file | OK | direct | See flag F2. |

## Transitive runtime dependencies

| Package | Version | License | Source of license information | MIT compatible | Role | Note |
|---|---|---|---|---|---|---|
| bibtexparser | 1.4.4 | LGPL-3.0 OR BSD-3-Clause (dual, user choice) | METADATA "License: LGPLv3 or BSD"; COPYING | OK, conditional | transitive runtime (via pymatgen-core) | Use under the BSD option. See flag F1. |
| certifi | 2026.2.25 | MPL-2.0 | METADATA; classifier | OK | transitive runtime (via requests) | See flag F2. |
| charset-normalizer | 3.4.7 | MIT | METADATA | OK | transitive runtime (via requests) |  |
| contourpy | 1.3.3 | BSD-3-Clause | LICENSE file | OK | transitive runtime (via matplotlib) |  |
| cycler | 0.12.1 | BSD-3-Clause | LICENSE file | OK | transitive runtime (via matplotlib) |  |
| flexcache | 0.3 | BSD (3-clause text) | METADATA; LICENSE file | OK | transitive runtime (via pint) |  |
| flexparser | 0.4 | BSD-3-Clause | METADATA | OK | transitive runtime (via pint) |  |
| fonttools | 4.62.1 | MIT | METADATA | OK | transitive runtime (via matplotlib) |  |
| idna | 3.11 | BSD-3-Clause | METADATA (License-Expression) | OK | transitive runtime (via requests) |  |
| joblib | 1.5.3 | BSD-3-Clause | METADATA (License-Expression) | OK | transitive runtime (via pymatgen-core) |  |
| kiwisolver | 1.5.0 | BSD-3-Clause | LICENSE file | OK | transitive runtime (via matplotlib) |  |
| monty | 2026.2.18 | MIT | classifier; LICENSE.rst | OK | transitive runtime (via pymatgen-core) |  |
| mpmath | 1.3.0 | BSD-3-Clause | METADATA "BSD"; LICENSE file | OK | transitive runtime (via sympy) |  |
| narwhals | 2.19.0 | MIT | METADATA | OK | transitive runtime (via plotly) |  |
| orjson | 3.11.8 | MPL-2.0 AND (Apache-2.0 OR MIT) | METADATA (License-Expression) | OK | transitive runtime (via pymatgen-core) | See flag F2. |
| packaging | 26.1 | Apache-2.0 OR BSD-2-Clause | METADATA (License-Expression) | OK | transitive runtime (via matplotlib, plotly, pytest) |  |
| palettable | 3.3.3 | MIT | classifier; METADATA text | OK | transitive runtime (via pymatgen-core) | Bundles ColorBrewer color data under its own license file. |
| pandas | 3.0.2 | BSD-3-Clause | METADATA text; LICENSE file | OK | transitive runtime (via pymatgen-core) |  |
| pillow | 12.2.0 | MIT-CMU (HPND) | METADATA (License-Expression) | OK | transitive runtime (via matplotlib) | Wheels bundle libraries under their own permissive licenses (see its LICENSE file). |
| platformdirs | 4.9.6 | MIT | METADATA | OK | transitive runtime (via pint) |  |
| plotly | 6.7.0 | MIT | METADATA | OK | transitive runtime (via pymatgen-core) |  |
| pymatgen-core | 2026.4.7 | MIT | METADATA; PyPI | OK | transitive runtime (via pymatgen) | Installed with pymatgen. Provides the pymatgen modules used here. |
| pyparsing | 3.3.2 | MIT | METADATA | OK | transitive runtime (via bibtexparser, matplotlib) |  |
| python-dateutil | 2.9.0.post0 | Apache-2.0 OR BSD-3-Clause (dual) | METADATA "Dual License"; classifiers; LICENSE file | OK | transitive runtime (via matplotlib, pandas) | Either license permits MIT redistribution. |
| requests | 2.33.1 | Apache-2.0 | METADATA | OK | transitive runtime (via pymatgen-core) |  |
| six | 1.17.0 | MIT | METADATA | OK | transitive runtime (via python-dateutil) |  |
| tabulate | 0.10.0 | MIT | METADATA | OK | transitive runtime (via pymatgen-core) |  |
| typing_extensions | 4.15.0 | PSF-2.0 | METADATA | OK | transitive runtime (via flexcache, flexparser, pint) |  |
| uncertainties | 3.2.3 | BSD-3-Clause | METADATA "Revised BSD" | OK | transitive runtime (via pymatgen-core) |  |
| urllib3 | 2.6.3 | MIT | METADATA | OK | transitive runtime (via requests) |  |

## Test-only dependencies

Needed only to run `pytest`. They are in the `dev` dependency group and are not needed to use SymDCT.

| Package | Version | License | Source of license information | MIT compatible | Role | Note |
|---|---|---|---|---|---|---|
| iniconfig | 2.3.0 | MIT | METADATA | OK | transitive test (via pytest) |  |
| pluggy | 1.6.0 | MIT | METADATA | OK | transitive test (via pytest) |  |
| Pygments | 2.20.0 | BSD-2-Clause | METADATA | OK | transitive test (via pytest) |  |
| pytest | 9.0.3 | MIT | METADATA | OK | test | Test only. |

## Platform-specific packages (in `uv.lock`, not installed on Linux)

| Package | Version in `uv.lock` | License | Source | MIT compatible | Role |
|---|---|---|---|---|---|
| colorama | 0.4.6 | BSD-3-Clause (copyright Jonathan Hartley and Arnon Yaari) | PyPI project page, https://github.com/tartley/colorama | OK | Windows only, via tqdm and pytest |
| tzdata | 2026.1 | Apache-2.0 | PyPI JSON API | OK | Windows only, via pandas |

The license of colorama is stated as "BSD 3-Clause license" on its PyPI project page.

## Flags

**F1. bibtexparser 1.4.4 (LGPL-3.0 or BSD, dual).**
The package states "The code is distributed under a dual license (at your choice)" in its `COPYING` file.
Choosing the BSD option makes it compatible with MIT. Choosing LGPL would also allow use as a library.
SymDCT only imports it indirectly through pymatgen and does not ship it. No action is required.

**F2. MPL-2.0 packages: certifi, orjson (MPL-2.0 AND Apache-2.0 OR MIT), tqdm (MPL-2.0 AND MIT).**
MPL-2.0 is a file-level (weak) copyleft license. It applies only to modifications of the MPL-licensed
files themselves. It does not extend to code that merely imports or uses the package. SymDCT does not copy
or modify these files. MIT distribution of SymDCT is therefore unaffected. If a future release bundles
these packages (for example in a frozen executable or container image), the MPL-2.0 text and a
pointer to the source must be provided for them.

**F3. Packages with dual licenses: python-dateutil (Apache-2.0 or BSD-3-Clause), packaging (Apache-2.0 or BSD-2-Clause).**
Both options are permissive. No action is required.

**F4. GCC runtime libraries inside binary wheels of numpy and scipy.**
The `LICENSE` files of the numpy and scipy wheels list a bundled GCC runtime library under
"GPL-3.0-or-later WITH GCC-exception-3.1". The GCC Runtime Library Exception allows linking with
non-GPL software. This matters only to those who redistribute the binary wheels. SymDCT does not.

No license was found with status "UNKNOWN". Every package above has a license statement in its metadata.

## Bundled non-package material in this repository

| Item | Content | Origin and license |
|---|---|---|
| `examples/*/input/*.poscar` | Crystal structures (bcc, hcp, fcc model structures) | Prepared by the author. |
| `examples/*/input/*energetics*.yaml` | Energy parameters used in the examples | Literature values and illustrative values. The source is stated in the header of each file (hcp_vac: Ferry et al., J. Nucl. Mater. 524 (2019) 208-215). |
| `examples/*/assets/` | Generated results, plots, Mathematica files | Generated by SymDCT. |
| `examples/figures/*.tex` | TikZ figures | Authored in this project. |

No copied third-party source code, third-party images or third-party data files were found
in the search for copyright, license, author and source notices. The search was a text search
and cannot prove the absence of such material.
