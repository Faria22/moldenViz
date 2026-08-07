# PySCF Molden development fixtures

These fixtures are generated from selected geometries in
`src/moldenViz/examples/molden_files` for
[#116](https://github.com/Faria22/moldenViz/issues/116).

The committed CO, water, and benzene fixtures provide small, medium, and large
spherical-basis (`5D`/`7F`/`9G`) cases for parser tests and high-angular-momentum
benchmarks. The generator uses PySCF 2.14.0 and the `cc-pVQZ` basis, so every
fixture has shells from `l=0` through `l=4`. To limit fixture size, each file
contains at most 12 orbitals around the occupied/virtual frontier, with
coefficients spanning every atomic orbital.

The orbitals come from diagonalizing PySCF's one-electron core Hamiltonian
rather than running an SCF calculation. Parser and tabulator fixtures need
deterministic, orthonormal coefficient vectors in the relevant AO layout, not a
converged electronic structure. Avoiding SCF also avoids constructing large
four- or three-index electron-repulsion tensors for the bigger examples.

Regenerate the three committed fixtures from the repository root:

```console
uv run --with pyscf==2.14.0 python scripts/generate_pyscf_molden.py
```

Generate one molecule, every bundled molecule, or development-only Cartesian
data explicitly:

```console
uv run --with pyscf==2.14.0 python scripts/generate_pyscf_molden.py co
uv run --with pyscf==2.14.0 python scripts/generate_pyscf_molden.py \
  all
uv run --with pyscf==2.14.0 python scripts/generate_pyscf_molden.py \
  --representation cartesian h2o
```

Existing files are left untouched unless `--overwrite` is passed.
