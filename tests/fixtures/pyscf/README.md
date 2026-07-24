# PySCF Molden development fixtures

These fixtures are generated from the geometries in
`src/moldenViz/examples/molden_files` as preliminary test data for issues
[#116](https://github.com/Faria22/moldenViz/issues/116) and
[#117](https://github.com/Faria22/moldenViz/issues/117).

The generator uses PySCF 2.14.0 and the `cc-pVQZ` basis. Every bundled geometry
then has shells from `l=0` through `l=4`. It writes both spherical
(`5D`/`7F`/`9G`) and Cartesian (`6D`/`10F`/`15G`) representations. To limit
fixture size, each file contains at most 12 orbitals around the occupied/virtual
frontier, with coefficients spanning every atomic orbital.

The orbitals come from diagonalizing PySCF's one-electron core Hamiltonian
rather than running an SCF calculation. Parser and tabulator fixtures need
deterministic, orthonormal coefficient vectors in the relevant AO layout, not a
converged electronic structure. Avoiding SCF also avoids constructing large
four- or three-index electron-repulsion tensors for the bigger examples.

Generate all files from the repository root:

```console
uv run --with pyscf==2.14.0 python scripts/generate_pyscf_molden.py
```

Generate one molecule or one representation:

```console
uv run --with pyscf==2.14.0 python scripts/generate_pyscf_molden.py co
uv run --with pyscf==2.14.0 python scripts/generate_pyscf_molden.py \
  --representation cartesian h2o
```

Existing files are left untouched unless `--overwrite` is passed.
