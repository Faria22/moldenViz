Performance Benchmarks
======================

The performance suite uses `Airspeed Velocity (ASV) <https://asv.readthedocs.io/>`_
instead of pytest. Keeping correctness tests and performance measurements separate
lets ordinary test and coverage runs work without benchmark plugins while ASV
retains comparable timing, peak-memory, commit, environment, and machine metadata.

Benchmark Matrix
----------------

The suite covers each bundled example. The molecule sizes below provide context
for scaling results:

.. list-table::
   :header-rows: 1

   * - Molecule
     - Atoms
     - Basis functions
     - Molecular orbitals
   * - acrolein
     - 8
     - 64
     - 64
   * - benzene
     - 12
     - 96
     - 96
   * - CO
     - 2
     - 26
     - 26
   * - CO2
     - 3
     - 39
     - 39
   * - furan
     - 9
     - 77
     - 77
   * - H2O
     - 3
     - 19
     - 19
   * - O2
     - 2
     - 26
     - 26
   * - prismane
     - 12
     - 96
     - 96
   * - pyridine
     - 11
     - 93
     - 93

Runtime coverage includes:

* parsing every bundled molecule;
* Cartesian grid construction at edge sizes 10, 25, 50, and 100;
* GTO tabulation for every molecule on 10³, 25³, 50³, and 100³ point grids;
* sequential and four-worker GTO tabulation for H2O, furan, and benzene on
  50³ and 100³ point grids;
* single, five, and all-MO contractions for every molecule on those grids; and
* Cartesian real solid harmonics for 10,000, 125,000, and 1,000,000 points
  through angular momenta 0, 2, and 4.

Peak resident memory is measured for H2O, furan, and benzene at 50³ and 75³
points. These cases represent small, medium, and large bundled basis sets and
cover both GTO tabulation and all-MO contraction. A focused benzene 75³ case also
compares one and four GTO workers so the concurrency-dependent memory
amplification remains visible. ASV's peak-memory measurement includes the process
and setup baseline, so compare results from the same machine and environment
rather than treating them as allocation-only byte counts.

The default GTO policy uses up to four workers through 125,000 points and switches
to sequential tabulation for larger grids. Worker-scaling benchmarks pass an
explicit limit to compare the runtime and memory tradeoff on both sides of that
threshold.

Authoritative Environment
-------------------------

Routine benchmark results use CPython 3.13 and NumPy 2.2 in an ASV-managed
``virtualenv``. This single environment keeps the full molecule and grid matrix
tractable and makes commit-to-commit results comparable. The regular CI matrix
continues to test every supported Python version independently.

Install the development dependencies, then validate benchmark discovery and run
a short isolated-environment smoke test:

.. code-block:: bash

   just sync
   just bench-smoke

The smoke command runs ``asv check`` and quick grid-construction cases. It verifies
configuration, discovery, project building, isolated installation, and benchmark
execution, but its one-shot timings are not performance evidence. It also records
ASV's detected metadata for the current machine in ``~/.asv-machine.json`` if that
machine has not run ASV before.

Running and Comparing Results
-----------------------------

Run the complete suite for the current commit:

.. code-block:: bash

   just bench

Pass an ASV revision and optional benchmark filters through the Just recipe when a
focused run is more useful:

.. code-block:: bash

   just bench main..HEAD
   just bench HEAD^! --bench TimeGTOTabulation

For a statistically sampled side-by-side comparison between a base revision and
the current branch, use:

.. code-block:: bash

   uv run --locked --group dev asv continuous --factor 1.1 origin/main HEAD

The default 10% reporting factor is a triage signal, not an automatic correctness
threshold. Re-run suspicious cases on an idle, stable machine and inspect raw
samples before concluding that a change regressed. Peak-memory changes use a
wider 20% publication threshold because process RSS is coarser and more
platform-dependent.

Results, Metadata, and Retention
--------------------------------

ASV writes environments, raw JSON results, and the generated site beneath the
ignored ``.asv/`` directory:

* ``.asv/env/`` contains reusable isolated benchmark environments;
* ``.asv/results/`` stores results by machine, commit, Python, and dependency
  metadata; and
* ``.asv/html/`` contains the published static report.

Inspect saved measurements on the command line or build and preview the site:

.. code-block:: bash

   uv run --locked --group dev asv show HEAD
   uv run --locked --group dev asv publish
   uv run --locked --group dev asv preview

Local results remain untracked because comparisons across unrelated machines are
usually misleading. Results intended as a durable baseline should come from a
named, dedicated runner and retain ``.asv/results/`` as a versioned benchmark-data
branch or CI artifact together with the machine metadata. Do not combine results
from different machine names into a regression decision.

CI Strategy
-----------

Pull-request CI runs ``just bench-smoke`` on the authoritative Python
version. This catches broken configurations, imports, builds, isolated installs,
and benchmark calls without asserting noisy timing limits on shared GitHub
runners.

Meaningful regression monitoring should run the full ``asv continuous`` comparison
on dedicated hardware with fixed power settings and minimal background load.
Treat the configured 10% runtime and 20% peak-memory factors as review triggers,
retain the raw samples, and confirm a reported regression with another run before
blocking or reverting a change.
