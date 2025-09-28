# moldenViz Documentation Improvement Checklist

Use this checklist to guide the next round of documentation updates. Tackle the items in order; each top-level checkbox represents a deliverable.

## 1. Give Readers Context for Molden Files
- [x] Add a "What is a Molden file?" subsection in the Getting Started guide.
  1. Briefly explain the role of the `[Atoms]`, `[GTO]`, and `[MO]` sections and which quantum-chemistry packages export Molden files.
  2. State that moldenViz currently accepts only spherical Gaussian functions and rejects Cartesian (5D/9G) functions.
  3. Link to the official Molden format description: https://www.theochem.ru.nl/molden/molden_format.html.

## 2. Plan Visual Examples
- [x] Reserve image slots in the usage/tutorial pages for GUI walkthroughs.
  1. Insert figure placeholders (e.g., ``.. image:: _static/placeholder.png``) at the points where screenshots should appear.
  2. Add descriptive alt text for each placeholder ("Screenshot of orbital selection panel", "Iso-surface rendering with custom palette", etc.).
  3. Note in the caption which configuration or action the future screenshot must illustrate.

## 3. Expand CLI Coverage
- [x] Build a complete CLI options reference page.
  1. Extract every flag/option from ``moldenViz/cli.py`` and ``moldenViz/__main__.py``.
  2. Present them in a table with the flag, default, and a one-line description.

## 4. Document VTK/CUBE Export (v1.1)
- [x] Add a dedicated subsection on exporting volumes to VTK and Gaussian cube files.
  1. Describe the new v1.1 commands/flags for exporting from the CLI (expected filenames, resolution flags, and how to combine with orbital selection).
  2. Show a Python API example that writes both formats, including required grid/tabulator configuration.
  3. Clarify where outputs are saved, note any dependencies, and list known limitations (e.g., file size, supported grids).
  4. Mention the feature in the release notes or changelog and cross-link from the Getting Started and CLI pages.

## 5. Enrich Configuration Guidance
- [x] Add worked configuration examples for atoms and bonds in the configuration reference.
  1. Provide TOML snippets demonstrating per-atom radii overrides, custom bond colours, and reduced bond radii for large systems.
  2. Describe when to use each override and call out performance implications if relevant.
  3. Cross-link these examples from the troubleshooting section for quick discovery.

## 6. Document Error Handling
- [x] Summarize known parser/tabulator/plotter exceptions in the Troubleshooting section.
  1. List each exception raised for invalid shells, unsupported basis functions, and out-of-bounds orbital indices.
  2. Offer suggested fixes or next steps so users know how to resolve each error.

## 7. Deprecate Old Documentation
- [x] Draft clear retirement steps for legacy Read the Docs pages.
  1. Document the deprecation process in ``deprecating_olds_docs.md`` (banner text, redirect configuration, and owner responsible).
  2. Add numbered instructions for removing or redirecting the v0.3.x Usage page.
  3. Confirm the README points users only to the maintained docs once redirects are in place.

## 8. Highlight Project Roadmap and Contributions
- [x] Create a "Contributing" section in the docs and refresh the README.
  1. In the docs, cover how to file issues, submit PRs, run tests (``pytest``, ``ruff``, ``basedpyright``), and build docs locally.
  2. Move the roadmap detail from the README into a new "Roadmap" page and link to it from both README and docs index.
  3. Leave a concise roadmap summary in the README (high-level bullet list with link to the full page).

## 9. Reuse README Content Thoughtfully
- [x] Align README and docs to avoid duplication.
  1. Transfer installation, tkinter setup, and CLI quickstart instructions from README to the docs' Getting Started page.
  2. Trim the README to a short overview, installation snippet, and link to detailed docs.
  3. Verify that every README section has a richer counterpart in the docs or is intentionally unique.

## 10. Enhance API Usage Examples
- [x] Enrich the Parser/Tabulator/Plotter API pages with runnable snippets.
  1. Add code blocks that create Cartesian grids, tabulate orbitals, and feed them into plotters, noting RAM/performance caveats.
  2. Show both spherical and Cartesian grid workflows and explain when to choose each.
  3. Include notes about sharing tabulators between plotters and reusing computed grids.
