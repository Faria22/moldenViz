# Review of the moldenViz documentation and recommendations for improvement

## Key gaps and limitations in the documentation
1. Minimal high‑level explanation of the Molden format and limitations.  The docs assume readers already know what a
Molden file is.  There is no explanation of the file format, the meaning of the [Atoms], [GTO] and [MO] sections or
which quantum‑chemistry programs can output Molden files.  The parser currently only supports spherical Gaussian
functions; the code raises an error when cartesian functions (5D or 9G) appear, yet this limitation is mentioned only in
a code comment rather than in the user guide.
- Instructions: Briefly mention what the molden format is, and add a link for more details. The link is: https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.theochem.ru.nl/molden/molden_format.html&ved=2ahUKEwjw6u-tmfqPAxWnQzABHVDCO9EQFnoECB4QAQ&usg=AOvVaw0fShQ_5YiCTzxQCMCzw-7W
2. No diagrams or screenshots.  The documentation does not show what the interactive plot window looks like, what the
orbital selection interface contains or how different configuration settings affect the rendered molecule.  Screenshots
and annotated images would help users understand the visual output and the effect of configuration options.
- Instructions: Add spaces where I can attach screenshots. Make sure to add alt text so that I know what should go into each image.
3. Sparse CLI documentation.  The CLI guide lists only the basic commands and example names.  It does not enumerate all
available flags and options (e.g., output image export, resolution settings or configuration overrides).  Users must run
moldenViz -h to discover other flags; the documentation should instead summarise them in a table.
4. Missing configuration examples for atoms and bonds.  While the configuration reference lists the available keys, it
lacks examples demonstrating how to set per‑atom radii, adjust bond colouring or reduce bond radius for large molecules.
Users might not realise they can override element‑specific colours or maximum bond counts until they read the raw TOML
snippet.
6. Version confusion and outdated pages.  The v0.3.x documentation remains on Read‑the‑Docs but contains no content in
the Usage page.  Users searching for the package might land on these old pages and assume the project is poorly
documented.  The README file includes useful information not mirrored in earlier docs, such as example molecules and
configuration hints.
- Instructions: Please give me more details on how to do this. Add these details to `deprecating_olds_docs.md`
7. No contribution guidelines or roadmap discussion.  Although the GitHub README lists a roadmap with planned features
and completed milestones, the documentation site does not mention these plans or provide guidance on how to contribute
to the project.
- Instructions: Add this to the documentation as well as the README

## Recommendations for improvement ### Consolidate and reorganise documentation
- Redirect or remove outdated versions.  Configure Read‑the‑Docs to redirect the v0.3.x Usage page to the latest version
or add a banner indicating that the page is outdated.  This will prevent confusion for users who see the empty Usage
page.
- Integrate README content into the Getting Started guide.  The README contains crucial information on installation,
tkinter setup, CLI usage and example molecules.  Incorporate these details into the Getting Started and CLI sections of
the documentation rather than referring users back to GitHub.
- Provide a high‑level overview of the Molden format.  Add a new section in the Getting Started or Python API guide
explaining the structure of Molden files, the purpose of the [Atoms], [GTO] and [MO] sections and which
quantum‑chemistry codes generate compatible outputs.  Explicitly state that the current parser supports only spherical
Gaussian basis functions and does not support Cartesian orbitals (e.g., 5D/9G), which are rejected during parsing.
- Add a detailed CLI reference.  Create a table enumerating all command‑line options (flags for selecting orbitals,
exporting images, disabling bonds, specifying configuration files, etc.), with descriptions and default values.  This
reduces the need for users to run moldenViz ‑h to discover features.
- Enhance the API reference with examples.  For each class and key method in the Parser, Tabulator and Plotter APIs,
include code snippets demonstrating typical usage.  For instance, show how to create a Cartesian grid with custom
resolution and tabulate a set of orbitals, then pass the same Tabulator to a Plotter.  Explain performance
considerations (e.g., large grids may consume significant RAM) and how to choose between spherical and Cartesian grids.
- Document configuration options comprehensively.  Build upon the configuration reference by adding narrative
explanations and examples.  For example:
- Show how to customise bond colours and radii to emphasise specific bonds, with before/after screenshots.
- Explain the difference between color_type='uniform' and color_type='split' for bonds.
- Provide examples of per‑atom overrides using Atom.<atomic number> tables to change colour and radius.
- Describe how the min_radius and max_radius_multiplier settings influence the default spherical grid and how increasing
num_r_points or num_x_points affects the resolution and computational cost.
- Show how to adjust molecular orbital contour levels and opacities in the configuration and illustrate how these
changes affect the visualisation.

### Include visuals and interactive demonstrations
- Screenshots of the GUI and interactive controls.  Add annotated images showing the PyVista viewer, the Tkinter orbital
selection window, and the effect of toggling the molecule visibility or adjusting contour levels.  Visuals can accompany
each major section (CLI, Python API, configuration reference) to orient new users.
- Demonstrate configuration effects.  Provide side‑by‑side images of a molecule rendered with default settings versus
with customised colours, bond radii and orbital contour levels.  This will help users understand the purpose of each
configuration key.
- Embed short videos or GIFs.  For features such as rotating the molecule, selecting orbitals, and changing isosurface
levels, consider embedding short animations (hosted on an external platform if necessary) to show the interactivity of
the application.

### Expand guidance for advanced workflows and performance
- Explain how to handle large systems.  Provide guidelines on selecting appropriate grid resolutions and memory
considerations when dealing with large molecules or high numbers of molecular orbitals.  Discuss how to sample only a
subset of orbitals to reduce computation.
- Headless or remote visualisation.  Document whether moldenViz can be used in headless environments (e.g., HPC clusters
or remote servers without a display).  If not currently supported, provide instructions for exporting the grid and using
external tools (such as VMD or PyVista offline rendering) to generate static images.
- Extend Tabulator to support custom grids.  The current implementation restricts the grid type to spherical or
Cartesian and raises an error otherwise.  Consider adding support for custom grids (e.g., radial grids with non‑uniform
spacing or adaptive grids) and document how to use them.
- Clarify how to supply pre‑computed grid data.  The advanced workflow section touches on passing a Tabulator with a
custom grid to Plotter, but it does not show how to share grid data between multiple calculations or how to reuse a grid
for different molecules.  Provide guidelines on serialising grid and orbital data for reuse.
- Document the v1.1 export feature.  Add a dedicated subsection in the user guide explaining the new export workflow,
covering CLI flags, Python API usage, supported formats and a short example so users can find and apply the feature.
- Document error handling and exceptions.  Beyond the Troubleshooting page, include a comprehensive list of potential
exceptions thrown by the parser, tabulator and plotter (e.g., invalid shell labels, unsupported basis functions,
out‑of‑bounds orbital indices).  Explain how to catch these errors and suggest corrective actions.

### Encourage contributions and community involvement
- Add a “Contributing” section.  Describe how users can contribute bug reports, feature requests or code contributions.
Provide a link to the GitHub repository and instruct potential contributors on how to run tests and build the
documentation locally.
- Maintain an updated roadmap.  Move the roadmap currently in the README to a “Roadmap” page in the documentation.
Include a changelog summarising major changes between versions and planned features (e.g., support for Cartesian basis
functions or new visualisation features).  This helps users understand the project’s evolution and future direction.
- Gather and showcase user examples.  Encourage users to share screenshots or notebooks using moldenViz and host a
gallery of interesting molecules and orbitals rendered with the tool.  Real‑world examples inspire new users and
demonstrate the capabilities of the software.

## Conclusion The 1.0.0 documentation for moldenViz represents a significant improvement over earlier versions,
providing a structured user guide, configuration reference and troubleshooting section.  However, there remain several
areas that could be enhanced.  The documentation should offer more context about the Molden format and the mathematical
foundations of the visualisation, include comprehensive CLI and API references with examples, and provide visuals to
demonstrate the tool’s capabilities.  Clearer guidance on advanced usage, performance considerations, citation, and
community involvement will make the package more accessible to scientists and developers.  Implementing the
recommendations outlined above will help users adopt moldenViz confidently and contribute to its growth.
