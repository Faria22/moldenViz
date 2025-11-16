# Dynamic Grid Plan
## Done
1. Start with low resolution cartesian grid
    This allows the program to tabulate the gtos and mos very quickly
2. Split each mo into its separate lobes
3. For each lobe, create a new grid with the bounds of the original (plus a margin of error) and a very fine grid
4. Merge all grids together and plot
5. Add settings/configuration toggles so users can switch between single and dynamic grids and choose the dynamic grid coordinate system (spherical or cartesian)

## TODO
1. Update documentation
2. Check for errors (when mos are too small)
