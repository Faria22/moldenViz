# Dynamic Grid Plan
## Done
1. Start with low resolution cartesian grid
    This allows the program to tabulate the gtos and mos very quickly
2. Split each mo into its separate lobes
3. For each lobe, create a new grid with the bounds of the original (plus a margin of error) and a very fine grid
4. Merge all grids together and plot

## TODO
1. See if spherical grid gives better results with less points
2. Update settings screen and configuration file options
3. Update documentation
