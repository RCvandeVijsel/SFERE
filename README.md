# SFERE
SFERE (Scale-dependent Feedback Recursion). Mathematical model to explain self-organization of complex spatial patterns from recursion of a scale-dependent feedback at multiple nested scales. First applied to self-organization of branching channel networks in coastal wetlands (van de Vijsel et al. 2023, Nature Communications).

## Files
#### clPy.ComplexChannelPatterns_Fig1f_RUN9.ipynb
This is the main model script. This is the default version of the SFERE model, i.e. the version used to create Fig. 1f in the manuscript (van de Vijsel et al. 2023, Nature Communications). The script can be opened and run in Jupyter Notebook.

#### HydroFunctions_iPy.cl
This script contains functions (e.g., derivatives and boundary conditions) which are called upon in the main model script. Make sure to place this script in the same directory as the main model script.

#### RUN9
A folder named "RUN9" should be created and placed in the same directory as clPy.ComplexChannelPatterns_Fig1f_RUN9.ipynb and HydroFunctions_iPy.cl. The model results are stored in this folder.

## More information
More information about the functioning of the SFERE model can be found in:
- The manuscript (van de Vijsel et al. 2023, Nature Communications)
- Commented lines in the model scripts (clPy.ComplexChannelPatterns_Fig1f_RUN9.ipynb and HydroFunctions_iPy.cl) themselves
- The full dataset and its readme file, underlying the manuscript (van de Vijsel et al. 2023, Nature Communications). This dataset can be found via https://doi.org/10.4121/8d361887-ec02-4472-a8eb-a9d0f3eacfd6.

## References
van de Vijsel, R.C., van Belzen, J., Bouma, T.J., van der Wal, D., Borsje, B.W., Temmerman, S., Cornacchia, L., Gourgue, O., van de Koppel, J. Vegetation controls on channel network complexity in coastal wetlands. Nature Commununications 14, 7158 (2023). https://doi.org/10.1038/s41467-023-42731-3.

![SupplFig2_dpi=250](https://github.com/RCvandeVijsel/SFERE/assets/130892539/54ec43a5-d2a1-48a6-ac18-43fd80ae0f91)
Supplementary Fig. 2 (van de Vijsel et al. 2023) shows channel network development in the default configuration of SFERE (with vegetation, right column), and idem without vegetation dynamics (left column). Vegetation clearly increases the extent and complexity of channel networks.
