# SFERE
SFERE (Scale-dependent Feedback Recursion). Mathematical model to explain self-organization of complex spatial patterns from recursion of a scale-dependent feedback at multiple nested scales. First applied to self-organization of branching channel networks in coastal wetlands (van de Vijsel et al. 2023, Nature Communications).

## Files
### clPy.ComplexChannelPatterns_Fig1f_RUN9.ipynb
This is the main model script. This is the default version of the SFERE model, i.e. the version used to create Fig. 1f in the manuscript (van de Vijsel et al. 2023, Nature Communications). The script can be opened and run in Jupyter Notebook.

### HydroFunctions_iPy.cl
This script contains functions (e.g., derivatives and boundary conditions) which are called upon in the main model script. Make sure to place this script in the same directory as the main model script.

### RUN9
A folder named "RUN9" should be created and placed in the same directory as clPy.ComplexChannelPatterns_Fig1f_RUN9.ipynb and HydroFunctions_iPy.cl. The model results are stored in this folder.

## References
van de Vijsel, R.C., van Belzen, J., Bouma, T.J., van der Wal, D., Borsje, B.W., Temmerman, S., Cornacchia, L., Gourgue, O., van de Koppel, J. (2023, Nature Communications). Vegetation controls on channel network complexity in coastal wetlands.
