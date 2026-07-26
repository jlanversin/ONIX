# ONIX Nuclear Depletion Code

ONIX is an open-source depletion code for reactor analysis, fuel-cycle studies, and nuclear archaeology. It is written in Python 3 and can operate either with user-provided one-group data or in coupled mode with [OpenMC](https://github.com/openmc-dev/openmc) for transport-driven depletion.

The full documentation is available on Read the Docs: https://onix-documentation.readthedocs.io/en/latest/index.html

## Core Capabilities

- Standalone depletion with decay, transmutation, and fission-yield data handling.
- Coupled OpenMC depletion with automated exchange of material densities, fluxes, spectra, and one-group reaction data.
- OpenMC-side transport parallelization through the `set_MPI` coupling interface.
- Flexible BUCell selection and per-cell tally nuclide control.
- Support for fixed-power and fixed-flux depletion workflows.
- Reaction-rate ranking, density summaries, and post-processing outputs for burnup analysis.

## MSR And MSRE Workflows

ONIX includes dedicated molten-salt-reactor coupling logic through `Couple_msr`, including support for circulating-fuel loop models built from OpenMC cells.

Supported loop terms include:

- `upstream_cells`
- `elements_removal`
- `nuclide_removal`
- `external_nuclide_addition`
- `holdup`
- `external_salt_stream`
- `recycling`

The repository also includes notebook examples for circulating-fuel studies, including:

- [MSRE.ipynb](examples/MSRE.ipynb)
- [MSFR_fixed_flux.ipynb](examples/MSFR_fixed_flux.ipynb)
- [MSFR_fixed_power.ipynb](examples/MSFR_fixed_power.ipynb)
- [MSFR_with_feed.ipynb](examples/MSFR_with_feed.ipynb)
- [post_processing.ipynb](examples/post_processing.ipynb)

## Installation

Installation instructions are documented here:

https://onix-documentation.readthedocs.io/en/latest/installation.html

## Recent Robustness Updates

This revision fixes two previously hidden MSR edge cases that could affect results in specific user setups:

- `loop_dict` and BUCell ordering edge case:
  Older MSR connectivity assembly relied on list-index lookups that implicitly assumed the `loop_dict` key order, the OpenMC/ONIX BUCell order, and the block layout of the global depletion system were all aligned. If those orders diverged, cross-cell circulation terms could be applied to the wrong cell block.
- Single-letter element removal edge case:
  Older `elements_removal` matching relied on a two-character slice of the nuclide name. That worked for two-letter element symbols but could miss valid one-letter symbols such as `F`, `N`, `O`, or `U`.

These are edge cases rather than general failures:

- The ordering issue matters when the BUCell ordering used by OpenMC/ONIX differs from the order implied by the user loop definition.
- The element-name issue matters when `elements_removal` targets an element with a one-letter symbol.

Other robustness fixes included in this update:

- Mixed ONIX, OpenMC, and ZAMID nuclide identifiers are now normalized consistently in the OpenMC coupling path.
- Nuclide list duplicate detection now works correctly before index dictionaries are built.
- MSR connectivity assembly now validates BUCell passlist ordering before forming the global matrix and checks the matrix/vector shapes before each burn step.
- A legacy indentation issue in `onix/data/script/build_text_mat.py` was corrected so the package compiles cleanly.

## Citing

If you use ONIX in your research, please cite:

- Julien de Troullioud de Lanversin, Moritz Kutt, and Alexander Glaser, "[ONIX: An open-source depletion code](https://doi.org/10.1016/j.anucene.2020.107903)," *Annals of Nuclear Energy* **151** (2021).

If you use the molten-salt-reactor module or the multi-point MSR depletion capability, please also cite:

- M. H. Elhareef and Z. Wu, "[A Novel Multi-Point Depletion Model for Molten Salt Reactors](https://doi.org/10.3390/jne7010017)," *Journal of Nuclear Engineering* **7**(1), 17 (2026).

## Contact

Questions and contributions are welcome. For project inquiries, please contact Julien de Troullioud de Lanversin at `j.detroullioud@gmail.com`.

## Acknowledgments

Special thanks to [Solchan Han](https://github.com/hsc91) for originally identifying the two hidden MSR edge-case bugs and for assisting with their resolution and validation.

## License

ONIX is distributed under the MIT/X license:

https://onix-documentation.readthedocs.io/en/latest/license.html
