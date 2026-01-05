# Changelog

All notable changes to **TEffectBayes** will be documented in this file.

## [1.1.0] - 2026-01-05

### Added
- User-configurable discretization parameters for Bayesian Network Modeling:
  - `--bnm_n_bins` to control the number of bins in KBinsDiscretizer (default: 2)
  - `--bnm_bin_strategy` to select discretization strategy (`uniform`, `quantile`, `kmeans`)
- Discretization parameters are now exposed at the pipeline level and passed explicitly to the `BNM_CALCULATION` module.

### Changed
- Bayesian Network discretization step no longer uses hard-coded parameters.
- Improved reproducibility and transparency of Bayesian model learning.
- test.config has been updated and both main.nf file and the nf files of some of the module errors based on test data is corrected.

## [1.0.0] - 2025-07-02
### Added
- Initial release of TEffectBayes pipeline.
- Integrated RNA-seq, ChIP-seq, and repeat expression analysis.
- Bayesian Network modeling of gene–repeat–histone interactions.
