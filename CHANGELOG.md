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
- test.config has been updated 
- main.nf of the pipeline and the main.nf files of Preprocessing stage module errors based on the new test.config were corrected.

## [1.0.0] - 2025-07-02
### Added
- Initial release of TEffectBayes pipeline.
- Integrated RNA-seq, ChIP-seq, and repeat expression analysis.
- Bayesian Network modeling of gene–repeat–histone interactions.
