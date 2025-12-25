# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Command-line argument parsing with argparse for flexible parameter configuration
- Progress bars using tqdm for better user feedback during long simulations
- Structured logging system with file and console outputs
- JSON metadata output for each simulation run (includes timestamp, parameters, simulation type)
- Comprehensive docstrings with type hints for all major functions
- Unit tests using pytest (38 tests across test_inferno.py, test_irr_inferno.py, test_performance.py)
- Performance benchmarks and stress tests (large lattices, long runs, validation mechanisms)
- Examples directory with sample scripts including benchmark_validation.py
- Type hints throughout codebase for better IDE support
- Makefile targets for testing
- Automatic directory creation for data output
- Logs directory for simulation logs
- **Performance optimization**: Configurable validation mode (`--validate` option)
  - `off` (default): No automatic validation, fastest performance
  - `periodic`: Validation every 100 sweeps for testing
  - `frequent`: Validation every sweep for debugging
- Validation overhead reduction: automatic checks are now optional

### Changed

- Replaced hardcoded absolute paths with relative paths from project root
- Replaced print statements with proper logging
- Simulation parameters now configurable via command-line (--n, --s, --r, --m)
- Default parameters remain in code but can be overridden
- Improved error handling with structured exception logging
- Energy conservation checks now logged rather than printed

### Removed

- Hostname-based path switching (replaced with portable relative paths)
- Old progress indicator print statements
- Redundant `sim_small.py` and `irr_sim_small.py` wrapper scripts

### Fixed

- Path portability issues across different machines and environments
- Consistent use of `os.path.join()` for cross-platform compatibility

## [1.0.0] - 2025-12-25

### Added

- Initial release
- Reversible simulation (Inferno class)
- Irreversible simulation (irrInferno class)
- Basic visualization scripts (sim_plot.py, Sk_comparison.py)
- SLURM batch job scripts for HPC clusters
- Virtual environment setup with requirements.txt
- Makefile for task automation
- Comprehensive README and ARCHITECTURE documentation
