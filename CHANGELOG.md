# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Comprehensive Test Coverage** (2026-01-02): Achieved 84% coverage with 337 tests

  - JIT functions: 28 comprehensive behavioral tests (coverage tool limitation documented)
  - sim_utils.py: 100% coverage with exception handling tests
  - parallel_sim.py & parallel_sim_irr.py: 65% coverage (35% is unmeasurable **main** blocks)
  - Visualization files: Main function tests added (99%/98% coverage)
  - Performance test stability: Fixed intermittent timing-based test failures
  - Test suite reliability: All tests pass consistently with proper thresholds

- **Numba JIT Compilation** (2025-01-18): ~3,881x speedup via machine code compilation
  - New files: `jit_functions.py`, `jit_inferno.py`, `jit_inferno_irr.py`
  - Benchmark script: `benchmark_jit.py` with comprehensive performance testing
  - Test suite: `test_jit_implementation.py` with 22 comprehensive tests
  - Combined with parallel processing: ~8,900x total speedup
  - Production impact: 2.48 hours → 1.0 second
  - Documentation: `JIT_OPTIMIZATION_SUMMARY.md`, `JIT_TEST_SUMMARY.md`
  - Makefile targets: `make benchmark-jit` (simplified from two targets)
  - New dependency: `numba>=0.63.1`
- **Parallel Processing** (2025-01-17): 13-14x speedup with multiprocessing
  - New files: `parallel_sim.py`, `parallel_sim_irr.py`
  - Auto-detects CPU cores, supports SLURM integration
  - 21 new tests in `test_parallel_execution.py`
  - Makefile targets: `make run-sim`, `make run-irr-sim` (parallel JIT by default)
- **Random Sign Generation Optimization** (2025-01-16): 1.7x speedup
  - Replaced array generation with `np.random.choice([-1, 1], size=N)`
  - Applied to both Inferno and irrInferno classes
- **Index Cycling Optimization** (2025-01-16): Replaced modulo with conditional reset
  - Branch prediction friendly (99.9% correct speculation)
  - Applied to 6 indices in Inferno, 2 in irrInferno
- **Pre-computed Neighbor Arrays** (2025-01-16): Eliminated modulo in hot path
  - `right_neighbor` and `left_neighbor` arrays computed at initialization
  - Direct array lookups replace modulo operations
- Command-line argument parsing with argparse for flexible parameter configuration
- Progress bars using tqdm for better user feedback during long simulations
- Structured logging system with file and console outputs
- JSON metadata output for each simulation run (includes timestamp, parameters, simulation type)
- Comprehensive docstrings with type hints for all major functions
- Unit tests using pytest (337 tests total with 84% coverage)
- Performance benchmarks and stress tests (large lattices, long runs, validation mechanisms)
- Examples directory with sample scripts including benchmark_validation.py
- Type hints throughout codebase for better IDE support
- Makefile targets for testing, benchmarking, and running simulations
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
- Updated `README.md` with performance metrics and JIT usage
- Updated `OPTIMIZATIONS.md` with Section 8 (Numba JIT) details
- Updated `requirements.txt` to include numba dependency

### Removed

- Hostname-based path switching (replaced with portable relative paths)
- Old progress indicator print statements
- `sim_small.py` and `irr_sim_small.py` (redundant with command-line args)
- Deprecated `my_venv/` directory (use `venv/` instead)

### Fixed

- Path portability issues across different machines and environments
- Consistent use of `os.path.join()` for cross-platform compatibility
- Energy conservation verified at machine precision (no drift)

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
