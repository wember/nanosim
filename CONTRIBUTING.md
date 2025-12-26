# Contributing to Nanosim

Thank you for your interest in contributing to Nanosim! This document provides guidelines for contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How to Contribute](#how-to-contribute)
- [Development Setup](#development-setup)
- [Testing](#testing)
- [Style Guidelines](#style-guidelines)
- [Submitting Changes](#submitting-changes)

## Code of Conduct

This project adheres to the Contributor Covenant [Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/nanosim.git
   cd nanosim
   ```
3. **Set up the development environment**:
   ```bash
   ./setup.sh
   make test-env
   ```

## How to Contribute

### Reporting Bugs

Before creating bug reports, please check existing issues to avoid duplicates. When creating a bug report, include:

- **Clear, descriptive title**
- **Detailed description** of the issue
- **Steps to reproduce** the behavior
- **Expected vs actual behavior**
- **Environment details** (OS, Python version, dependencies)
- **Relevant code samples or error messages**

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When creating an enhancement suggestion, include:

- **Clear, descriptive title**
- **Detailed description** of the proposed functionality
- **Use cases** and examples
- **Why this enhancement would be useful** to most users

### Pull Requests

1. **Create a topic branch** from `main`:

   ```bash
   git checkout -b feature/my-new-feature
   ```

2. **Make your changes**:

   - Follow the [style guidelines](#style-guidelines)
   - Add tests for new functionality
   - Update documentation as needed

3. **Test your changes**:

   ```bash
   make run-tests
   ```

4. **Commit your changes**:

   ```bash
   git commit -m "Add feature: brief description"
   ```

   - Use clear, descriptive commit messages
   - Reference issue numbers when applicable

5. **Push to your fork**:

   ```bash
   git push origin feature/my-new-feature
   ```

6. **Open a Pull Request** on GitHub

## Development Setup

### Environment Setup

```bash
# Create and activate virtual environment
./setup.sh

# Run tests
make run-tests

# Run specific test file
make run-test-file FILE=test_inferno.py

# Run with coverage
make coverage
```

### Project Structure

```
nanosim/
├── creutz-sim/          # Core simulation code
│   ├── inferno.py       # Reversible simulation
│   ├── inferno_irr.py   # Irreversible simulation
│   ├── jit_inferno.py   # JIT-optimized reversible
│   ├── jit_inferno_irr.py # JIT-optimized irreversible
│   ├── jit_functions.py # JIT-compiled functions
│   ├── sim_utils.py     # Shared utilities
│   └── parallel_sim*.py # Parallel execution
├── tests/               # Test suite
├── tools/               # Development tools
├── examples/            # Example scripts
└── docs/               # Documentation
```

## Testing

### Running Tests

```bash
# Full test suite (parallel, ~3x faster)
make run-tests

# Serial execution (for debugging)
make run-tests-serial

# With coverage report
make coverage
```

### Writing Tests

- Place tests in `tests/` directory
- Name test files `test_*.py`
- Use descriptive test names: `test_energy_conservation_during_sweep()`
- Include docstrings explaining what the test validates
- Test both reversible and irreversible versions when applicable

Example:

```python
def test_energy_conservation():
    """Verify total energy remains constant during simulation."""
    sim = Inferno(N=100, R=5)
    initial_energy = sim.E_total

    for _ in range(100):
        sim.demon_move()

    assert sim.E_total == initial_energy
```

## Style Guidelines

### Python Code Style

- **Follow PEP 8** for general Python style
- **Use meaningful variable names**:
  - Physics variables: Match paper notation (`N`, `R`, `K`, `U`)
  - Code variables: Descriptive (`lattice`, `bonds`, `energy_sum`)
- **Add docstrings** to all functions and classes
- **Type hints** are encouraged but not required
- **Comments**: Explain _why_, not _what_

### Physics Code Conventions

- **Energy calculations**: Use integer arithmetic only (prevents float drift)
- **Array operations**: Prefer NumPy operations over loops
- **JIT functions**: Keep pure (no side effects except array modifications)
- **Validation**: Use validation modes for debugging, disable for production

### Documentation

- Update README.md for user-facing changes
- Update ARCHITECTURE.md for structural changes
- Update CHANGELOG.md following [Keep a Changelog](https://keepachangelog.com/)
- Add docstrings to new functions/classes

## Submitting Changes

### Pull Request Checklist

Before submitting a PR, ensure:

- [ ] Code follows style guidelines
- [ ] All tests pass (`make run-tests`)
- [ ] New functionality includes tests
- [ ] Documentation is updated
- [ ] CHANGELOG.md is updated (if applicable)
- [ ] Commit messages are clear and descriptive
- [ ] Branch is up to date with `main`

### Review Process

1. Maintainers will review your PR
2. Address any feedback or requested changes
3. Once approved, your PR will be merged
4. Your contribution will be acknowledged in the release notes

## Performance Considerations

When contributing to performance-critical code:

- **Profile before optimizing**: Use `make profile`
- **Benchmark changes**: Use `make benchmark-jit`
- **Consider JIT compilation**: Hot loops benefit from Numba
- **Validate correctness**: Performance means nothing if results are wrong
- **Document trade-offs**: Explain why you chose one approach over another

## Questions?

- **Open an issue** for questions about the codebase
- **Start a discussion** for broader topics
- **Check existing issues** before asking

Thank you for contributing to Nanosim! 🚀
