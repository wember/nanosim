.PHONY: help setup clean activate test-env compile run-sim run-irr-sim run-sim-small run-irr-sim-small sbatch-sim sbatch-irr-sim plot plot-radii plot-comparison run-examples benchmark-jit profile profile-inferno profile-irr view-profile run-tests run-tests-serial run-test-file coverage coverage-html

help:  ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@printf '\033[1mEnvironment Setup:\033[0m\n'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; /^(setup|clean|activate|test-env|compile):/ {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ''
	@printf '\033[1mSimulations:\033[0m\n'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; /^(run-sim|run-irr-sim|run-sim-small|run-irr-sim-small):/ {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ''
	@printf '\033[1mHPC / SLURM:\033[0m\n'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; /^(sbatch-sim|sbatch-irr-sim):/ {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ''
	@printf '\033[1mAnalysis & Visualization:\033[0m\n'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; /^(plot|plot-radii|plot-comparison|run-examples):/ {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ''
	@printf '\033[1mPerformance & Profiling:\033[0m\n'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; /^(benchmark-jit|profile|profile-inferno|profile-irr|view-profile):/ {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ''
	@printf '\033[1mTesting & Coverage:\033[0m\n'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; /^(run-tests|run-tests-serial|run-test-file|coverage|coverage-html):/ {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ''
	@echo 'Notes:'
	@echo '  - All simulations use parallel JIT by default (fastest: ~1400x speedup)'
	@echo '  - Use --cores 1 for single-core execution if needed'
	@echo '  - Custom params: make run-sim ARGS="--n 100000 --s 5000 --r 8"'

# =============================================================================
# Environment Setup
# =============================================================================

setup:  ## Create virtual environment, install dependencies, and compile JIT functions
	@./setup.sh
	@echo ""
	@echo "Compiling JIT functions (one-time setup, ~15-20 seconds)..."
	@make compile

clean:  ## Remove virtual environment and cached files
	@echo "Cleaning up..."
	rm -rf venv/
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	@echo "✓ Cleanup complete"

activate:  ## Show command to activate virtual environment
	@echo "Run this command to activate the virtual environment:"
	@echo "    source venv/bin/activate"

test-env:  ## Test that the environment is properly configured
	@echo "Testing environment..."
	@./venv/bin/python -c "import numpy, scipy, pandas, plotly; print('✓ All dependencies imported successfully')"
	@cd creutz-sim && ../venv/bin/python -c "from inferno import Inferno; x = Inferno(100, 5); assert x.E_total == 200, 'Energy conservation failed'; print('✓ Inferno class works correctly')"
	@cd creutz-sim && ../venv/bin/python -c "from irr_inferno import irrInferno; x = irrInferno(100, 5); assert x.E_total == 200, 'Energy conservation failed'; print('✓ irrInferno class works correctly')"
	@echo "✅ All tests passed!"

compile:  ## Warm up JIT compilation (run once after installation)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@./venv/bin/python creutz-sim/parallel_sim.py --jit --n 100 --s 5 --r 3 --m 2 > /dev/null 2>&1
	@echo "✓ JIT compilation complete (cached for future runs)"

# =============================================================================
# Simulations
# =============================================================================

# =============================================================================
# Simulations
# =============================================================================

run-sim:  ## Run reversible simulation with parallel JIT (fastest, use ARGS="--cores 1" for single-core)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running reversible simulation (parallel JIT, ~1400x faster)..."
	@./venv/bin/python creutz-sim/parallel_sim.py --jit $(ARGS)

run-irr-sim:  ## Run irreversible simulation with parallel JIT (fastest, use ARGS="--cores 1" for single-core)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running irreversible simulation (parallel JIT, ~1400x faster)..."
	@./venv/bin/python creutz-sim/parallel_irr_sim.py --jit $(ARGS)

run-sim-small:  ## Quick small-scale reversible simulation (n=100, s=10, r=3, m=6)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@./venv/bin/python creutz-sim/parallel_sim.py --jit --n 100 --s 10 --r 3 --m 6

run-irr-sim-small:  ## Quick small-scale irreversible simulation (n=100, s=10, r=3, m=6)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@./venv/bin/python creutz-sim/parallel_irr_sim.py --jit --n 100 --s 10 --r 3 --m 6

# =============================================================================
# HPC / SLURM
# =============================================================================

sbatch-sim:  ## Submit reversible simulation to SLURM (parallel JIT, fastest)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Submitting reversible simulation to SLURM (parallel JIT)..."
	cd creutz-sim/batch_jobs && sbatch sim_sbatch.sh

sbatch-irr-sim:  ## Submit irreversible simulation to SLURM (parallel JIT, fastest)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Submitting irreversible simulation to SLURM (parallel JIT)..."
	cd creutz-sim/batch_jobs && sbatch irr_sim_sbatch.sh

# =============================================================================
# Analysis & Visualization
# =============================================================================

plot:  ## Run plotting script (requires simulation CSV output)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Note: Run a simulation first to generate CSV data (e.g., make run-sim-small)"
	@./venv/bin/python creutz-sim/sim_plot.py

plot-radii:  ## Plot results across all radii for single simulation type (reversible or irreversible)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Generating multi-radius comparison plot..."
	@echo "Note: Requires data in data/r{0-10}/ directories"
	@./venv/bin/python creutz-sim/sim_plot_r.py

plot-comparison:  ## Compare entropy between reversible and irreversible simulations (requires data from both)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Generating entropy comparison plot (reversible vs irreversible)..."
	@echo "Note: Requires data in data/r{0-10}/ and data/irr/r{0-10}/ directories"
	@./venv/bin/python creutz-sim/Sk_comparison.py

run-examples:  ## Run all example scripts
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running quick_test.py..."
	@./venv/bin/python examples/quick_test.py
	@echo "\nRunning custom_parameters.py..."
	@./venv/bin/python examples/custom_parameters.py
	@echo "\nRunning analysis_pipeline.py..."
	@./venv/bin/python examples/analysis_pipeline.py

# =============================================================================
# Performance & Profiling
# =============================================================================

benchmark-jit:  ## Benchmark JIT vs non-JIT performance (shows 70-107x speedup)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running JIT benchmark (N=10000, S=100, R=5)..."
	@./venv/bin/python benchmark_jit.py

profile:  ## Profile simulation to identify bottlenecks (usage: make profile MODE=inferno N=10000 S=100 R=5)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@mkdir -p profiling
	@echo "Profiling $(or $(MODE),inferno) simulation..."
	@./venv/bin/python profile_sim.py --mode $(or $(MODE),inferno) --n $(or $(N),10000) --s $(or $(S),100) --r $(or $(R),5)

profile-inferno:  ## Profile Inferno (reversible) simulation
	@make profile MODE=inferno N=10000 S=100 R=5

profile-irr:  ## Profile irrInferno (irreversible) simulation
	@make profile MODE=irr_inferno N=10000 S=100 R=5

view-profile:  ## View most recent profiling results (usage: make view-profile FILE=profile_inferno_n10000_s100.stats)
	@if [ -z "$(FILE)" ]; then \
		echo "Available profile files:"; \
		ls -lh profiling/*.stats 2>/dev/null || echo "No profile files found"; \
	else \
		python -c "import pstats; p = pstats.Stats('profiling/$(FILE)'); p.sort_stats('cumulative').print_stats(30)"; \
	fi

# =============================================================================
# Testing & Coverage
# =============================================================================

run-tests:  ## Run unit tests in parallel (pass ARGS="..." for custom options)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running unit tests in parallel..."
	@./venv/bin/python -m pytest tests/ -v -n auto $(ARGS)

run-tests-serial:  ## Run unit tests serially (for debugging race conditions)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running unit tests serially..."
	@./venv/bin/python -m pytest tests/ -v $(ARGS)

run-test-file:  ## Run tests from a single file (usage: make run-test-file FILE=test_inferno.py)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@if [ -z "$(FILE)" ]; then \
		echo "❌ Error: FILE variable not set. Usage: make run-test-file FILE=test_inferno.py"; \
		exit 1; \
	fi
	@echo "Running tests from $(FILE)..."
	@./venv/bin/python -m pytest tests/$(FILE) -v

coverage:  ## Run tests with coverage report (parallel execution)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running tests with coverage analysis (parallel)..."
	@COVERAGE_FILE=/tmp/.coverage ./venv/bin/python -m pytest tests/ --cov=creutz-sim --cov-report=term-missing --cov-report=html -n auto $(ARGS)
	@rm -f .coverage
	@echo ""
	@echo "✓ Coverage report generated at htmlcov/index.html"

coverage-html:  ## Run coverage and open HTML report in browser
	@make coverage
	@echo "Opening coverage report in browser..."
	@open htmlcov/index.html || xdg-open htmlcov/index.html 2>/dev/null || echo "Please open htmlcov/index.html manually"
