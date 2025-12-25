.PHONY: help setup clean test run-sim run-irr-sim run-sim-small run-irr-sim-small run-sim-test run-irr-sim-test sbatch-sim sbatch-irr-sim plot activate run-tests run-tests-serial run-test-file run-examples coverage coverage-html benchmark-jit benchmark-jit-quick

help:  ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Available targets:'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

setup:  ## Create virtual environment and install dependencies
	@./setup.sh

clean:  ## Remove virtual environment and cached files
	@echo "Cleaning up..."
	rm -rf venv/
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	@echo "✓ Cleanup complete"

test:  ## Test that the environment is properly configured
	@echo "Testing environment..."
	@./venv/bin/python -c "import numpy, scipy, pandas, plotly; print('✓ All dependencies imported successfully')"
	@cd creutz-sim && ../venv/bin/python -c "from inferno import Inferno; x = Inferno(100, 5); assert x.E_total == 200, 'Energy conservation failed'; print('✓ Inferno class works correctly')"
	@cd creutz-sim && ../venv/bin/python -c "from irr_inferno import irrInferno; x = irrInferno(100, 5); assert x.E_total == 200, 'Energy conservation failed'; print('✓ irrInferno class works correctly')"
	@echo "✅ All tests passed!"

run-sim:  ## Run reversible simulation (full: n=1000000, s=5000)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running reversible simulation (full parameters)..."
	@./venv/bin/python creutz-sim/sim.py

run-irr-sim:  ## Run irreversible simulation (full: n=1000000, s=5000)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running irreversible simulation (full parameters)..."
	@./venv/bin/python creutz-sim/irr_sim.py

run-parallel-sim:  ## Run parallel reversible simulation (auto-detect cores)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running parallel reversible simulation..."
	@./venv/bin/python creutz-sim/parallel_sim.py

run-parallel-irr-sim:  ## Run parallel irreversible simulation (auto-detect cores)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running parallel irreversible simulation..."
	@./venv/bin/python creutz-sim/parallel_irr_sim.py

run-sim-small:  ## Run reversible simulation (small: n=1000, s=100, r=3, m=2)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@./venv/bin/python creutz-sim/sim.py --n 1000 --s 100 --r 3 --m 2

run-irr-sim-small:  ## Run irreversible simulation (small: n=1000, s=100, r=3, m=2)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@./venv/bin/python creutz-sim/irr_sim.py --n 1000 --s 100 --r 3 --m 2

run-sim-test:  ## Run reversible simulation (test: n=100, s=10, r=2, m=1)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@./venv/bin/python creutz-sim/sim.py --n 100 --s 10 --r 2 --m 1

run-irr-sim-test:  ## Run irreversible simulation (test: n=100, s=10, r=2, m=1)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@./venv/bin/python creutz-sim/irr_sim.py --n 100 --s 10 --r 2 --m 1

run-parallel-sim-test:  ## Run parallel reversible simulation (test: n=100, s=10, r=3, m=6)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@./venv/bin/python creutz-sim/parallel_sim.py --n 100 --s 10 --r 3 --m 6

run-parallel-irr-sim-test:  ## Run parallel irreversible simulation (test: n=100, s=10, r=3, m=6)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@./venv/bin/python creutz-sim/parallel_irr_sim.py --n 100 --s 10 --r 3 --m 6

sbatch-sim:  ## Submit reversible simulation to SLURM cluster
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Submitting reversible simulation to SLURM..."
	cd creutz-sim/batch_jobs && sbatch sim_sbatch.sh

sbatch-irr-sim:  ## Submit irreversible simulation to SLURM cluster
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Submitting irreversible simulation to SLURM..."
	cd creutz-sim/batch_jobs && sbatch irr_sim_sbatch.sh

sbatch-parallel-sim:  ## Submit parallel reversible simulation to SLURM cluster
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Submitting parallel reversible simulation to SLURM..."
	cd creutz-sim/batch_jobs && sbatch parallel_sim_sbatch.sh

sbatch-parallel-irr-sim:  ## Submit parallel irreversible simulation to SLURM cluster
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Submitting parallel irreversible simulation to SLURM..."
	cd creutz-sim/batch_jobs && sbatch parallel_irr_sim_sbatch.sh

plot:  ## Run plotting script (requires sim_data.csv)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@./venv/bin/python creutz-sim/sim_plot.py

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

run-examples:  ## Run all example scripts
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running quick_test.py..."
	@./venv/bin/python examples/quick_test.py
	@echo "\nRunning custom_parameters.py..."
	@./venv/bin/python examples/custom_parameters.py
	@echo "\nRunning analysis_pipeline.py..."
	@./venv/bin/python examples/analysis_pipeline.py

benchmark-jit:  ## Benchmark JIT vs non-JIT performance (usage: make benchmark-jit N=10000 S=100 R=5)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running JIT benchmark..."
	@./venv/bin/python benchmark_jit.py $(N) $(S) $(R)

benchmark-jit-quick:  ## Quick JIT benchmark (N=1000, S=50, R=3)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running quick JIT benchmark..."
	@./venv/bin/python benchmark_jit.py 1000 50 3

activate:  ## Show command to activate virtual environment
	@echo "Run this command to activate the virtual environment:"
	@echo "    source venv/bin/activate"
