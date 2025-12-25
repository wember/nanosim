.PHONY: help setup clean test run-sim run-irr-sim run-sim-small run-irr-sim-small sbatch-sim sbatch-irr-sim plot activate

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

run-sim-small:  ## Run reversible simulation (small: n=1000, s=100, r=3, m=2)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running reversible simulation (small parameters: n=1000, s=100, r=3, m=2)..."
	@./venv/bin/python creutz-sim/sim_small.py

run-irr-sim-small:  ## Run irreversible simulation (small: n=1000, s=100, r=3, m=2)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Running irreversible simulation (small parameters: n=1000, s=100, r=3, m=2)..."
	@./venv/bin/python creutz-sim/irr_sim_small.py

sbatch-sim:  ## Submit reversible simulation to SLURM cluster
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Submitting reversible simulation to SLURM..."
	cd creutz-sim/batch_jobs && sbatch sim_sbatch.sh

sbatch-irr-sim:  ## Submit irreversible simulation to SLURM cluster
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@echo "Submitting irreversible simulation to SLURM..."
	cd creutz-sim/batch_jobs && sbatch irr_sim_sbatch.sh

plot:  ## Run plotting script (requires sim_data.csv)
	@if [ ! -d "venv" ]; then echo "❌ Virtual environment not found. Run 'make setup' first."; exit 1; fi
	@./venv/bin/python creutz-sim/sim_plot.py

activate:  ## Show command to activate virtual environment
	@echo "Run this command to activate the virtual environment:"
	@echo "    source venv/bin/activate"
