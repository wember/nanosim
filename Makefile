.PHONY: setup activate setup-clean sim sim-i sim-archive sim-test sim-test-clean plot plot-test browse help

# Variables
VENV_NAME = venv
PYTHON = python3
VENV_BIN = $(VENV_NAME)/bin
VENV_PIP = $(VENV_BIN)/pip
VENV_PYTHON = $(VENV_BIN)/python

help:
	@echo "Environment:"
	@echo "  make setup       - Create virtual environment and install dependencies"
	@echo "  make activate    - Print command to activate the virtual environment"
	@echo "  make setup-clean - Remove virtual environment"
	@echo ""
	@echo "Simulation:"
	@echo "  make sim            - Run simulation (archives existing data/ first)"
	@echo "                        With args: make sim ARGS=\"-n 500\""
	@echo "  make sim-i          - Run simulation in interactive mode"
	@echo "  make sim-test       - Run test simulation (cleans test_data/ first)"
	@echo "  make sim-test-clean - Remove test simulation data"
	@echo ""
	@echo "Visualization:"
	@echo "  make browse    - Open web browser to view all runs in data/"
	@echo "  make plot      - Generate plots from last simulation in data/"
	@echo "  make plot-test - Generate plots from test simulation in test_data/"

# ============================================================================
# Environment
# ============================================================================

setup: $(VENV_NAME)/bin/activate

$(VENV_NAME)/bin/activate: requirements.txt
	@echo "Creating virtual environment..."
	$(PYTHON) -m venv $(VENV_NAME)
	@echo "Installing dependencies..."
	$(VENV_PIP) install --upgrade pip
	$(VENV_PIP) install -r requirements.txt
	@echo ""
	@echo "Setup complete! Run 'make activate' to see how to activate the environment."

activate:
	@echo "To activate the virtual environment, run:"
	@echo "  source $(VENV_NAME)/bin/activate"

setup-clean:
	@echo "Removing virtual environment..."
	@rm -rf $(VENV_NAME)

# ============================================================================
# Simulation
# ============================================================================

sim:
	@echo "Running simulation..."
	@$(VENV_PYTHON) creutz-sim/sim.py $(ARGS)

sim-i: sim-archive
	@echo "Running simulation in interactive mode..."
	@$(VENV_PYTHON) creutz-sim/sim.py -i $(ARGS)

sim-test: sim-test-clean
	@echo "Running test simulation..."
	@$(VENV_PYTHON) creutz-sim/sim.py -i --data-dir test_data
	@echo "Flattening test data structure..."
	@if [ -d test_data ]; then \
		TIMESTAMP_DIR=$$(ls -d test_data/*/ 2>/dev/null | head -1); \
		if [ -n "$$TIMESTAMP_DIR" ]; then \
			mv "$$TIMESTAMP_DIR"/* test_data/ 2>/dev/null || true; \
			rm -rf "$$TIMESTAMP_DIR"; \
		fi; \
	fi

plot:
	@echo "Generating plots..."
	@$(VENV_PYTHON) creutz-sim/Sk_comparison.py

plot-test:
	@echo "Generating test plots..."
	@$(VENV_PYTHON) creutz-sim/Sk_comparison.py --data-dir test_data

sim-test-clean:
	@echo "Removing test simulation data..."
	@rm -rf test_data/

browse:
	@echo "Starting archive browser at http://127.0.0.1:5001"
	@(sleep 1.5 && open http://127.0.0.1:5001 2>/dev/null || xdg-open http://127.0.0.1:5001 2>/dev/null) &
	@$(VENV_PYTHON) tools/browse_plots.py
