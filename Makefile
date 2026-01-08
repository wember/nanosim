.PHONY: setup activate clean-venv sim plot clean help

# Variables
VENV_NAME = venv
PYTHON = python3
VENV_BIN = $(VENV_NAME)/bin
VENV_PIP = $(VENV_BIN)/pip
VENV_PYTHON = $(VENV_BIN)/python

help:
	@echo "Available targets:"
	@echo ""
	@echo "Environment:"
	@echo "  make setup      - Create virtual environment and install dependencies"
	@echo "  make activate   - Print command to activate the virtual environment"
	@echo "  make clean-venv - Remove virtual environment"
	@echo ""
	@echo "Simulation:"
	@echo "  make sim        - Run simulation"
	@echo "  make plot       - Generate plots from simulation data"
	@echo "  make clean      - Remove simulation data and output files"

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

clean-venv:
	@echo "Removing virtual environment..."
	@rm -rf $(VENV_NAME)

# ============================================================================
# Simulation
# ============================================================================

sim: clean
	@echo "Running simulation..."
	@$(VENV_PYTHON) creutz-sim/sim.py

plot:
	@echo "Generating plots..."
	@$(VENV_PYTHON) creutz-sim/Sk_comparison.py

clean:
	@echo "Removing simulation data..."
	@rm -rf data/ init-fin/
