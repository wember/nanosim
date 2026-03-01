.PHONY: setup activate setup-clean sim sim-i sim-repeat sim-archive sim-test sim-test-clean plot plot-test browse status remote-status restart remote-restart help

# Variables
VENV_NAME = venv
PYTHON = python3
VENV_BIN = $(VENV_NAME)/bin
VENV_PIP = $(VENV_BIN)/pip
VENV_PYTHON = $(VENV_BIN)/python

# Check if caffeinate exists, use it if available, otherwise use empty command
CAFFEINATE := $(shell command -v caffeinate 2>/dev/null)
ifdef CAFFEINATE
	PREVENT_SLEEP = caffeinate -i
else
	PREVENT_SLEEP = 
endif

help:
	@echo "Environment:"
	@echo "  make setup       - Create virtual environment and install dependencies"
	@echo "  make activate    - Print command to activate the virtual environment"
	@echo "  make setup-clean - Remove virtual environment"
	@echo ""
	@echo "Simulation:"
	@echo "  make sim            - Run simulation (archives existing data/ first)"
	@echo "                        With args: make sim ARGS=\"-n 500\""
	@echo "                        Disable progress bar: make sim ARGS=\"--no-pbar\""
	@echo "  make sim-i          - Run simulation in interactive mode"
	@echo "  make sim-repeat     - Repeat last simulation with same parameters"
	@echo "  make sim-test       - Run test simulation (cleans test_data/ first)"
	@echo "  make sim-test-clean - Remove test simulation data"
	@echo ""
	@echo "Visualization:"
	@echo "  make browse    - Open web browser to view all runs in data/"
	@echo "  make plot      - Generate plots from last simulation in data/"
	@echo "  make plot-test - Generate plots from test simulation in test_data/"
	@echo ""
	@echo "Remote (Lightsail):"
	@echo "  make status         - Check status of nanosim service when run locally"
	@echo "  make remote-status  - Check status of nanosim service on Lightsail"
	@echo "  make restart        - Restart nanosim service when run locally"
	@echo "  make remote-restart - Restart nanosim service on Lightsail"

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
	@$(PREVENT_SLEEP) $(VENV_PYTHON) creutz-sim/sim.py $(ARGS)

sim-i: sim-archive
	@echo "Running simulation in interactive mode..."
	@$(PREVENT_SLEEP) $(VENV_PYTHON) creutz-sim/sim.py -i $(ARGS)

sim-repeat:
	@SIM_INFO=$$($(VENV_PYTHON) -c "from pathlib import Path; import re, sys; files=list(Path('data').resolve().rglob('sim_started.txt')); cands=[]; \
	[(cands.append((ts.group(0), f)) if (ts:=next((re.fullmatch(r'\\d{8}_\\d{6}', part) for part in reversed(f.parts) if re.fullmatch(r'\\d{8}_\\d{6}', part)), None)) else None) for f in files]; \
	sys.exit(0) if not cands else None; latest=max(cands, key=lambda t: t[0])[1]; text=latest.read_text(); m=re.search(r'Parameters:\\s*(.*)', text); \
sys.exit(2) if m is None else None; params=dict(part.strip().split('=', 1) for part in m.group(1).split(',') if '=' in part); \
sys.exit(3) if any(k not in params for k in ('n','sweeps','flag','radius','runs')) else None; \
print(str(latest)+'|-n '+params['n']+' -s '+params['sweeps']+' -f '+params['flag']+' -r '+params['radius']+' -m '+params['runs'])" 2>/dev/null); \
	if [ -z "$$SIM_INFO" ]; then \
		echo "No previous simulation found. Running in interactive mode..."; \
		$(MAKE) sim-i; \
	else \
		LATEST_SIM=$${SIM_INFO%%|*}; \
		SIM_ARGS=$${SIM_INFO#*|}; \
		echo "Found previous simulation: $$LATEST_SIM"; \
		echo "Repeating simulation with: $$SIM_ARGS"; \
		$(PREVENT_SLEEP) $(VENV_PYTHON) creutz-sim/sim.py $$SIM_ARGS $(ARGS); \
	fi

sim-test: sim-test-clean
	@echo "Running test simulation..."
	@$(PREVENT_SLEEP) $(VENV_PYTHON) creutz-sim/sim.py -i --data-dir test_data
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

# ============================================================================
# Remote (Lightsail)
# ============================================================================

status:
	@echo "Checking nanosim service status..."
	@sudo supervisorctl status nanosim

remote-status:
	@echo "Checking nanosim service status on Lightsail..."
	@ssh plots.myember.org "sudo supervisorctl status nanosim"

restart:
	@echo "Restarting nanosim service..."
	@sudo supervisorctl restart nanosim
	@echo "Waiting for service to start..."
	@sleep 2
	@sudo supervisorctl status nanosim

remote-restart:
	@echo "Restarting nanosim service on Lightsail..."
	@ssh plots.myember.org "sudo supervisorctl restart nanosim"
	@echo "Waiting for service to start..."
	@sleep 2
	@ssh plots.myember.org "sudo supervisorctl status nanosim"