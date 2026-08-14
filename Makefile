.PHONY: setup activate setup-clean sim sim-i sim-repeat sim-archive sim-test sim-test-clean plot plot-test plot-cache plot-cache-all plot-cache-sync-remote publish-run browse status remote-status restart remote-restart debug remote-debug help

# Variables
VENV_NAME = venv
PYTHON = python3
VENV_BIN = $(VENV_NAME)/bin
VENV_PIP = $(VENV_BIN)/pip
VENV_PYTHON = $(VENV_BIN)/python
REMOTE_HOST ?= plots.myember.org
REMOTE_DATA_DIR ?= /var/www/nanosim/data
REMOTE_SSH_TARGET ?= $(REMOTE_HOST)

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
	@echo "  make plot-cache - Build dark+light plot cache for a run locally"
	@echo "                     Example: make plot-cache ARGS=\"--run 20260725_143642\""
	@echo "  make plot-cache-all - Rebuild dark+light plot cache for every run in data/"
	@echo ""
	@echo "Remote (Lightsail):"
	@echo "  make status         - Check status of nanosim service when run locally"
	@echo "  make remote-status  - Check status of nanosim service on Lightsail"
	@echo "  make restart        - Restart nanosim service when run locally"
	@echo "  make remote-restart - Restart nanosim service on Lightsail"
	@echo "  make publish-run RUN=<run_id>"
	@echo "                      - Upload data/<run_id>/sim_* and plot_cache/ to Lightsail"
	@echo "                        and normalize permissions for nginx"
	@echo "                        Override target/path: make publish-run RUN=... REMOTE_SSH_TARGET=... REMOTE_DATA_DIR=..."
	@echo "  make plot-cache-sync-remote - Copy all data/*/plot_cache to Lightsail"
	@echo "                                Override host/path: make plot-cache-sync-remote REMOTE_HOST=... REMOTE_DATA_DIR=..."
	@echo "  make debug          - Show local nanosim status + recent logs"
	@echo "  make remote-debug   - Show remote nanosim status + recent logs"

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
	@SIM_INFO=$$($(VENV_PYTHON) tools/sim_repeat_info.py 2>/dev/null); \
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

plot-cache:
	@echo "Building plot cache..."
	@$(VENV_PYTHON) tools/build_plot_cache.py $(ARGS)

plot-cache-all:
	@echo "Rebuilding plot cache for all runs in data/..."
	@set -e; \
	runs=$$(find data -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort); \
	if [ -z "$$runs" ]; then \
		echo "No run directories found under data/."; \
		exit 0; \
	fi; \
	for run in $$runs; do \
		echo ""; \
		echo "=== $$run ==="; \
		$(VENV_PYTHON) tools/build_plot_cache.py --run "$$run" || exit $$?; \
	done; \
	echo ""; \
	echo "Finished rebuilding plot cache for all runs."

plot-cache-sync-remote:
	@echo "Syncing plot_cache directories to $(REMOTE_HOST):$(REMOTE_DATA_DIR)..."
	@set -e; \
	cache_dirs=$$(find data -mindepth 2 -maxdepth 2 -type d -name plot_cache | sort); \
	if [ -z "$$cache_dirs" ]; then \
		echo "No plot_cache directories found under data/."; \
		exit 0; \
	fi; \
	echo "Planned copies:"; \
	for cache_dir in $$cache_dirs; do \
		run_dir=$$(basename $$(dirname "$$cache_dir")); \
		echo "  $$cache_dir/ -> $(REMOTE_HOST):$(REMOTE_DATA_DIR)/$$run_dir/plot_cache/"; \
	done; \
	echo ""; \
	printf "Proceed with sync? [y/N] "; \
	read confirm; \
	case "$$confirm" in \
		y|Y|yes|YES) ;; \
		*) echo "Aborted."; exit 0;; \
	esac; \
	for cache_dir in $$cache_dirs; do \
		run_dir=$$(basename $$(dirname "$$cache_dir")); \
		echo ""; \
		echo "=== $$run_dir ==="; \
		ssh "$(REMOTE_HOST)" "mkdir -p '$(REMOTE_DATA_DIR)/$$run_dir/plot_cache'"; \
		rsync -az "$$cache_dir/" "$(REMOTE_HOST):$(REMOTE_DATA_DIR)/$$run_dir/plot_cache/"; \
	done; \
	echo ""; \
	echo "Finished syncing plot_cache directories to $(REMOTE_HOST)."

publish-run:
	@set -e; \
	if [ -z "$(RUN)" ]; then \
		echo "RUN is required. Example: make publish-run RUN=20260809_100023"; \
		exit 2; \
	fi; \
	run_id="$(patsubst %/,%,$(RUN))"; \
	run_dir="data/$$run_id"; \
	if [ ! -d "$$run_dir" ]; then \
		echo "Run directory not found: $$run_dir"; \
		exit 2; \
	fi; \
	if [ ! -d "$$run_dir/plot_cache" ]; then \
		echo "Missing plot cache directory: $$run_dir/plot_cache"; \
		exit 2; \
	fi; \
	if ! find "$$run_dir" -maxdepth 1 -type f -name 'sim_*' | grep -q .; then \
		echo "No sim_* files found under $$run_dir"; \
		exit 2; \
	fi; \
	flag=""; \
	if [ -f "$$run_dir/sim_started.txt" ]; then \
		flag=$$(grep -oE 'flag=[a-z]' "$$run_dir/sim_started.txt" | head -n 1 | cut -d= -f2 || true); \
	fi; \
	if [ -z "$$flag" ] && [ -f "$$run_dir/sim_status.txt" ]; then \
		flag=$$(grep -oE 'flag=[a-z]' "$$run_dir/sim_status.txt" | head -n 1 | cut -d= -f2 || true); \
	fi; \
	if [ -z "$$flag" ]; then \
		flag="c"; \
	fi; \
	case "$$flag" in \
		c) mkdir -p "$$run_dir/rev" "$$run_dir/irr" ;; \
		r) mkdir -p "$$run_dir/rev" ; mkdir -p "$$run_dir/irr" ;; \
		i) mkdir -p "$$run_dir/irr" ; mkdir -p "$$run_dir/rev" ;; \
		*) mkdir -p "$$run_dir/rev" "$$run_dir/irr" ;; \
	esac; \
	remote_run_dir="$(REMOTE_DATA_DIR)/$$run_id"; \
	echo "Publishing $$run_dir to $(REMOTE_SSH_TARGET):$$remote_run_dir"; \
	echo "  - sim_* files"; \
	echo "  - plot_cache/"; \
	echo "  - dynamics placeholders: flag=$$flag (rev/irr)"; \
	echo ""; \
	printf "Proceed with publish? [y/N] "; \
	read confirm; \
	case "$$confirm" in \
		y|Y|yes|YES) ;; \
		*) echo "Aborted."; exit 0;; \
	esac; \
	echo "Creating remote directories..."; \
	ssh "$(REMOTE_SSH_TARGET)" "mkdir -p '$$remote_run_dir/plot_cache' && case '$$flag' in c) mkdir -p '$$remote_run_dir/rev' '$$remote_run_dir/irr' ;; r) mkdir -p '$$remote_run_dir/rev' ; mkdir -p '$$remote_run_dir/irr' ;; i) mkdir -p '$$remote_run_dir/irr' ; mkdir -p '$$remote_run_dir/rev' ;; *) mkdir -p '$$remote_run_dir/rev' '$$remote_run_dir/irr' ;; esac"; \
	echo "Syncing files via rsync..."; \
	rsync -az -v --progress --include='sim_*' --include='plot_cache/***' --exclude='*' "$$run_dir/" "$(REMOTE_SSH_TARGET):$$remote_run_dir/"; \
	echo "Normalizing permissions..."; \
	ssh "$(REMOTE_SSH_TARGET)" "find '$$remote_run_dir' -type d -exec chmod 755 {} + && find '$$remote_run_dir' -type f -exec chmod 644 {} +"; \
	echo ""; \
	echo "Publish complete: $(REMOTE_SSH_TARGET):$$remote_run_dir"

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
	@ssh "$(REMOTE_HOST)" "sudo supervisorctl status nanosim"

debug:
	@echo "Local nanosim status:"
	@sudo supervisorctl status nanosim || true
	@echo ""
	@echo "Recent local error log (/var/log/nanosim/error.log):"
	@sudo tail -n 80 /var/log/nanosim/error.log || true
	@echo ""
	@echo "Recent local access log (/var/log/nanosim/access.log):"
	@sudo tail -n 80 /var/log/nanosim/access.log || true

remote-debug:
	@echo "Remote nanosim status on Lightsail:"
	@ssh "$(REMOTE_HOST)" "sudo supervisorctl status nanosim" || true
	@echo ""
	@echo "Recent remote error log (/var/log/nanosim/error.log):"
	@ssh "$(REMOTE_HOST)" "sudo tail -n 80 /var/log/nanosim/error.log" || true
	@echo ""
	@echo "Recent remote access log (/var/log/nanosim/access.log):"
	@ssh "$(REMOTE_HOST)" "sudo tail -n 80 /var/log/nanosim/access.log" || true

restart:
	@echo "Restarting nanosim service..."
	@sudo supervisorctl restart nanosim
	@echo "Waiting for service to start..."
	@attempt=1; \
	max_attempts=12; \
	while [ $$attempt -le $$max_attempts ]; do \
		status_line=$$(sudo supervisorctl status nanosim); \
		echo "[$$attempt/$$max_attempts] $$status_line"; \
		echo "$$status_line" | grep -q "RUNNING" && exit 0; \
		sleep 2; \
		attempt=$$((attempt + 1)); \
	done; \
	echo "nanosim failed to reach RUNNING state after $$max_attempts attempts"; \
	echo "Recent service error log:"; \
	sudo tail -n 80 /var/log/nanosim/error.log || true; \
	exit 7

remote-restart:
	@echo "Restarting nanosim service on Lightsail..."
	@ssh "$(REMOTE_HOST)" "sudo supervisorctl restart nanosim"
	@echo "Waiting for service to start..."
	@attempt=1; \
	max_attempts=12; \
	while [ $$attempt -le $$max_attempts ]; do \
		status_line=$$(ssh "$(REMOTE_HOST)" "sudo supervisorctl status nanosim"); \
		echo "[$$attempt/$$max_attempts] $$status_line"; \
		echo "$$status_line" | grep -q "RUNNING" && exit 0; \
		sleep 2; \
		attempt=$$((attempt + 1)); \
	done; \
	echo "nanosim failed to reach RUNNING state on Lightsail after $$max_attempts attempts"; \
	echo "Recent remote error log:"; \
	ssh "$(REMOTE_HOST)" "sudo tail -n 80 /var/log/nanosim/error.log" || true; \
	exit 7