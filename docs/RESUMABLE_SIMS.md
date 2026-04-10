# Design Spec: Resumable Simulations

**Status**: Proposed  
**Context**: Long runs (30+ days) on a MacBook M3 Max where interruption is highly likely.

---

## Table of Contents

1. [Physics Safety Guarantee](#physics-safety-guarantee)
1. [Problem Statement](#problem-statement)
2. [How the Sim Currently Works](#how-the-sim-currently-works)
3. [The Natural Resume Unit](#the-natural-resume-unit)
4. [Proposed Approach: Skip-Completed-Unit Resume](#proposed-approach-skip-completed-unit-resume)
5. [Detailed Design](#detailed-design)
   - [Completeness Check](#completeness-check)
   - [sim.py Changes](#simpy-changes)
   - [Makefile Changes](#makefile-changes)
   - [Status File Updates](#status-file-updates)
6. [What Gets Lost on a Crash](#what-gets-lost-on-a-crash)
7. [macOS-Specific Hardening](#macos-specific-hardening)
8. [Alternatives Considered](#alternatives-considered)

---

## Physics Safety Guarantee

Nothing in `inferno.py` is touched. The JIT functions (`spin_flip_jit`, `bond_change_jit`),
the `Inferno` class, `demon_move`, `demon_reverse`, `reset`, and the entropy formula are all
completely unchanged. The spec does not require any modifications to them.

The changes are entirely outside the physics loop:

| Change | Where | Touches physics? |
|---|---|---|
| `--resume` CLI arg | `get_params()` startup | No |
| Completeness check | Before pool starts | No |
| Refactor worker into `run_single_unit` | Same loop body, new function name | No — identical code, moved |
| `tqdm initial=already_done_sweeps` | Progress display only | No |
| Append "Resumed:" to `sim_started.txt` | Status file I/O | No |
| `caffeinate -d -i -s` | OS process flag | No |

**The one thing worth scrutinizing:** when a partial CSV is deleted and re-run, the unit
starts with a fresh `Inferno(n, R)` — same as it would in any fresh start. The random
shuffle and energy distribution are regenerated. That is correct: the runs are *designed*
to be statistically independent, and a re-run unit is simply a new independent sample for
that `(R, M)` slot. The data is valid.

The promise holds because the implementation strategy was specifically chosen to require
zero changes to the simulation code. The resume logic sits entirely in the orchestration
layer above it.

---

## Problem Statement

A full run with large parameters (e.g. `n=1000, s=100000, r=10, m=5`) on a MacBook M3 Max
is expected to take 30+ days. Over that window, interruption is near-certain:

- System sleep (lid close, idle timeout)
- macOS automatic updates requiring a restart
- Power loss
- Accidental terminal close
- Kernel panic

Currently, any interruption loses all in-progress work. There is no way to resume; the user
must start from scratch.

---

## How the Sim Currently Works

The simulation is parallelized over **radius values** (`R = 0` to `r`, inclusive).
A `multiprocessing.Pool` assigns each `R` value to a CPU core.

Each pool worker (`run_radius_simulations`) does the following for its assigned `R`:

```
for M in range(m):                          # m independent runs per radius
    for dynamics_flag in [0, 1]:            # rev=0, irr=1 (if flag=='c')
        write CSV header to new file
        for i in range(s//2):               # forward sweeps
            run n demon_move calls
            write one row to CSV
        for i in range(s//2):               # reverse sweeps
            run n demon_reverse calls
            write one row to CSV
        x.reset()                           # new shuffle order for next run
```

Each `(R, M, dynamics_flag)` combination writes to a single independent CSV file:

| dynamics_flag | file path |
|---|---|
| 0 (rev) | `data/<timestamp>/rev/r{R}/sim_data_r{R}_{M}.csv` |
| 1 (irr) | `data/<timestamp>/irr/r{R}/irr_sim_data_r{R}_{M}.csv` |

A complete CSV has exactly `s + 1` lines (1 header + `s` data rows: `s//2` forward + `s//2` reverse).

---

## The Natural Resume Unit

The `(R, M, dynamics_flag)` triple maps 1:1 to a single CSV file. This is the perfect
resume unit because:

1. **Each file is independent.** Workers share no state except the `pbar_queue`.
2. **`Inferno` is stateless between units.** `x.reset()` regenerates a fresh random order;
   there is nothing to carry forward.
3. **Completeness is trivially verifiable.** Count lines in the CSV. A complete file has
   exactly `s + 1` lines.
4. **The parameters needed to rerun any unit are already written to `sim_started.txt`.**

A crash means some in-flight CSVs are partial (fewer than `s + 1` lines). On resume, those
partials are removed and their units are re-queued.

---

## Proposed Approach: Skip-Completed-Unit Resume

On restart, `sim.py` is told to target an existing timestamped directory. It:

1. Reads `sim_started.txt` to recover original parameters (`n, s, flag, r, m`).
2. Scans all expected CSV files. Classifies each as **complete**, **partial**, or **missing**.
3. Deletes partials (they contain data from an interrupted run that cannot be appended to
   cleanly, because `Inferno` state is not recoverable mid-run).
4. Builds the work queue from the **incomplete** `(R, M, dynamics_flag)` units only.
5. Submits that reduced work queue to the pool.
6. Reuses the same data directory and status files.

At the end, the directory looks identical to a run that was never interrupted.

---

## Detailed Design

### Completeness Check

```python
def count_csv_rows(path: Path) -> int:
    """Return number of data rows in CSV (excluding header). Returns 0 if missing."""
    if not path.exists():
        return 0
    with open(path) as f:
        return sum(1 for _ in f) - 1  # subtract header
    
def is_unit_complete(path: Path, expected_rows: int) -> bool:
    return count_csv_rows(path) >= expected_rows
```

`expected_rows` is always `s` (the `sweeps` parameter): `s//2` forward + `s//2` reverse.

---

### sim.py Changes

#### New CLI argument

```
--resume TIMESTAMP    Resume an interrupted run. TIMESTAMP is the directory name
                      under data/ (e.g. 20260409_142300). All original parameters
                      are read from sim_started.txt; any parameters passed on the
                      command line are ignored.
```

#### Updated `get_params()` return value

Add `resume_dir: str | None` to the returned tuple.

#### New function: `build_work_queue()`

```python
def build_work_queue(r, m, flag, s, file_names, irr_files):
    """
    Return list of (R, M, dynamics_flag) units that are not yet complete.
    Removes partial CSV files in-place.
    """
    pending = []
    if flag == 'c':
        dynamics_flags = [0, 1]
    elif flag == 'r':
        dynamics_flags = [0]
    else:
        dynamics_flags = [1]
    
    for R in range(r + 1):
        for M in range(m):
            for df in dynamics_flags:
                if df == 0:
                    path = file_names[R].parent / f"{file_names[R].name}_{M}.csv"
                else:
                    path = irr_files[R].parent / f"{irr_files[R].name}_{M}.csv"
                
                if is_unit_complete(path, s):
                    continue   # already done, skip
                
                if path.exists():
                    path.unlink()  # remove partial, will be rerun cleanly
                
                pending.append((R, M, df))
    
    return pending
```

#### Updated `run_radius_simulations()`

Currently the worker loops over `for M in range(m)`. With resumability, the pool must be
able to submit individual `(R, M, dynamics_flag)` units, not whole `R` batches.

The worker signature changes to:

```python
def run_single_unit(R, M, dynamics_flag, n, s, file_names, irr_files, init_files, pbar_queue):
    """Run one (R, M, dynamics_flag) unit."""
    ...
```

The inner loop body of the current `run_radius_simulations` moves here verbatim.
The pool maps over the `pending` list rather than `range(r+1)`.

This is a pure refactor of the existing loop — no logic changes to the simulation itself.

#### Updated `__main__` block

```python
# Detect resume mode
if resume_dir:
    data_folder = data_root / resume_dir
    # read original parameters from sim_started.txt
    n, s, flag, r, m = parse_started_txt(data_folder / ... / 'sim_started.txt')
    timestamp = resume_dir
    # append resume event to start marker
    with open(start_marker, 'a') as f:
        f.write(f"Resumed: {datetime.now().isoformat()}\n")
else:
    # existing fresh-start logic (create timestamped dir, write sim_started.txt)
    ...

# Build work queue (works for both fresh and resume)
pending_units = build_work_queue(r, m, flag, s, file_names, irr_files)
completed_units = total_units - len(pending_units)  # units already fully written to disk

# Total sweeps across the entire run (original scope, not just remaining)
# Used to initialise tqdm so the bar feels like the run was never interrupted.
sweeps_per_unit = s  # s//2 fwd + s//2 rev
multiplier = 2 if flag == 'c' else 1  # combined runs both dynamics per (R, M)
total_sweeps = total_units * multiplier * sweeps_per_unit
already_done_sweeps = completed_units * multiplier * sweeps_per_unit

# Progress bar starts at already_done_sweeps so it picks up where it left off
# tqdm: pbar = tqdm(total=total_sweeps, initial=already_done_sweeps, ...)

# Pool maps over pending_units instead of range(r+1)
worker_func = partial(run_single_unit, ...)
results = pool.map_async(worker_func, pending_units)
```

---

### Makefile Changes

#### New target: `sim-resume`

```makefile
sim-resume:
	@LAST=$$(ls -1t data/ | grep -E '^[0-9]{8}_[0-9]{6}$$' | head -1); \
	if [ -z "$$LAST" ]; then \
		echo "No previous run found in data/"; exit 1; \
	fi; \
	echo "Resuming $$LAST"; \
	$(PREVENT_SLEEP) $(VENV_PYTHON) creutz-sim/sim.py --resume $$LAST $(ARGS)
```

This finds the most recent timestamped directory automatically. The user can also specify
one explicitly: `make sim-resume ARGS="--resume 20260409_142300"`.

---

### Status File Updates

The existing `sim_status.txt` pattern is extended to support resume state.

**On fresh start** (unchanged):
```
Status: RUNNING
Started: 2026-04-09T14:23:00
Parameters: n=1000, sweeps=100000, flag=c, radius=10, runs=5
```

**On resume** (appended section):
```
Resumed: 2026-04-09T19:05:00
Completed units at resume: 67 / 110
Completed sweeps at resume: 6700000 / 11000000
```

Logging `completed_sweeps` here means the browser can display accurate overall progress
even while the sim is mid-resume, and the final `sim_completed.txt` can record the true
total elapsed wall-clock time across all sessions.

**On completion** (unchanged):
```
Status: COMPLETED
...
```

---

## What Gets Lost on a Crash

With the skip-completed-unit approach, only **in-flight units** at the moment of the crash
are lost. An "in-flight unit" is a `(R, M, dynamics_flag)` triple currently being executed
by a worker.

- Number of in-flight units at any moment = `num_cores` (at most)
- `num_cores = min(mp.cpu_count(), r + 1)`
- On an M3 Max (16 performance cores), with `r=10`, `num_cores = 11`
- Maximum lost work = 11 units × `s` sweeps each

For a 30-day run with `s=100,000` and `r=10, m=5` (110 total units), worst case is losing
~10% of total work no matter when the crash happens. In practice, crashes later in the run
lose only the remaining fraction. Average expected re-work: ~0.5 units per resume event.

---

## macOS-Specific Hardening

These are recommended regardless of the resume implementation.

### Current: `caffeinate -i`

Prevents idle sleep only. Does **not** prevent:
- Manual sleep (lid close, Apple menu → Sleep)
- Sleep on low battery
- Automatic restart for macOS updates

### Recommended additions

| Risk | Mitigation |
|---|---|
| Accidental terminal close | Run inside `tmux`: `tmux new -s nanosim` then `make sim` |
| Manual/automatic sleep | Add `-d` (prevent display sleep) and `-s` (prevent AC system sleep) to `caffeinate` flags in Makefile |
| macOS auto-update restart | System Settings → General → Software Update → disable "Install macOS updates" for the run duration |
| Battery depleted | Keep plugged in; `caffeinate -s` only works on AC power anyway |

**Updated Makefile `PREVENT_SLEEP`:**

```makefile
CAFFEINATE := $(shell command -v caffeinate 2>/dev/null)
ifdef CAFFEINATE
    PREVENT_SLEEP = caffeinate -d -i -s
else
    PREVENT_SLEEP =
endif
```

**tmux workflow** (no code changes required today):

```bash
tmux new -s nanosim        # start named session
make sim ARGS="-n 1000 ..."
# detach: Ctrl-B then D
# reattach later:
tmux attach -t nanosim
```

If the terminal app crashes or the SSH session drops, the sim keeps running and the progress
bar is visible when you reattach.

---

## Alternatives Considered

### Mid-sweep state checkpointing

Serialize the full `Inferno` state (lattice, bonds, E_demon arrays, counters) every N sweeps
using `numpy.save`. On resume, reload and continue from the last checkpoint.

**Pros:** Loses at most one checkpoint interval of work (could be minutes, not hours).  
**Cons:** Complex with `multiprocessing.Pool` — workers would need a safe moment to flush
state to disk, adding synchronization overhead. The partial CSV problem still exists (a
partially-written checkpoint sweep row must be truncated). For runs where each unit takes
hours, the added complexity is not worth the marginal recovery improvement.

**Verdict:** Not recommended at this time. Revisit if individual units grow to be multi-day.

### Separate process per unit (orchestrator pattern)

A shell script or Python orchestrator launches each `(R, M, dynamics_flag)` unit as a
separate `subprocess`. Completed units are tracked in a manifest file. The orchestrator
is restartable.

**Pros:** Very robust; each unit is completely isolated.  
**Cons:** Loses the `multiprocessing.Pool` parallelism and shared progress bar. More
infrastructure to build and maintain.

**Verdict:** Not recommended. The pool-based approach already provides the needed isolation
at the `(R, M)` level.
