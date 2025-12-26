"""
Shared utilities for Monte Carlo simulations.

This module contains common functionality shared between Inferno and irrInferno
classes, including entropy calculations, validation methods, and initialization
helpers.
"""

import numpy as np
import argparse
import logging
import os
from scipy.special import loggamma as logg


# =============================================================================
# ENTROPY CALCULATIONS
# =============================================================================

# Legacy lambda functions (kept for backward compatibility)
Sk = lambda N, K: logg(K + N) - logg(K+1) - logg(N)
Su = lambda N, N0, Nx: logg(N+1) + np.log(2**(N0)) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))
Su0 = lambda N, N0, Nx: logg(N+1) + np.log(2**(N0+1)) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))


def Sk_stable(N: int, K: int) -> float:
    """Calculate kinetic entropy using stable logarithm.
    
    Args:
        N: Number of sites in lattice
        K: Total demon energy (sum of all oscillator energies)
        
    Returns:
        Kinetic entropy in units of k_B
        
    Note:
        Uses loggamma for numerical stability with large factorials.
        Calculates log(Ω_k) where Ω_k = (K+N-1)! / (K! * (N-1)!)
    """
    return logg(K + N) - logg(K + 1) - logg(N)


def Su_stable(N: int, N0: int, Nx: int, N0_exp: int) -> float:
    """Calculate configurational entropy using stable logarithm.
    
    Args:
        N: Number of sites in lattice
        N0: Number of broken bonds
        Nx: Number of anti-aligned neighbor pairs
        N0_exp: N0 value for exponential term (max(N0, 1) to avoid log(0))
        
    Returns:
        Configurational entropy in units of k_B
        
    Note:
        Uses N0_exp * log(2) instead of log(2^N0_exp) to avoid overflow.
        Calculates log(Ω_u) where Ω_u = N! * 2^N0 / ((N-N0-Nx)! * N0! * Nx!)
    """
    log_2_term = N0_exp * np.log(2)
    return logg(N + 1) + log_2_term - (logg(N - N0 - Nx + 1) + logg(N0 + 1) + logg(Nx + 1))


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def create_argument_parser(sim_type: str) -> argparse.ArgumentParser:
    """Create command-line argument parser for parallel simulations.
    
    Args:
        sim_type: Either 'reversible' or 'irreversible'
        
    Returns:
        Configured ArgumentParser instance
    """
    parser = argparse.ArgumentParser(
        description=f'Run parallel {sim_type} Creutz demon simulation'
    )
    parser.add_argument('--n', type=int, default=1000000, 
                       help='Lattice size (default: 1000000)')
    parser.add_argument('--s', type=int, default=5000, 
                       help='Number of sweeps per phase (default: 5000)')
    parser.add_argument('--r', type=int, default=11, 
                       help='Max demon-coupling radius, tests 1 to r-1 (default: 11)')
    parser.add_argument('--m', type=int, default=5, 
                       help='Number of independent runs (default: 5)')
    parser.add_argument('--cores', type=int, default=None, 
                       help='Number of CPU cores to use (default: auto-detect)')
    parser.add_argument('--validate', type=str, default='off', 
                       choices=['off', 'periodic', 'frequent'],
                       help='Validation mode: off (fastest), periodic (every 100 sweeps), frequent (every sweep)')
    parser.add_argument('--jit', action='store_true',
                       help='Use JIT-compiled version for 70x speedup (requires numba)')
    return parser


# =============================================================================
# SIMULATION SETUP & CONFIGURATION
# =============================================================================

def print_simulation_info(sim_type: str, n: int, s: int, r: int, m: int, 
                         num_cores: int, validate_mode: str, use_jit: bool):
    """Print simulation configuration information.
    
    Args:
        sim_type: Either 'reversible' or 'irreversible'
        n: Lattice size
        s: Number of sweeps
        r: Max radius value
        m: Number of runs per radius
        num_cores: Number of CPU cores being used
        validate_mode: Validation mode string
        use_jit: Whether JIT compilation is enabled
    """
    total_sims = r * m
    print(f"Parallel {sim_type} simulation:")
    print(f"  Lattice size: n={n}")
    print(f"  Sweeps: s={s}")
    print(f"  Radii: R=0 to R={r-1}")
    print(f"  Runs per radius: m={m}")
    print(f"  Total simulations: {total_sims}")
    print(f"  CPU cores: {num_cores}")
    print(f"  Validation mode: {validate_mode}")
    print(f"  JIT compilation: {'ENABLED' if use_jit else 'disabled'}")
    print(f"  Parallelization: {num_cores} cores")


def setup_logging(project_root: str, sim_type: str) -> str:
    """Set up logging directory and return log file path.
    
    Args:
        project_root: Project root directory path
        sim_type: Either 'reversible' or 'irreversible'
        
    Returns:
        Path to log file
    """
    log_dir = os.path.join(project_root, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    return log_dir


def create_data_directories(project_root: str, r: int, sim_type: str):
    """Create data output directories for simulations.
    
    Args:
        project_root: Project root directory path
        r: Max radius value (creates directories for R=0 to R=r-1)
        sim_type: Either 'reversible' or 'irreversible'
    """
    for R in range(r):
        if sim_type == 'irreversible':
            os.makedirs(os.path.join(project_root, 'data', 'irr', f'r{R}'), exist_ok=True)
        else:
            os.makedirs(os.path.join(project_root, 'data', f'r{R}'), exist_ok=True)


def build_simulation_parameters(r: int, m: int, n: int, s: int, validate_mode: str, 
                                project_root: str, use_jit: bool):
    """Build list of all simulation parameter tuples.
    
    Args:
        r: Max radius value
        m: Number of runs per radius
        n: Lattice size
        s: Number of sweeps
        validate_mode: Validation mode string
        project_root: Project root directory
        use_jit: Whether to use JIT compilation
        
    Returns:
        List of parameter tuples for parallel execution
    """
    sim_params = []
    sim_num = 1
    for M in range(m):
        for R in range(r):
            sim_params.append((R, M, n, s, validate_mode, project_root, use_jit, sim_num, None))
            sim_num += 1
    return sim_params


# =============================================================================
# RESULTS REPORTING
# =============================================================================

def print_final_results(results, start_time, end_time, sim_type: str):
    """Print and log final simulation results.
    
    Args:
        results: List of simulation result dictionaries
        start_time: Simulation start timestamp
        end_time: Simulation end timestamp  
        sim_type: Either 'reversible' or 'irreversible'
    """
    import logging
    from datetime import datetime
    
    elapsed = (end_time - start_time).total_seconds()
    
    # Log results
    logging.info(f"Parallel simulation complete!")
    logging.info(f"Total time: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    logging.info(f"Time per simulation: {elapsed/len(results):.1f} seconds")
    
    print(f"\nCompleted {len(results)} simulations in {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"Average time per simulation: {elapsed/len(results):.1f}s")
    
    # Verify energy conservation
    energy_errors = []
    for result in results:
        if abs(result['E_total'] - result['E_initial']) > 1e-10:
            energy_errors.append(result)
            logging.warning(f"Energy conservation error: R={result['R']}, M={result['M']}")
    
    if energy_errors:
        print(f"\nWARNING: {len(energy_errors)} simulations had energy conservation errors")
    else:
        print(f"\n✓ All simulations passed energy conservation check")


# =============================================================================
# PROGRESS DISPLAY & TIME FORMATTING
# =============================================================================

def format_time(minutes):
    """
    Format time duration in minutes to human-readable string.
    
    Args:
        minutes: Time duration in minutes
        
    Returns:
        str: Formatted time string (e.g., "45 min", "1h 52min", "2d 5h 30min")
    """
    if minutes >= 1440:  # 24 hours
        days = int(minutes // 1440)
        hours = int((minutes % 1440) // 60)
        mins = int(minutes % 60)
        return f"{days}d {hours}h {mins} min"
    elif minutes >= 60:
        hours = int(minutes // 60)
        mins = int(minutes % 60)
        return f"{hours}h {mins} min"
    else:
        return f"{int(minutes)} min"


def calculate_progress_description(active_sims, sim_times, s, total_sims, num_cores, completed):
    """
    Calculate progress bar description with ETA based on current simulation state.
    
    Args:
        active_sims: Dictionary of currently running simulations
        sim_times: List of completion times for finished simulations
        s: Number of sweeps per simulation
        total_sims: Total number of simulations to run
        num_cores: Number of CPU cores being used
        completed: Number of completed simulations
        
    Returns:
        str: Formatted description string for progress bar
    """
    import time
    import numpy as np
    
    active_count = len(active_sims)
    if active_count == 0:
        return "Overall Progress"
    
    # Show first active simulation as example with percentage
    first_sim = list(active_sims.values())[0]
    pct = (first_sim['progress'] / s) * 100
    
    # Determine what to display based on available data
    if sim_times:
        # Use actual completion times once available (most accurate)
        avg_time = np.mean(sim_times)
        remaining = total_sims - completed
        eta_minutes = (remaining * avg_time / num_cores) / 60
        return f"Overall Progress (Running {active_count} | Ex: R={first_sim['R']}, M={first_sim['M']}, {first_sim['phase']} {first_sim['progress']}/{s} {pct:.0f}% | ETA: ~{format_time(eta_minutes)} remaining)"
    elif first_sim['progress'] >= s * 0.01 and 'phase_start_time' in first_sim:
        # Calculate ETA if we have 1%+ progress
        phase_elapsed = time.time() - first_sim['phase_start_time']
        phase_time_estimate = (phase_elapsed / first_sim['progress']) * s
        sim_time_estimate = phase_time_estimate * 2  # forward + reverse
        total_eta_minutes = (total_sims * sim_time_estimate / num_cores) / 60
        return f"Overall Progress (Running {active_count} | Ex: R={first_sim['R']}, M={first_sim['M']}, {first_sim['phase']} {first_sim['progress']}/{s} {pct:.0f}% | Est: ~{format_time(total_eta_minutes)} remaining)"
    else:
        # Show running status with example (before 1% threshold)
        return f"Overall Progress (Running {active_count} | Ex: R={first_sim['R']}, M={first_sim['M']}, {first_sim['phase']} {first_sim['progress']}/{s} {pct:.0f}%)"


def handle_progress_message(msg, active_sims, s, pbar, last_refresh, sim_times, total_sims, num_cores, completed):
    """
    Handle progress update message and update display if throttle allows.
    
    Args:
        msg: Progress message dictionary
        active_sims: Dictionary of currently running simulations
        s: Number of sweeps per simulation
        pbar: tqdm progress bar object
        last_refresh: Timestamp of last display refresh
        sim_times: List of completion times
        total_sims: Total number of simulations
        num_cores: Number of CPU cores
        completed: Number of completed simulations
        
    Returns:
        float: Updated last_refresh timestamp
    """
    import time
    
    if msg['sim_num'] not in active_sims:
        return last_refresh
    
    # Track phase transitions
    old_phase = active_sims[msg['sim_num']]['phase']
    new_phase = msg['phase']
    if old_phase != new_phase:
        active_sims[msg['sim_num']]['phase_start_time'] = time.time()
    
    active_sims[msg['sim_num']]['phase'] = new_phase
    active_sims[msg['sim_num']]['progress'] = msg['sweep']
    
    # Only update display if throttle allows (1 second intervals)
    now = time.time()
    if now - last_refresh >= 1.0:
        desc = calculate_progress_description(active_sims, sim_times, s, total_sims, num_cores, completed)
        pbar.set_description(desc)
        return now
    
    return last_refresh


def process_message_queue(progress_queue, active_sims, s, pbar, last_refresh, 
                          sim_times, total_sims, num_cores, completed, result, pool):
    """
    Process messages from the progress queue and update tracking state.
    
    Args:
        progress_queue: Queue containing progress messages
        active_sims: Dictionary tracking active simulations
        s: Number of sweeps per simulation
        pbar: tqdm progress bar
        last_refresh: Last display refresh timestamp
        sim_times: List of completion times
        total_sims: Total number of simulations
        num_cores: Number of CPU cores
        completed: Current completion count
        result: AsyncResult from pool.map_async
        pool: Multiprocessing pool instance
        
    Returns:
        tuple: (completed count, last_refresh timestamp, results list or None)
    """
    import time
    import sys
    
    try:
        while not result.ready() or not progress_queue.empty():
            try:
                msg = progress_queue.get(timeout=0.1)
                    
                if msg['type'] == 'start':
                    active_sims[msg['sim_num']] = {
                        'R': msg['R'],
                        'M': msg['M'],
                        'start_time': time.time(),
                        'phase': 'forward',
                        'progress': 0,
                        'phase_start_time': time.time()
                    }
                    
                elif msg['type'] == 'progress':
                    last_refresh = handle_progress_message(
                        msg, active_sims, s, pbar, last_refresh, 
                        sim_times, total_sims, num_cores, completed
                    )
                
                elif msg['type'] == 'complete':
                    if msg['sim_num'] in active_sims:
                        elapsed = time.time() - active_sims[msg['sim_num']]['start_time']
                        sim_times.append(elapsed)
                        del active_sims[msg['sim_num']]
                        
                        completed += 1
                        pbar.update(1)
                            
            except Exception:
                pass
            
            # Force refresh display every second to show current status
            if time.time() - last_refresh > 1.0:
                pbar.refresh()
                last_refresh = time.time()
    
    except KeyboardInterrupt:
        pbar.close()
        pool.terminate()
        pool.join()
        print("\n\nSimulation interrupted by user (Ctrl-C)")
        print(f"Completed {completed}/{total_sims} simulations before interruption")
        sys.exit(130)
    
    pbar.close()
    return completed, last_refresh, result.get()


# =============================================================================
# BASE SIMULATION CLASS
# =============================================================================

class SimulationBase:
    """
    Base class for Monte Carlo simulations with common validation and
    initialization methods.
    
    This class provides shared functionality for energy conservation,
    bond count validation, and state retrieval that is common to both
    reversible (Inferno) and irreversible (irrInferno) simulations.
    """
    
    @staticmethod
    def initialize_order_arrays(N):
        """
        Initialize randomized order arrays for lattice traversal.
        
        Creates a shuffled array for forward traversal and its reverse
        for backward traversal (used in reversibility tests).
        
        Args:
            N: Number of lattice sites
            
        Returns:
            tuple: (order, rev_order) - forward and reverse traversal arrays
        """
        order = np.arange(N)
        np.random.shuffle(order)
        rev_order = np.flip(order)
        
        return order, rev_order
    
    @staticmethod
    def initialize_lattice(N):
        """
        Initialize a lattice with half spins up (+1) and half down (-1).
        
        Args:
            N: Number of lattice sites
            
        Returns:
            tuple: (lattice, bonds, bond_count) arrays
        """
        lattice = np.concatenate((np.ones(N//2, dtype=np.int8),
                                 (-1)*np.ones(N//2, dtype=np.int8)))
        bonds = np.ones(N, dtype=np.int8)*(-1)
        bonds[[N//2-1, -1]] = 1
        bond_count = np.array([N-2, 0, 2], dtype=np.int64)
        
        return lattice, bonds, bond_count
    
    @staticmethod
    def initialize_demon_energy(N, total_energy, E_lattice):
        """
        Initialize demon energy array by distributing remaining energy randomly.
        
        Args:
            N: Number of lattice sites (demons)
            total_energy: Total energy to conserve
            E_lattice: Current lattice energy
            
        Returns:
            tuple: (E_demon array, E_demon_sum, d_energy)
        """
        d_energy = total_energy - E_lattice
        result = np.zeros(N, dtype=np.int64)
        indices = np.random.randint(0, N, size=d_energy)
        np.add.at(result, indices, 1)
        
        return result, np.int64(d_energy), d_energy
    
    @staticmethod
    def setup_neighbor_arrays(N):
        """
        Pre-compute neighbor index arrays for periodic boundaries.
        
        This optimization replaces modulo operations with direct array lookups.
        
        Args:
            N: Number of lattice sites
            
        Returns:
            tuple: (right_neighbor, left_neighbor) arrays
        """
        right_neighbor = np.arange(1, N+1, dtype=np.int32) % N
        left_neighbor = np.arange(-1, N-1, dtype=np.int32) % N
        
        return right_neighbor, left_neighbor
    
    def setup_validation_mode(self, N, validate_mode, initial_total_energy):
        """
        Configure validation parameters based on mode.
        
        Args:
            N: Number of lattice sites (used to calculate check intervals)
            validate_mode: 'off', 'periodic', or 'frequent'
            initial_total_energy: Total energy for conservation checks
        """
        self._initial_total_energy = np.int64(initial_total_energy)
        self._validate_mode = validate_mode
        self._check_counter = 0
        
        # Configure validation interval based on mode
        if validate_mode == 'off':
            self._check_interval = float('inf')  # Never validate automatically
        elif validate_mode == 'frequent':
            self._check_interval = N  # Every sweep (debug mode)
        else:  # 'periodic' or default
            self._check_interval = 100 * N  # Every 100 sweeps
    
    def perform_periodic_validation(self):
        """
        Perform periodic validation check if enabled.
        
        Should be called at the end of demon_move() and demon_reverse().
        Increments counter and validates when interval is reached.
        """
        if self._validate_mode != 'off':
            self._check_counter += 1
            if self._check_counter >= self._check_interval:
                self._check_counter = 0
                self.validate_energy_conservation()
                self.validate_bond_counts()
    
    def validate_energy_conservation(self):
        """
        Periodic energy conservation check to catch drift.
        
        Recalculates total energy from scratch and compares with tracked value.
        If drift is detected, corrects the cached values and prints warning.
        
        Returns:
            bool: True if energy is conserved, False if drift was detected and corrected
        """
        current_total = self.E_lattice + self.E_demon_sum
        if current_total != self._initial_total_energy:
            # Recalculate from scratch
            actual_lattice = np.sum(self.bonds, dtype=np.int64)
            actual_demon = np.sum(self.E_demon, dtype=np.int64)

            print(f"WARNING: Energy drift detected!")
            print(f"  Expected total: {self._initial_total_energy}")
            print(f"  Tracked total: {current_total}")
            print(f"  Actual total: {actual_lattice + actual_demon}")
            print(f"  Drift: {current_total - self._initial_total_energy}")

            # Correct the cached values
            self.E_lattice = actual_lattice
            self.E_demon_sum = actual_demon
            self.d_energy = self.E_demon_sum

            return False
        return True
    
    def validate_bond_counts(self):
        """
        Periodic bond count validation to catch drift.
        
        Recalculates bond counts from scratch and compares with tracked values.
        If drift is detected, corrects the cached values and prints warning.
        
        Returns:
            bool: True if bond counts match, False if drift was detected and corrected
        """
        actual_counts = np.bincount(self.bonds + 1, minlength=3).astype(np.int64)
        if not np.array_equal(actual_counts, self.bond_count):
            print(f"WARNING: Bond count drift detected!")
            print(f"  Tracked: {self.bond_count}")
            print(f"  Actual: {actual_counts}")

            # Correct the cached values
            self.bond_count = actual_counts
            return False
        return True
    
    def get_validated_state(self):
        """
        Return current state with validation.
        
        Always recalculates values from scratch rather than using cached values.
        Use this method when you need guaranteed accurate values.
        
        Returns:
            dict: Dictionary with keys:
                - E_lattice: Lattice energy
                - E_demon_sum: Total demon energy
                - E_total: Total energy (lattice + demon)
                - bond_count: Array of bond counts [broken, neutral, anti-aligned]
                - d_energy: Demon energy (same as E_demon_sum)
        """
        # Force validation - recalculate from scratch
        actual_lattice = np.sum(self.bonds, dtype=np.int64)
        actual_demon = np.sum(self.E_demon, dtype=np.int64)
        actual_bond_counts = np.bincount(self.bonds + 1, minlength=3).astype(np.int64)

        # Return validated values
        return {
            'E_lattice': actual_lattice,
            'E_demon_sum': actual_demon,
            'E_total': actual_lattice + actual_demon,
            'bond_count': actual_bond_counts,
            'd_energy': actual_demon
        }
    
    def spin_flip(self, a, i):
        """
        Attempt to flip the spin at site a using demon i.
        
        Uses all integer arithmetic to prevent floating point drift.
        
        Args:
            a: Lattice site index
            i: Demon index
            
        Returns:
            bool: True if spin was flipped (bonds may have changed), False otherwise
        """
        s = self.lattice[a]
        d = self.E_demon[i]

        # Calculate energy cost (integer arithmetic)
        nb = (self.lattice[self.right_neighbor[a]] * abs(self.bonds[a]) +
              self.lattice[self.left_neighbor[a]] * abs(self.bonds[self.left_neighbor[a]]))

        cost = 2 * s * nb  # Always an integer
        bonds_changed = False

        if cost < 0 or cost <= d:
            s *= -1
            # Update energies using integer arithmetic only
            self.E_demon[i] -= cost
            self.E_demon_sum -= cost
            self.d_energy -= cost
            self.E_lattice += cost
            self.lattice[a] = s
            bonds_changed = True

        return bonds_changed
    
    def update_bonds_incremental(self, a):
        """
        Update bond states after spin flip at site a.
        
        Updates both the right bond (at site a) and left bond (at left_neighbor[a])
        using careful integer counting to maintain bond_count accuracy.
        
        Args:
            a: Lattice site index where spin was flipped
        """
        # Update right bond
        if self.bonds[a] != 0:
            old_bond = self.bonds[a]
            new_bond = np.int8(-1 if self.lattice[a] == self.lattice[self.right_neighbor[a]] else 1)

            if old_bond != new_bond:
                self.bonds[a] = new_bond
                # Update counts with explicit indexing
                old_idx = 0 if old_bond == -1 else 2
                new_idx = 0 if new_bond == -1 else 2
                self.bond_count[old_idx] -= 1
                self.bond_count[new_idx] += 1

        # Update left bond
        left_idx = self.left_neighbor[a]
        if self.bonds[left_idx] != 0:
            old_bond = self.bonds[left_idx]
            new_bond = np.int8(-1 if self.lattice[a] == self.lattice[left_idx] else 1)

            if old_bond != new_bond:
                self.bonds[left_idx] = new_bond
                old_idx = 0 if old_bond == -1 else 2
                new_idx = 0 if new_bond == -1 else 2
                self.bond_count[old_idx] -= 1
                self.bond_count[new_idx] += 1
    
    def bond_change(self, a, i):
        """
        Attempt to change the bond state at site a using demon i.
        
        Can either break an existing bond (0 → ±1) or create a bond (±1 → 0).
        Uses integer arithmetic only to prevent drift.
        
        Args:
            a: Lattice site index
            i: Demon index
        """
        s = self.lattice[a]
        b = self.bonds[a]
        d = self.E_demon[i]
        n = self.lattice[self.right_neighbor[a]]

        cost = -1 if s == n else 1  # Always ±1

        if b == 0 and d - cost >= 0:
            # Update energies
            self.E_lattice += cost
            self.E_demon[i] -= cost
            self.E_demon_sum -= cost
            self.d_energy -= cost
            self.bonds[a] = np.int8(cost)

            # Update bond_count: broken → aligned/misaligned
            self.bond_count[1] -= 1
            self.bond_count[0 if cost == -1 else 2] += 1

        elif b != 0 and d + cost >= 0:
            # Update energies
            self.E_lattice -= cost
            self.E_demon[i] += cost
            self.E_demon_sum += cost
            self.d_energy += cost

            # Update bond_count
            old_idx = 0 if self.bonds[a] == -1 else 2
            self.bond_count[old_idx] -= 1
            self.bond_count[1] += 1

            self.bonds[a] = 0

        # Update left neighbor bond
        left_idx = self.left_neighbor[a]
        if self.bonds[left_idx] != 0:
            old_bond = self.bonds[left_idx]
            new_bond = np.int8(-1 if self.lattice[a] == self.lattice[left_idx] else 1)

            if old_bond != new_bond:
                self.bonds[left_idx] = new_bond
                old_idx = 0 if old_bond == -1 else 2
                new_idx = 0 if new_bond == -1 else 2
                self.bond_count[old_idx] -= 1
                self.bond_count[new_idx] += 1
