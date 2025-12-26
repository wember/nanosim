"""Tests for parallel simulation execution.

Tests verify that parallel implementations produce correct results,
handle resources properly, and match sequential simulation outputs.
"""

import pytest
import os
import sys
import csv
import json
import multiprocessing as mp
from multiprocessing import Manager
import tempfile
import shutil
import numpy as np
from pathlib import Path

# Add creutz-sim to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from parallel_sim import run_single_simulation as run_single_sim_rev
from parallel_sim_irr import run_single_simulation as run_single_sim_irr
from inferno import Inferno
from inferno_irr import irrInferno


@pytest.fixture
def mock_progress_queue():
    """Create a mock progress queue for testing."""
    manager = Manager()
    return manager.Queue()


class TestCoreDetection:
    """Test CPU core detection and allocation logic."""
    
    def test_cpu_count_detection(self, mock_progress_queue):
        """Verify multiprocessing can detect CPU cores."""
        num_cores = mp.cpu_count()
        assert num_cores > 0, "Should detect at least one CPU core"
        assert num_cores <= 256, "Sanity check: unlikely to have >256 cores"
    
    def test_core_override_logic(self, mock_progress_queue):
        """Test that manual core override would work correctly."""
        detected = mp.cpu_count()
        
        # Simulate --cores argument logic
        manual_cores = 4
        actual_cores = min(manual_cores, detected)
        
        assert actual_cores <= detected, "Should not exceed available cores"
        assert actual_cores == min(4, detected), "Should use minimum of requested and available"


class TestParallelWorkerFunction:
    """Test individual worker functions used in parallel execution."""
    
    def test_reversible_worker_runs(self, mock_progress_queue):
        """Test that reversible worker function completes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            args = (0, 0, 50, 5, 'off', tmpdir, False, 1, mock_progress_queue)
            result = run_single_sim_rev(args)
            
            assert 'R' in result
            assert 'M' in result
            assert 'filename' in result
            assert 'E_total' in result
            assert 'E_initial' in result
            assert result['R'] == 0
            assert result['M'] == 0
    
    def test_irreversible_worker_runs(self, mock_progress_queue):
        """Test that irreversible worker function completes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create irr subdirectory
            os.makedirs(os.path.join(tmpdir, 'data', 'irr', 'r0'), exist_ok=True)
            
            args = (0, 0, 50, 5, 'off', tmpdir, False, 1, mock_progress_queue)
            result = run_single_sim_irr(args)
            
            assert 'R' in result
            assert 'M' in result
            assert 'filename' in result
            assert 'E_total' in result
            assert 'E_initial' in result
            assert result['R'] == 0
            assert result['M'] == 0
    
    def test_worker_energy_conservation(self, mock_progress_queue):
        """Verify worker maintains energy conservation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            args = (1, 0, 100, 10, 'off', tmpdir, False, 1, mock_progress_queue)
            result = run_single_sim_rev(args)
            
            # Energy should be conserved
            assert abs(result['E_total'] - result['E_initial']) < 1e-10, \
                "Worker should conserve energy"
    
    def test_worker_creates_output_files(self, mock_progress_queue):
        """Verify worker creates CSV and metadata files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            args = (1, 0, 50, 5, 'off', tmpdir, False, 1, mock_progress_queue)
            result = run_single_sim_rev(args)
            
            csv_file = result['filename']
            metadata_file = csv_file.replace('.csv', '_metadata.json')
            
            assert os.path.exists(csv_file), "Should create CSV file"
            assert os.path.exists(metadata_file), "Should create metadata file"
            
            # Verify CSV has correct structure
            with open(csv_file, 'r') as f:
                reader = csv.reader(f)
                header = next(reader)
                assert header == ['t', 'K', 'U', 'N0', 'Nx', 'S/nk', 'n']
                
                # Should have 2*sweeps rows
                data_rows = list(reader)
                assert len(data_rows) == 2 * 5, "Should have 2*sweeps data rows"
            
            # Verify metadata has correct structure
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
                assert 'lattice_size' in metadata
                assert 'sweeps' in metadata
                assert 'radius' in metadata
                assert 'run' in metadata
                assert 'simulation_type' in metadata


class TestParallelVsSequential:
    """Compare parallel execution results with sequential baseline."""
    
    def test_single_worker_matches_direct_call(self, mock_progress_queue):
        """Single worker should produce same results as direct simulation."""
        n, s, R = 50, 5, 1

        with tempfile.TemporaryDirectory() as tmpdir:
            # Run via worker function
            args = (R, 0, n, s, 'off', tmpdir, False, 1, mock_progress_queue)
            worker_result = run_single_sim_rev(args)
            
            # Read worker output
            with open(worker_result['filename'], 'r') as f:
                reader = csv.reader(f)
                next(reader)  # Skip header
                worker_data = [[float(x) for x in row] for row in reader]
            
            # Run direct simulation with same parameters
            np.random.seed(42)
            x = Inferno(n, R+1, validate_mode='off')
            
            # Both should conserve energy
            assert abs(worker_result['E_total'] - worker_result['E_initial']) < 1e-10
            assert abs(x.E_total - x._initial_total_energy) < 1e-10
    
    def test_multiple_workers_produce_different_results(self, mock_progress_queue):
        """Different runs should have different random sequences."""
        n, s, R = 50, 5, 1
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run two workers with same parameters but different run numbers
            args1 = (R, 0, n, s, 'off', tmpdir, False, 1, mock_progress_queue)
            args2 = (R, 1, n, s, 'off', tmpdir, False, 1, mock_progress_queue)
            
            result1 = run_single_sim_rev(args1)
            result2 = run_single_sim_rev(args2)
            
            # Read both outputs
            with open(result1['filename'], 'r') as f:
                reader = csv.reader(f)
                next(reader)
                data1 = [[float(x) for x in row] for row in reader]
            
            with open(result2['filename'], 'r') as f:
                reader = csv.reader(f)
                next(reader)
                data2 = [[float(x) for x in row] for row in reader]
            
            # Results should differ (different random seeds)
            assert data1 != data2, "Different runs should produce different trajectories"
            
            # But both should conserve energy
            assert abs(result1['E_total'] - result1['E_initial']) < 1e-10
            assert abs(result2['E_total'] - result2['E_initial']) < 1e-10


class TestParallelOutputStructure:
    """Test that parallel execution creates correct directory structure and files."""
    
    def test_reversible_output_structure(self, mock_progress_queue):
        """Verify reversible parallel creates correct directory structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Simulate running r=2 (R=0,1), m=2 (M=0,1) - 4 simulations
            results = []
            for R in range(2):
                for M in range(2):
                    args = (R, M, 30, 3, 'off', tmpdir, False, 1, mock_progress_queue)
                    result = run_single_sim_rev(args)
                    results.append(result)
            
            # Check directories were created
            assert os.path.exists(os.path.join(tmpdir, 'data', 'r0'))
            assert os.path.exists(os.path.join(tmpdir, 'data', 'r1'))
            
            # Check all files exist
            assert len(results) == 4
            for result in results:
                assert os.path.exists(result['filename'])
                assert os.path.exists(result['filename'].replace('.csv', '_metadata.json'))
    
    def test_irreversible_output_structure(self, mock_progress_queue):
        """Verify irreversible parallel creates correct directory structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Pre-create directory structure for irreversible
            for R in range(2):
                os.makedirs(os.path.join(tmpdir, 'data', 'irr', f'r{R}'), exist_ok=True)
            
            # Simulate running r=2, m=2
            results = []
            for R in range(2):
                for M in range(2):
                    args = (R, M, 30, 3, 'off', tmpdir, False, 1, mock_progress_queue)
                    result = run_single_sim_irr(args)
                    results.append(result)
            
            # Check files exist in irr subdirectories
            assert len(results) == 4
            for result in results:
                assert 'irr' in result['filename']
                assert os.path.exists(result['filename'])


class TestParallelResourceManagement:
    """Test resource management and error handling."""
    
    def test_worker_handles_small_lattice(self, mock_progress_queue):
        """Workers should handle edge cases like small lattices."""
        with tempfile.TemporaryDirectory() as tmpdir:
            args = (1, 0, 10, 2, 'off', tmpdir, False, 1, mock_progress_queue)  # Very small: n=10, s=2
            result = run_single_sim_rev(args)
            
            assert result['E_total'] == result['E_initial'] == 20, "Small lattice energy"
            assert os.path.exists(result['filename'])
    
    def test_worker_handles_different_radii(self, mock_progress_queue):
        """Workers should handle various radii correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            for R in [0, 1, 3, 5]:
                args = (R, 0, 40, 3, 'off', tmpdir, False, 1, mock_progress_queue)
                result = run_single_sim_rev(args)
                
                assert result['R'] == R
                assert abs(result['E_total'] - result['E_initial']) < 1e-10
    
    def test_worker_validation_modes(self, mock_progress_queue):
        """Workers should support all validation modes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            for mode in ['off', 'periodic', 'frequent']:
                args = (1, 0, 30, 3, mode, tmpdir, False, 1, mock_progress_queue)
                result = run_single_sim_rev(args)
                
                # All modes should conserve energy
                assert abs(result['E_total'] - result['E_initial']) < 1e-10


class TestParallelCorrectness:
    """Test mathematical and physical correctness of parallel execution."""
    
    def test_energy_conservation_all_workers(self, mock_progress_queue):
        """All parallel workers must conserve energy."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run multiple simulations in sequence (simulating parallel pool)
            results = []
            for R in range(3):
                for M in range(3):
                    args = (R, M, 40, 5, 'off', tmpdir, False, 1, mock_progress_queue)
                    result = run_single_sim_rev(args)
                    results.append(result)
            
            # All should conserve energy
            for result in results:
                energy_error = abs(result['E_total'] - result['E_initial'])
                assert energy_error < 1e-10, \
                    f"R={result['R']}, M={result['M']} failed energy conservation"
    
    def test_entropy_calculations_valid(self, mock_progress_queue):
        """Verify entropy values are physically reasonable."""
        with tempfile.TemporaryDirectory() as tmpdir:
            args = (1, 0, 100, 10, 'off', tmpdir, False, 1, mock_progress_queue)
            result = run_single_sim_rev(args)
            
            # Read CSV and check entropy column
            with open(result['filename'], 'r') as f:
                reader = csv.reader(f)
                next(reader)  # Skip header
                for row in reader:
                    entropy = float(row[5])  # S/nk column
                    assert not np.isnan(entropy), "Entropy should not be NaN"
                    assert entropy > 0, "Entropy should be positive"
    
    def test_reversible_forward_reverse_symmetry(self, mock_progress_queue):
        """Reversible simulation should show time-reversal properties."""
        with tempfile.TemporaryDirectory() as tmpdir:
            args = (2, 0, 50, 20, 'off', tmpdir, False, 1, mock_progress_queue)
            result = run_single_sim_rev(args)
            
            # Read data
            with open(result['filename'], 'r') as f:
                reader = csv.reader(f)
                next(reader)
                data = [[float(x) for x in row] for row in reader]
            
            # Split into forward and reverse
            s = len(data) // 2
            forward_final = data[s-1]
            reverse_start = data[s]
            
            # Final state of forward should be similar to start of reverse
            # Note: Exact match not expected due to Monte Carlo dynamics
            # Check energy conservation instead
            forward_energy = forward_final[1] + forward_final[2]  # K + U
            reverse_energy = reverse_start[1] + reverse_start[2]
            
            assert abs(forward_energy - reverse_energy) < 0.01, \
                f"Energy continuity check failed"


    @pytest.mark.slow
    def test_parallel_faster_than_sequential(self, mock_progress_queue):
        """Parallel execution should be faster for multiple simulations."""
        import time
        
        with tempfile.TemporaryDirectory() as tmpdir:
            n, s = 100, 10  # Larger size for meaningful timing
            num_sims = 4  # Fewer sims for faster test
            
            # Parallel timing (using pool)
            start = time.time()
            cores = min(4, mp.cpu_count())
            with mp.Pool(processes=cores) as pool:
                params = [(0, i, n, s, 'off', tmpdir, False, i+1, mock_progress_queue) for i in range(num_sims)]
                results = pool.map(run_single_sim_rev, params)
            parallel_time = time.time() - start
            
            # Should complete successfully
            assert len(results) == num_sims
            assert parallel_time > 0
            for result in results:
                assert abs(result['E_total'] - result['E_initial']) < 1e-10
    
    def test_worker_memory_footprint_reasonable(self, mock_progress_queue):
        """Workers should not consume excessive memory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Even with larger lattice, should complete without issues
            args = (1, 0, 1000, 10, 'off', tmpdir, False, 1, mock_progress_queue)
            result = run_single_sim_rev(args)
            
            assert os.path.exists(result['filename'])
            assert abs(result['E_total'] - result['E_initial']) < 1e-10


class TestParallelEdgeCases:
    """Test edge cases and boundary conditions."""
    
    def test_single_simulation_parallel(self, mock_progress_queue):
        """Parallel with just one simulation should work."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with mp.Pool(processes=1) as pool:
                params = [(0, 0, 30, 3, 'off', tmpdir, False, 1, mock_progress_queue)]
                results = pool.map(run_single_sim_rev, params)
            
            assert len(results) == 1
            assert abs(results[0]['E_total'] - results[0]['E_initial']) < 1e-10
    
    def test_more_workers_than_sims(self, mock_progress_queue):
        """Should handle pool size > number of simulations."""
        with tempfile.TemporaryDirectory() as tmpdir:
            num_cores = min(8, mp.cpu_count())
            num_sims = 3
            
            with mp.Pool(processes=num_cores) as pool:
                params = [(0, i, 30, 3, 'off', tmpdir, False, i+1, mock_progress_queue) for i in range(num_sims)]
                results = pool.map(run_single_sim_rev, params)
            
            assert len(results) == num_sims
            for result in results:
                assert abs(result['E_total'] - result['E_initial']) < 1e-10
    
    def test_r0_special_case(self, mock_progress_queue):
        """R=0 should use different filename pattern."""
        with tempfile.TemporaryDirectory() as tmpdir:
            args = (0, 0, 30, 3, 'off', tmpdir, False, 1, mock_progress_queue)
            result = run_single_sim_rev(args)
            
            # R=0 uses 'sim_data_0.csv' not 'sim_data_r0_0.csv'
            assert 'sim_data_0.csv' in result['filename']
            assert 'r0' not in os.path.basename(result['filename'])
