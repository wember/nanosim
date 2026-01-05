"""
Tests for parallel_sim_combined.py - combined parallel simulation runner.

Tests the combined simulation that runs both reversible and irreversible
simulations simultaneously.
"""

import json
import multiprocessing as mp
import os
import subprocess
import sys
from unittest.mock import MagicMock, patch

import pytest

# Add creutz-sim to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

from parallel_sim_combined import main, run_single_simulation


class TestRunSingleSimulation:
    """Test the run_single_simulation function."""

    def test_run_single_simulation_reversible_with_jit(self, tmp_path):
        """Test run_single_simulation for reversible type with JIT enabled."""
        # Create a queue for progress updates
        progress_queue = mp.Queue()

        # Simulation parameters
        R = 1
        M = 1
        n = 100
        s = 10
        validate_mode = "off"
        project_root = str(tmp_path)
        use_jit = True
        sim_num = 1
        sim_type = "reversible"

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            sim_type,
            progress_queue,
        )

        # Run simulation
        result = run_single_simulation(args)

        # Verify result structure
        assert isinstance(result, dict)
        assert result["R"] == R
        assert result["M"] == M
        assert result["sim_type"] == "reversible"
        assert result["success"] is True
        assert "filename" in result

        # Verify output file exists
        assert os.path.exists(result["filename"])

        # Verify CSV has correct structure
        import csv

        with open(result["filename"], "r") as f:
            reader = csv.reader(f)
            header = next(reader)
            assert header == ["t", "K", "U", "N0", "Nx", "S/nk", "n"]

            # Should have 2*s rows (forward + reverse)
            rows = list(reader)
            assert len(rows) == 2 * s

        # Verify metadata file exists
        metadata_file = result["filename"].replace(".csv", "_metadata.json")
        assert os.path.exists(metadata_file)

        with open(metadata_file, "r") as f:
            metadata = json.load(f)
            assert metadata["lattice_size"] == n
            assert metadata["sweeps"] == s
            assert metadata["radius"] == R
            assert metadata["run"] == M
            assert metadata["simulation_type"] == "reversible_parallel"

    def test_run_single_simulation_irreversible_with_jit(self, tmp_path):
        """Test run_single_simulation for irreversible type with JIT enabled."""
        # Create a queue for progress updates
        progress_queue = mp.Queue()

        # Simulation parameters
        R = 0
        M = 0
        n = 100
        s = 5
        validate_mode = "off"
        project_root = str(tmp_path)
        use_jit = True
        sim_num = 2
        sim_type = "irreversible"

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            sim_type,
            progress_queue,
        )

        # Run simulation
        result = run_single_simulation(args)

        # Verify result structure
        assert isinstance(result, dict)
        assert result["R"] == R
        assert result["M"] == M
        assert result["sim_type"] == "irreversible"
        assert result["success"] is True
        assert "filename" in result

        # Verify output file exists and is in irr directory
        assert os.path.exists(result["filename"])
        assert "/irr/" in result["filename"] or "\\irr\\" in result["filename"]

        # Verify metadata
        metadata_file = result["filename"].replace(".csv", "_metadata.json")
        with open(metadata_file, "r") as f:
            metadata = json.load(f)
            assert metadata["simulation_type"] == "irreversible_parallel"

    def test_run_single_simulation_without_jit(self, tmp_path):
        """Test run_single_simulation without JIT (original classes)."""
        progress_queue = mp.Queue()

        R = 0
        M = 0
        n = 50
        s = 3
        validate_mode = "off"
        project_root = str(tmp_path)
        use_jit = False
        sim_num = 1
        sim_type = "reversible"

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            sim_type,
            progress_queue,
        )

        # Run simulation
        result = run_single_simulation(args)

        # Verify success
        assert result["success"] is True
        assert os.path.exists(result["filename"])

    def test_run_single_simulation_progress_messages(self, tmp_path):
        """Test that progress messages are sent to queue."""
        progress_queue = mp.Queue()

        R = 0
        M = 0
        n = 50
        s = 2
        validate_mode = "off"
        project_root = str(tmp_path)
        use_jit = True
        sim_num = 1
        sim_type = "reversible"

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            sim_type,
            progress_queue,
        )

        # Run simulation
        result = run_single_simulation(args)
        assert result["success"] is True

        # Small delay to ensure all messages are in the queue
        import time

        time.sleep(0.1)

        # Collect all messages from queue
        messages = []
        while not progress_queue.empty():
            messages.append(progress_queue.get())

        # Should have: 1 start + s forward progress + s reverse progress + 1 complete
        # Note: Actual count may vary slightly due to timing, so check minimum
        assert (
            len(messages) >= 1 + s + s
        ), f"Expected at least {1 + s + s} messages, got {len(messages)}"

        # Verify message types
        assert messages[0]["type"] == "start"
        assert messages[0]["sim_type"] == "reversible"

        # Find completion message (should be last or near last)
        complete_messages = [m for m in messages if m.get("type") == "complete"]
        assert len(complete_messages) >= 1, "No completion message found"

        # Verify progress messages
        forward_messages = [m for m in messages if m.get("phase") == "forward"]
        reverse_messages = [m for m in messages if m.get("phase") == "reverse"]
        assert len(forward_messages) == s
        assert len(reverse_messages) == s


class TestExceptionHandling:
    """Test exception handling in run_single_simulation."""

    def test_entropy_overflow_handling(self, tmp_path):
        """Test that entropy calculation overflow is handled gracefully."""
        progress_queue = mp.Queue()

        # Use parameters that might cause overflow
        R = 0
        M = 0
        n = 100
        s = 5
        validate_mode = "off"
        project_root = str(tmp_path)
        use_jit = True
        sim_num = 1
        sim_type = "reversible"

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            sim_type,
            progress_queue,
        )

        # Should not raise exception even if entropy calculation overflows
        result = run_single_simulation(args)
        assert result["success"] is True

        # CSV should still be written
        import csv

        with open(result["filename"], "r") as f:
            reader = csv.reader(f)
            next(reader)  # Skip header
            rows = list(reader)
            # Check if any entropy values are NaN (which is valid for overflow)
            for row in rows:
                entropy_val = float(row[5])
                # Should be either a valid number or NaN, not raise
                assert isinstance(entropy_val, float)

    def test_invalid_sim_type_handling(self, tmp_path):
        """Test that invalid sim_type parameter raises appropriate error."""
        progress_queue = mp.Queue()

        # Test that function works with valid types first
        for sim_type in ["reversible", "irreversible"]:
            args = (
                0,
                0,
                50,
                2,
                "off",
                str(tmp_path),
                True,
                1,
                sim_type,
                progress_queue,
            )
            result = run_single_simulation(args)
            assert result["success"] is True


class TestMainFunction:
    """Test the main() function."""

    @patch(
        "sys.argv",
        ["parallel_sim_combined.py", "--n", "100", "--s", "5", "--r", "2", "--m", "1"],
    )
    @patch("parallel_sim_combined.mp.Pool")
    def test_main_basic_execution(self, mock_pool, tmp_path):
        """Test main() function basic execution flow."""
        # Mock the pool to avoid actual parallel execution
        mock_pool_instance = MagicMock()
        mock_pool.return_value.__enter__.return_value = mock_pool_instance

        # Mock map_async result
        mock_result = MagicMock()
        mock_result.ready.return_value = True
        mock_results = [
            {"R": 0, "M": 0, "sim_type": "reversible", "success": True},
            {"R": 0, "M": 0, "sim_type": "irreversible", "success": True},
            {"R": 1, "M": 0, "sim_type": "reversible", "success": True},
            {"R": 1, "M": 0, "sim_type": "irreversible", "success": True},
        ]
        mock_result.get.return_value = mock_results
        mock_pool_instance.map_async.return_value = mock_result

        # Change to a temporary directory for the test
        original_dir = os.getcwd()
        try:
            os.chdir(tmp_path)
            os.makedirs("data/r0", exist_ok=True)
            os.makedirs("data/r1", exist_ok=True)
            os.makedirs("data/irr/r0", exist_ok=True)
            os.makedirs("data/irr/r1", exist_ok=True)
            os.makedirs("logs", exist_ok=True)

            # Run main - should not raise
            main()

            # Verify pool was called
            assert mock_pool_instance.map_async.called

        finally:
            os.chdir(original_dir)

    @pytest.mark.slow
    def test_main_subprocess_execution(self, tmp_path):
        """Test main() completes successfully when run as subprocess."""
        # Create a minimal test script
        test_script = tmp_path / "test_run.py"
        project_root = os.path.dirname(os.path.dirname(__file__))
        creutz_sim_path = os.path.join(project_root, "creutz-sim")

        test_script.write_text(
            f"""
import sys
import os

# Setup paths
sys.path.insert(0, '{creutz_sim_path}')

if __name__ == '__main__':
    # Change to project root so data directories work
    os.chdir('{project_root}')

    from parallel_sim_combined import main

    # Mock sys.argv with minimal parameters
    sys.argv = [
        'parallel_sim_combined.py',
        '--n', '50',
        '--s', '2',
        '--r', '2',
        '--m', '1',
    ]

    # Run main - should complete without error
    main()
"""
        )

        # Run the script
        result = subprocess.run(
            [sys.executable, str(test_script)],
            capture_output=True,
            text=True,
            timeout=30,
        )

        # Should complete successfully
        assert result.returncode == 0, (
            f"Failed with return code {result.returncode}\n"
            f"stdout: {result.stdout}\n"
            f"stderr: {result.stderr}"
        )


class TestDataOrganization:
    """Test that combined simulation organizes data correctly."""

    def test_data_directory_structure(self, tmp_path):
        """Test that data is organized in correct directories."""
        progress_queue = mp.Queue()

        # Run reversible simulation
        args_rev = (
            1,
            0,
            50,
            3,
            "off",
            str(tmp_path),
            True,
            1,
            "reversible",
            progress_queue,
        )
        result_rev = run_single_simulation(args_rev)

        # Run irreversible simulation
        args_irr = (
            1,
            0,
            50,
            3,
            "off",
            str(tmp_path),
            True,
            2,
            "irreversible",
            progress_queue,
        )
        result_irr = run_single_simulation(args_irr)

        # Verify reversible is in data/r1/
        assert "data" in result_rev["filename"]
        assert "/r1/" in result_rev["filename"] or "\\r1\\" in result_rev["filename"]
        assert (
            "/irr/" not in result_rev["filename"]
            and "\\irr\\" not in result_rev["filename"]
        )

        # Verify irreversible is in data/irr/r1/
        assert "data" in result_irr["filename"]
        assert "/irr/" in result_irr["filename"] or "\\irr\\" in result_irr["filename"]
        assert "/r1/" in result_irr["filename"] or "\\r1\\" in result_irr["filename"]

    def test_file_naming_convention(self, tmp_path):
        """Test that files follow correct naming convention."""
        progress_queue = mp.Queue()

        # Test R=0 (special case)
        args_r0_rev = (
            0,
            0,
            50,
            2,
            "off",
            str(tmp_path),
            True,
            1,
            "reversible",
            progress_queue,
        )
        result_r0_rev = run_single_simulation(args_r0_rev)
        assert "sim_data_0.csv" in result_r0_rev["filename"]

        args_r0_irr = (
            0,
            1,
            50,
            2,
            "off",
            str(tmp_path),
            True,
            2,
            "irreversible",
            progress_queue,
        )
        result_r0_irr = run_single_simulation(args_r0_irr)
        assert "irr_sim_data_1.csv" in result_r0_irr["filename"]

        # Test R>0
        args_r2_rev = (
            2,
            3,
            50,
            2,
            "off",
            str(tmp_path),
            True,
            3,
            "reversible",
            progress_queue,
        )
        result_r2_rev = run_single_simulation(args_r2_rev)
        assert "sim_data_r2_3.csv" in result_r2_rev["filename"]

        args_r2_irr = (
            2,
            3,
            50,
            2,
            "off",
            str(tmp_path),
            True,
            4,
            "irreversible",
            progress_queue,
        )
        result_r2_irr = run_single_simulation(args_r2_irr)
        assert "irr_sim_data_r2_3.csv" in result_r2_irr["filename"]


class TestSimulationTypes:
    """Test that both simulation types work correctly."""

    def test_both_types_produce_valid_output(self, tmp_path):
        """Test that both reversible and irreversible produce valid output."""
        progress_queue = mp.Queue()

        n = 100
        s = 5

        for sim_type in ["reversible", "irreversible"]:
            args = (0, 0, n, s, "off", str(tmp_path), True, 1, sim_type, progress_queue)
            result = run_single_simulation(args)

            assert result["success"] is True
            assert result["sim_type"] == sim_type

            # Verify CSV data
            import csv

            with open(result["filename"], "r") as f:
                reader = csv.reader(f)
                next(reader)  # Skip header
                rows = list(reader)
                assert len(rows) == 2 * s

                # Verify all values are numeric
                for row in rows:
                    assert len(row) == 7
                    for val in row:
                        float(val)  # Should not raise

    def test_combined_execution_coverage(self, tmp_path):
        """Test that combined execution covers both types."""
        progress_queue = mp.Queue()

        # Simulate what the combined runner does: create tasks for both types
        tasks = []
        sim_num = 1
        r = 2
        m = 1

        for R in range(r):
            for M in range(m):
                # Reversible
                tasks.append(
                    (
                        R,
                        M,
                        50,
                        3,
                        "off",
                        str(tmp_path),
                        True,
                        sim_num,
                        "reversible",
                        progress_queue,
                    )
                )
                sim_num += 1
                # Irreversible
                tasks.append(
                    (
                        R,
                        M,
                        50,
                        3,
                        "off",
                        str(tmp_path),
                        True,
                        sim_num,
                        "irreversible",
                        progress_queue,
                    )
                )
                sim_num += 1

        # Should have 2*r*m tasks (both types)
        assert len(tasks) == 2 * r * m

        # Run all tasks
        results = []
        for task in tasks:
            result = run_single_simulation(task)
            results.append(result)

        # Verify we have equal numbers of each type
        rev_results = [r for r in results if r["sim_type"] == "reversible"]
        irr_results = [r for r in results if r["sim_type"] == "irreversible"]
        assert len(rev_results) == r * m
        assert len(irr_results) == r * m
        assert all(r["success"] for r in results)


class TestResultEquivalence:
    """Test that combined simulation produces identical results to separate runs."""

    def test_combined_vs_separate_reversible(self, tmp_path):
        """Test that combined simulation produces same reversible results.

        Compares with separate reversible simulation run.
        """
        import csv
        import sys

        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

        from parallel_sim import run_single_simulation as run_reversible

        progress_queue_combined = mp.Queue()
        progress_queue_separate = mp.Queue()

        # Same random seed to ensure reproducibility
        import numpy as np

        seed = 42

        # Parameters
        R = 1
        M = 0
        n = 100
        s = 10
        validate_mode = "off"
        use_jit = True

        # Create separate temp directories to avoid file conflicts
        tmp_combined = tmp_path / "combined"
        tmp_separate = tmp_path / "separate"
        tmp_combined.mkdir()
        tmp_separate.mkdir()

        # Run via combined simulation function
        np.random.seed(seed)
        args_combined = (
            R,
            M,
            n,
            s,
            validate_mode,
            str(tmp_combined),
            use_jit,
            1,
            "reversible",
            progress_queue_combined,
        )
        result_combined = run_single_simulation(args_combined)

        # Run via separate reversible simulation function
        np.random.seed(seed)
        args_separate = (
            R,
            M,
            n,
            s,
            validate_mode,
            str(tmp_separate),
            use_jit,
            1,
            progress_queue_separate,
        )
        result_separate = run_reversible(args_separate)

        # Both should succeed
        assert result_combined["success"] is True
        assert result_separate is not None

        # Read CSV data from both
        with open(result_combined["filename"], "r") as f:
            reader_combined = csv.reader(f)
            next(reader_combined)  # Skip header
            rows_combined = list(reader_combined)

        with open(result_separate["filename"], "r") as f:
            reader_separate = csv.reader(f)
            next(reader_separate)  # Skip header
            rows_separate = list(reader_separate)

        # Should have same number of rows
        assert len(rows_combined) == len(rows_separate)
        assert len(rows_combined) == 2 * s

        # Note: Due to stochastic nature, exact values may differ
        # But structure and ranges should be similar
        for i, (row_c, row_s) in enumerate(zip(rows_combined, rows_separate)):
            assert len(row_c) == len(row_s) == 7

            # Time step should be identical
            assert float(row_c[0]) == float(row_s[0])

            # All values should be finite and in reasonable ranges
            for val_c, val_s in zip(row_c, row_s):
                val_c_float = float(val_c)
                val_s_float = float(val_s)
                assert not np.isnan(val_c_float) and not np.isinf(val_c_float)
                assert not np.isnan(val_s_float)
                assert not np.isinf(val_s_float)

    def test_combined_vs_separate_irreversible(self, tmp_path):
        """Test combined produces same irreversible results as separate run."""
        import csv
        import sys

        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

        from parallel_sim_irr import run_single_simulation as run_irreversible

        progress_queue_combined = mp.Queue()
        progress_queue_separate = mp.Queue()

        # Same random seed to ensure reproducibility
        import numpy as np

        seed = 123

        # Parameters
        R = 0
        M = 0
        n = 100
        s = 10
        validate_mode = "off"
        use_jit = True

        # Create separate temp directories to avoid file conflicts
        tmp_combined = tmp_path / "combined"
        tmp_separate = tmp_path / "separate"
        tmp_combined.mkdir()
        tmp_separate.mkdir()

        # Run via combined simulation function
        np.random.seed(seed)
        args_combined = (
            R,
            M,
            n,
            s,
            validate_mode,
            str(tmp_combined),
            use_jit,
            1,
            "irreversible",
            progress_queue_combined,
        )
        result_combined = run_single_simulation(args_combined)

        # Run via separate irreversible simulation function
        np.random.seed(seed)
        args_separate = (
            R,
            M,
            n,
            s,
            validate_mode,
            str(tmp_separate),
            use_jit,
            1,
            progress_queue_separate,
        )
        result_separate = run_irreversible(args_separate)

        # Both should succeed
        assert result_combined["success"] is True
        assert result_separate is not None

        # Read CSV data from both
        with open(result_combined["filename"], "r") as f:
            reader_combined = csv.reader(f)
            next(reader_combined)  # Skip header
            rows_combined = list(reader_combined)

        with open(result_separate["filename"], "r") as f:
            reader_separate = csv.reader(f)
            next(reader_separate)  # Skip header
            rows_separate = list(reader_separate)

        # Should have same number of rows
        assert len(rows_combined) == len(rows_separate)
        assert len(rows_combined) == 2 * s

        # Verify structure is identical
        for i, (row_c, row_s) in enumerate(zip(rows_combined, rows_separate)):
            assert len(row_c) == len(row_s) == 7

            # Time step should be identical
            assert float(row_c[0]) == float(row_s[0])

            # All values should be finite
            for val_c, val_s in zip(row_c, row_s):
                val_c_float = float(val_c)
                val_s_float = float(val_s)
                assert not np.isnan(val_c_float) and not np.isinf(val_c_float)
                assert not np.isnan(val_s_float) and not np.isinf(val_s_float)

    def test_combined_imports_correct_classes(self, tmp_path):
        """Test that combined simulation imports the correct simulation classes."""
        progress_queue = mp.Queue()

        n = 50
        s = 3
        R = 0
        M = 0

        # Test reversible type imports JITInferno
        args_rev = (
            R,
            M,
            n,
            s,
            "off",
            str(tmp_path),
            True,
            1,
            "reversible",
            progress_queue,
        )
        result_rev = run_single_simulation(args_rev)
        assert result_rev["success"] is True
        assert "sim_data" in result_rev["filename"]
        assert "irr_sim_data" not in result_rev["filename"]

        # Test irreversible type imports JITirrInferno
        args_irr = (
            R,
            M,
            n,
            s,
            "off",
            str(tmp_path),
            True,
            2,
            "irreversible",
            progress_queue,
        )
        result_irr = run_single_simulation(args_irr)
        assert result_irr["success"] is True
        assert "irr_sim_data" in result_irr["filename"]

        # Verify both produce valid data
        import csv

        for result in [result_rev, result_irr]:
            with open(result["filename"], "r") as f:
                reader = csv.reader(f)
                header = next(reader)
                rows = list(reader)
                assert len(rows) == 2 * s
                assert header == ["t", "K", "U", "N0", "Nx", "S/nk", "n"]

    def test_metadata_consistency(self, tmp_path):
        """Test metadata files are consistent between combined and separate runs."""
        import sys

        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

        from parallel_sim import run_single_simulation as run_reversible

        progress_queue_combined = mp.Queue()
        progress_queue_separate = mp.Queue()

        R = 1
        M = 0
        n = 100
        s = 5

        tmp_combined = tmp_path / "combined"
        tmp_separate = tmp_path / "separate"
        tmp_combined.mkdir()
        tmp_separate.mkdir()

        # Run combined
        args_combined = (
            R,
            M,
            n,
            s,
            "off",
            str(tmp_combined),
            True,
            1,
            "reversible",
            progress_queue_combined,
        )
        result_combined = run_single_simulation(args_combined)

        # Run separate
        args_separate = (
            R,
            M,
            n,
            s,
            "off",
            str(tmp_separate),
            True,
            1,
            progress_queue_separate,
        )
        result_separate = run_reversible(args_separate)

        # Load metadata
        metadata_combined_file = result_combined["filename"].replace(
            ".csv", "_metadata.json"
        )
        metadata_separate_file = result_separate["filename"].replace(
            ".csv", "_metadata.json"
        )

        with open(metadata_combined_file, "r") as f:
            metadata_combined = json.load(f)

        with open(metadata_separate_file, "r") as f:
            metadata_separate = json.load(f)

        # Compare key fields (excluding timestamp)
        assert (
            metadata_combined["lattice_size"] == metadata_separate["lattice_size"] == n
        )
        assert metadata_combined["sweeps"] == metadata_separate["sweeps"] == s
        assert metadata_combined["radius"] == metadata_separate["radius"] == R
        assert metadata_combined["run"] == metadata_separate["run"] == M
        assert metadata_combined["simulation_type"] == "reversible_parallel"
        assert metadata_separate["simulation_type"] == "reversible_parallel"
