"""
Tests for parallel_sim_irr.py - irreversible parallel simulation runner.

Focuses on testing exception handling, CLI entry point, and integration
with the run_single_simulation function.
"""

import json
import multiprocessing as mp
import os
import sys

# Add creutz-sim to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

from parallel_sim_irr import run_single_simulation


class TestRunSingleSimulationIrr:
    """Test the run_single_simulation function for irreversible simulation."""

    def test_run_single_simulation_with_jit(self, tmp_path):
        """Test run_single_simulation with JIT enabled."""
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

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            progress_queue,
        )

        # Run simulation
        result = run_single_simulation(args)

        # Verify result structure
        assert "R" in result
        assert "M" in result
        assert "filename" in result
        assert "E_total" in result
        assert "E_initial" in result
        assert result["R"] == R
        assert result["M"] == M

        # Verify output file was created in irr/ subdirectory
        assert os.path.exists(result["filename"])
        assert "irr" in result["filename"]

        # Verify CSV content
        import csv

        with open(result["filename"], "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            assert len(rows) == 2 * s
            assert all(
                key in rows[0] for key in ["t", "K", "U", "N0", "Nx", "S/nk", "n"]
            )

        # Verify metadata file was created
        metadata_file = result["filename"].replace(".csv", "_metadata.json")
        assert os.path.exists(metadata_file)

        with open(metadata_file, "r") as f:
            metadata = json.load(f)
            assert metadata["lattice_size"] == n
            assert metadata["sweeps"] == s
            assert metadata["radius"] == R
            assert metadata["run"] == M
            assert metadata["simulation_type"] == "irreversible_parallel"

        # Verify progress messages were sent
        # Add small delay to ensure queue is flushed
        import time

        time.sleep(0.01)

        messages = []
        while not progress_queue.empty():
            messages.append(progress_queue.get())

        # Should have: 1 start + 2*s progress + 1 complete
        assert len(messages) >= 3
        assert messages[0]["type"] == "start"

        # Check that complete message exists (might not be last due to queue ordering)
        message_types = [m["type"] for m in messages]
        assert "complete" in message_types

        # Check for both forward and reverse progress messages
        forward_msgs = [m for m in messages if m.get("phase") == "forward"]
        reverse_msgs = [m for m in messages if m.get("phase") == "reverse"]
        assert len(forward_msgs) == s
        assert len(reverse_msgs) == s

    def test_run_single_simulation_without_jit(self, tmp_path):
        """Test run_single_simulation with JIT disabled."""
        progress_queue = mp.Queue()

        R = 0
        M = 1
        n = 50
        s = 5
        validate_mode = "off"
        project_root = str(tmp_path)
        use_jit = False
        sim_num = 1

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            progress_queue,
        )

        result = run_single_simulation(args)

        # Verify result
        assert result["R"] == R
        assert result["M"] == M
        assert os.path.exists(result["filename"])

        # Verify R=0 directory structure in irr/
        expected_dir = os.path.join(project_root, "data", "irr", "r0")
        assert result["filename"].startswith(expected_dir)

    def test_run_single_simulation_with_validation(self, tmp_path):
        """Test run_single_simulation with validation enabled."""
        progress_queue = mp.Queue()

        R = 2
        M = 1
        n = 100
        s = 5
        validate_mode = "periodic"
        project_root = str(tmp_path)
        use_jit = True
        sim_num = 1

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            progress_queue,
        )

        result = run_single_simulation(args)

        # Verify energy conservation
        assert abs(result["E_total"] - result["E_initial"]) < 1e-10

    def test_run_single_simulation_entropy_overflow_handling(self, tmp_path):
        """Test that entropy overflow is handled gracefully with NaN."""
        progress_queue = mp.Queue()

        # Use parameters that might cause overflow
        R = 5
        M = 1
        n = 50
        s = 3
        validate_mode = "off"
        project_root = str(tmp_path)
        use_jit = True
        sim_num = 1

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            progress_queue,
        )

        # Should complete without raising exception
        result = run_single_simulation(args)

        # Verify CSV contains valid data (even if some entropy values are NaN)
        import csv

        with open(result["filename"], "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            assert len(rows) == 2 * s
            # At least some entropy values should be present
            entropy_values = [row["S/nk"] for row in rows]
            assert len(entropy_values) > 0

    def test_irreversible_output_directory_structure(self, tmp_path):
        """Test that irreversible simulation creates proper irr/ subdirectories."""
        progress_queue = mp.Queue()

        R = 3
        M = 2
        n = 50
        s = 3
        validate_mode = "off"
        project_root = str(tmp_path)
        use_jit = True
        sim_num = 1

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            progress_queue,
        )

        result = run_single_simulation(args)

        # Verify irr/r{R} directory structure
        expected_path = os.path.join(project_root, "data", "irr", f"r{R}")
        assert os.path.exists(expected_path)
        assert result["filename"].startswith(expected_path)


class TestCLIComponentsIrr:
    """Test CLI components and argument parsing for irreversible simulation."""

    def test_argument_parser_basic_args(self):
        """Test that argument parser correctly parses basic arguments."""
        from sim_utils import create_argument_parser

        parser = create_argument_parser("irreversible")
        args = parser.parse_args(["--n", "100", "--s", "5", "--r", "2", "--m", "1"])

        assert args.n == 100
        assert args.s == 5
        assert args.r == 2
        assert args.m == 1
        assert args.validate == "off"
        assert args.no_jit is False

    def test_argument_parser_with_no_jit(self):
        """Test argument parser with --no-jit flag."""
        from sim_utils import create_argument_parser

        parser = create_argument_parser("irreversible")
        args = parser.parse_args(
            ["--n", "50", "--s", "3", "--r", "1", "--m", "1", "--no-jit"]
        )

        assert args.no_jit is True
        assert args.n == 50

    def test_argument_parser_with_cores(self):
        """Test argument parser with --cores argument."""
        from sim_utils import create_argument_parser

        parser = create_argument_parser("irreversible")
        args = parser.parse_args(
            ["--n", "50", "--s", "2", "--r", "1", "--m", "2", "--cores", "2"]
        )

        assert args.cores == 2

    def test_argument_parser_with_validation(self):
        """Test argument parser with --validate argument."""
        from sim_utils import create_argument_parser

        parser = create_argument_parser("irreversible")
        args = parser.parse_args(
            ["--n", "50", "--s", "2", "--r", "1", "--m", "1", "--validate", "periodic"]
        )

        assert args.validate == "periodic"

    def test_simulation_parameters_build_for_irreversible(self, tmp_path):
        """Test that simulation parameters work for irreversible simulations."""
        from sim_utils import build_simulation_parameters

        project_root = str(tmp_path)
        r = 2
        m = 2
        n = 100
        s = 10
        validate_mode = "off"
        use_jit = True

        params = build_simulation_parameters(
            r, m, n, s, validate_mode, project_root, use_jit
        )

        # Should have r * m simulations
        assert len(params) == r * m

        # Parameters should be identical for both reversible and irreversible
        # (the difference is in which Inferno class is used)
        assert params[0][2] == n
        assert params[0][3] == s


class TestEntropyExceptionHandlingIrr:
    """Test exception handling in entropy calculations for irreversible simulations."""

    def test_entropy_calculation_with_overflow_irr(self, tmp_path, monkeypatch):
        """Test that entropy overflow errors are caught and handled gracefully."""
        import multiprocessing as mp

        # Mock Sk_stable to raise OverflowError
        def mock_Sk_overflow(n, E_demon_sum):
            raise OverflowError("Entropy calculation overflow")

        monkeypatch.setattr("parallel_sim_irr.Sk_stable", mock_Sk_overflow)

        progress_queue = mp.Queue()
        R = 0
        M = 0
        n = 50
        s = 5
        validate_mode = "off"
        project_root = str(tmp_path)
        use_jit = True
        sim_num = 1

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            progress_queue,
        )

        # Should handle overflow gracefully and set entropy to NaN
        result = run_single_simulation(args)

        # Verify simulation completed despite overflow
        assert result is not None
        assert os.path.exists(result["filename"])

        # Check CSV contains NaN for entropy
        import csv

        with open(result["filename"], "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            # At least some entropy values should be NaN
            nan_count = sum(1 for row in rows if row["S/nk"] == "nan")
            assert nan_count > 0

    def test_entropy_calculation_with_value_error_irr(self, tmp_path, monkeypatch):
        """Test that entropy ValueError is caught and handled gracefully."""
        import multiprocessing as mp

        # Mock Su_stable to raise ValueError
        def mock_Su_error(n, N0, Nx, N0_exp):
            raise ValueError("Invalid entropy parameters")

        monkeypatch.setattr("parallel_sim_irr.Su_stable", mock_Su_error)

        progress_queue = mp.Queue()
        R = 0
        M = 0
        n = 50
        s = 5
        validate_mode = "off"
        project_root = str(tmp_path)
        use_jit = True
        sim_num = 1

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            progress_queue,
        )

        # Should handle ValueError gracefully
        result = run_single_simulation(args)

        # Verify simulation completed
        assert result is not None
        assert os.path.exists(result["filename"])

    def test_entropy_both_phases_with_exceptions_irr(self, tmp_path, monkeypatch):
        """Test that exceptions in both forward and reverse phases are handled."""
        import multiprocessing as mp

        call_count = [0]

        def mock_Sk_alternating(n, E_demon_sum):
            call_count[0] += 1
            # Raise error on some calls to test both phases
            if call_count[0] % 3 == 0:
                raise OverflowError("Periodic overflow")
            return 1.0

        monkeypatch.setattr("parallel_sim_irr.Sk_stable", mock_Sk_alternating)

        progress_queue = mp.Queue()
        R = 0
        M = 0
        n = 30
        s = 10  # Enough sweeps to trigger multiple exceptions
        validate_mode = "off"
        project_root = str(tmp_path)
        use_jit = True
        sim_num = 1

        args = (
            R,
            M,
            n,
            s,
            validate_mode,
            project_root,
            use_jit,
            sim_num,
            progress_queue,
        )

        result = run_single_simulation(args)

        # Should complete with some NaN values
        assert result is not None
        assert os.path.exists(result["filename"])


class TestMainBlockExecutionIrr:
    """Test the __main__ block execution flow for irreversible simulations."""

    def test_main_block_setup_functions_irr(self, tmp_path):
        """Test that main block setup functions work correctly for irreversible."""
        import multiprocessing as mp

        from sim_utils import (
            build_simulation_parameters,
            create_argument_parser,
            create_data_directories,
            print_simulation_info,
            setup_logging,
        )

        # Test argument parser with typical args
        parser = create_argument_parser("irreversible")
        args = parser.parse_args(["--n", "100", "--s", "10", "--r", "2", "--m", "3"])

        n = args.n
        s = args.s
        r = args.r
        m = args.m
        validate_mode = args.validate
        use_jit = not args.no_jit

        # Test core count logic
        if args.cores is None:
            num_cores = mp.cpu_count()
        else:
            num_cores = min(args.cores, mp.cpu_count())

        assert num_cores > 0
        assert num_cores <= mp.cpu_count()

        # Test logging setup
        project_root = str(tmp_path)
        log_dir = setup_logging(project_root, "irreversible")
        assert os.path.exists(log_dir)

        # Test directory creation for irreversible (irr/ subdirectory)
        create_data_directories(project_root, r, "irreversible")
        for i in range(r):
            assert os.path.exists(os.path.join(project_root, "data", "irr", f"r{i}"))

        # Test parameter building
        sim_params = build_simulation_parameters(
            r, m, n, s, validate_mode, project_root, use_jit
        )
        assert len(sim_params) == r * m

        # Test that print_simulation_info runs without error
        print_simulation_info(
            "irreversible", n, s, r, m, num_cores, validate_mode, use_jit
        )

    def test_manager_and_queue_creation_irr(self):
        """Test that Manager and Queue can be created for parallel processing."""
        from multiprocessing import Manager

        manager = Manager()
        progress_queue = manager.Queue()

        # Verify queue works
        progress_queue.put({"type": "test", "data": 456})
        msg = progress_queue.get(timeout=1.0)
        assert msg["type"] == "test"
        assert msg["data"] == 456

    def test_script_execution_minimal_irr(self, tmp_path):
        """Test that the irreversible script can be executed directly with minimal
        parameters."""
        import subprocess

        script_path = os.path.join(
            os.path.dirname(__file__), "..", "creutz-sim", "parallel_sim_irr.py"
        )

        # Create data directory structure for irreversible
        data_dir = tmp_path / "data" / "irr" / "r0"
        data_dir.mkdir(parents=True)

        # Run with minimal parameters (very small simulation)
        result = subprocess.run(
            [
                sys.executable,
                script_path,
                "--n",
                "10",
                "--s",
                "2",
                "--r",
                "1",
                "--m",
                "1",
                "--cores",
                "1",
            ],
            cwd=str(tmp_path),
            capture_output=True,
            text=True,
            timeout=30,
        )

        # Check that script ran successfully
        assert result.returncode == 0, f"Script failed with: {result.stderr}"

        # Verify output was created in irr/ subdirectory
        assert (tmp_path / "data" / "irr" / "r0").exists()

        # Check for expected output messages
        output = result.stdout + result.stderr
        assert "simulation" in output.lower() or "complete" in output.lower()

    def test_cores_argument_capping_irr(self):
        """Test that cores argument is capped at cpu_count for irreversible."""
        import multiprocessing as mp

        from sim_utils import create_argument_parser

        parser = create_argument_parser("irreversible")

        # Try to request more cores than available
        max_cores = mp.cpu_count()
        args = parser.parse_args(
            [
                "--n",
                "50",
                "--s",
                "2",
                "--r",
                "1",
                "--m",
                "1",
                "--cores",
                str(max_cores + 10),
            ]
        )

        # The parallel_sim_irr main block should cap this at cpu_count
        if args.cores is None:
            num_cores = max_cores
        else:
            num_cores = min(args.cores, max_cores)

        assert num_cores == max_cores

    def test_directory_setup_and_logging_irr(self, tmp_path):
        """Test that directories and logging are set up correctly for irreversible."""
        from sim_utils import create_data_directories, setup_logging

        project_root = str(tmp_path)

        # Test logging setup
        log_dir = setup_logging(project_root, "irreversible")
        assert os.path.exists(log_dir)
        assert os.path.basename(log_dir) == "logs"

        # Test data directory creation for irreversible
        create_data_directories(project_root, 2, "irreversible")
        data_dir = tmp_path / "data" / "irr"
        assert data_dir.exists()
        assert (data_dir / "r0").exists()
        assert (data_dir / "r1").exists()
        assert (data_dir / "r0").exists()
        assert (data_dir / "r1").exists()
