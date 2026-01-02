"""
Tests for parallel_sim.py - reversible parallel simulation runner.

Focuses on testing exception handling, CLI entry point, and integration
with the run_single_simulation function.
"""

import json
import multiprocessing as mp
import os
import sys

# Add creutz-sim to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

from parallel_sim import run_single_simulation


class TestRunSingleSimulation:
    """Test the run_single_simulation function."""

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

        # Verify output file was created
        assert os.path.exists(result["filename"])

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
            assert metadata["simulation_type"] == "reversible_parallel"

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

        # Verify R=0 directory structure
        expected_dir = os.path.join(project_root, "data", "r0")
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


class TestCLIComponents:
    """Test CLI components and argument parsing."""

    def test_argument_parser_basic_args(self):
        """Test that argument parser correctly parses basic arguments."""
        from sim_utils import create_argument_parser

        parser = create_argument_parser("reversible")
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

        parser = create_argument_parser("reversible")
        args = parser.parse_args(
            ["--n", "50", "--s", "3", "--r", "1", "--m", "1", "--no-jit"]
        )

        assert args.no_jit is True
        assert args.n == 50

    def test_argument_parser_with_cores(self):
        """Test argument parser with --cores argument."""
        from sim_utils import create_argument_parser

        parser = create_argument_parser("reversible")
        args = parser.parse_args(
            ["--n", "50", "--s", "2", "--r", "1", "--m", "2", "--cores", "2"]
        )

        assert args.cores == 2

    def test_argument_parser_with_validation(self):
        """Test argument parser with --validate argument."""
        from sim_utils import create_argument_parser

        parser = create_argument_parser("reversible")
        args = parser.parse_args(
            ["--n", "50", "--s", "2", "--r", "1", "--m", "1", "--validate", "periodic"]
        )

        assert args.validate == "periodic"

    def test_simulation_parameters_build(self, tmp_path):
        """Test that simulation parameters are built correctly."""
        from sim_utils import build_simulation_parameters

        project_root = str(tmp_path)
        r = 2  # Will test R=1 and R=2
        m = 3  # 3 runs each
        n = 100
        s = 10
        validate_mode = "off"
        use_jit = True

        params = build_simulation_parameters(
            r, m, n, s, validate_mode, project_root, use_jit
        )

        # Should have r * m simulations
        assert len(params) == r * m

        # Check first and last param structure
        first_param = params[0]
        assert first_param[0] == 0  # R starts at 0
        assert first_param[1] == 0  # M starts at 0
        assert first_param[2] == n
        assert first_param[3] == s
        assert first_param[4] == validate_mode
        assert first_param[5] == project_root
        assert first_param[6] == use_jit
        assert first_param[7] == 1  # sim_num starts at 1

    def test_data_directory_creation(self, tmp_path):
        """Test that data directories are created correctly."""
        from sim_utils import create_data_directories

        project_root = str(tmp_path)
        r = 3

        create_data_directories(project_root, r, "reversible")

        # Check that directories were created
        data_dir = tmp_path / "data"
        assert data_dir.exists()

        # For reversible, should create r0, r1, r2 (R=0 to R=r-1)
        for i in range(r):
            assert (data_dir / f"r{i}").exists()

    def test_data_directory_creation_irreversible(self, tmp_path):
        """Test that irreversible directories are created in irr/ subdirectory."""
        from sim_utils import create_data_directories

        project_root = str(tmp_path)
        r = 2

        create_data_directories(project_root, r, "irreversible")

        # Check that irr/ directories were created
        data_dir = tmp_path / "data" / "irr"
        assert data_dir.exists()

        # Should create r0, r1 (R=0 to R=r-1)
        for i in range(r):
            assert (data_dir / f"r{i}").exists()
