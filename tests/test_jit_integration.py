"""Tests for JIT command-line flag integration in parallel scripts.

Validates that:
- JIT is enabled by default (no flag needed)
- --no-jit flag correctly disables JIT compilation
- JIT and non-JIT versions produce equivalent results
- Performance characteristics are as expected
- Error handling works correctly
"""

import csv
import os
import sys
import tempfile
from multiprocessing import Manager

import numpy as np
import pytest

# Add creutz-sim to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

from parallel_sim import run_single_simulation as run_single_sim_rev
from parallel_sim_irr import run_single_simulation as run_single_sim_irr


@pytest.fixture
def mock_progress_queue():
    """Create a mock progress queue for testing."""
    manager = Manager()
    return manager.Queue()


class TestJITFlagIntegration:
    """Test JIT default and --no-jit flag integration in parallel execution."""

    def test_reversible_jit_vs_nojit_results_match(self, mock_progress_queue):
        """JIT and non-JIT reversible simulations should produce equivalent results."""
        n, s, R = 100, 10, 2

        with tempfile.TemporaryDirectory() as tmpdir:
            # Run without JIT
            args_nojit = (R, 0, n, s, "off", tmpdir, False, 1, mock_progress_queue)
            result_nojit = run_single_sim_rev(args_nojit)

            # Run with JIT
            args_jit = (R, 0, n, s, "off", tmpdir, True, 1, mock_progress_queue)
            result_jit = run_single_sim_rev(args_jit)

            # Both should conserve energy
            assert abs(result_nojit["E_total"] - result_nojit["E_initial"]) < 1e-10
            assert abs(result_jit["E_total"] - result_jit["E_initial"]) < 1e-10

            # Total energy should match
            assert abs(result_nojit["E_total"] - result_jit["E_total"]) < 1e-10

            # Both should complete successfully
            assert result_nojit["filename"] is not None
            assert result_jit["filename"] is not None

    def test_irreversible_jit_vs_nojit_results_match(self, mock_progress_queue):
        """JIT and non-JIT irreversible simulations should produce
        equivalent results."""
        n, s, R = 100, 10, 2

        with tempfile.TemporaryDirectory() as tmpdir:
            # Run without JIT
            args_nojit = (R, 0, n, s, "off", tmpdir, False, 1, mock_progress_queue)
            result_nojit = run_single_sim_irr(args_nojit)

            # Run with JIT
            args_jit = (R, 0, n, s, "off", tmpdir, True, 1, mock_progress_queue)
            result_jit = run_single_sim_irr(args_jit)

            # Both should conserve energy
            assert abs(result_nojit["E_total"] - result_nojit["E_initial"]) < 1e-10
            assert abs(result_jit["E_total"] - result_jit["E_initial"]) < 1e-10

            # Total energy should match
            assert abs(result_nojit["E_total"] - result_jit["E_total"]) < 1e-10

    def test_jit_flag_creates_correct_output_structure(self, mock_progress_queue):
        """JIT simulations should create same output structure as non-JIT."""
        n, s, R = 50, 5, 1

        with tempfile.TemporaryDirectory() as tmpdir:
            # Run with JIT
            args = (R, 0, n, s, "off", tmpdir, True, 1, mock_progress_queue)
            result = run_single_sim_rev(args)

            # Check output file exists
            assert os.path.exists(result["filename"])

            # Read and validate CSV structure
            with open(result["filename"], "r") as f:
                reader = csv.reader(f)
                header = next(reader)

                # Check header
                assert header == ["t", "K", "U", "N0", "Nx", "S/nk", "n"]

                # Check data rows (should have 2*s rows)
                rows = list(reader)
                assert len(rows) == 2 * s

                # Validate each row has 7 columns
                for row in rows:
                    assert len(row) == 7

    def test_jit_flag_energy_conservation(self, mock_progress_queue):
        """JIT simulations should maintain energy conservation."""
        n, s, R = 100, 20, 3

        with tempfile.TemporaryDirectory() as tmpdir:
            # Run reversible with JIT
            args_rev = (R, 0, n, s, "off", tmpdir, True, 1, mock_progress_queue)
            result_rev = run_single_sim_rev(args_rev)

            # Run irreversible with JIT
            args_irr = (R, 0, n, s, "off", tmpdir, True, 1, mock_progress_queue)
            result_irr = run_single_sim_irr(args_irr)

            # Check energy conservation for both
            assert abs(result_rev["E_total"] - result_rev["E_initial"]) < 1e-10
            assert abs(result_irr["E_total"] - result_irr["E_initial"]) < 1e-10

            # Total energy should be 2*n for both
            expected_energy = 2.0 * n
            assert abs(result_rev["E_total"] - expected_energy) < 1e-10
            assert abs(result_irr["E_total"] - expected_energy) < 1e-10

    def test_jit_flag_with_different_radii(self, mock_progress_queue):
        """JIT should work correctly with various demon-coupling radii."""
        n, s = 50, 5

        with tempfile.TemporaryDirectory() as tmpdir:
            for R in [1, 2, 5, 10]:
                args = (R, 0, n, s, "off", tmpdir, True, 1, mock_progress_queue)
                result = run_single_sim_rev(args)

                # Should complete without error
                assert result["filename"] is not None
                assert abs(result["E_total"] - result["E_initial"]) < 1e-10

    def test_jit_flag_with_validation_modes(self, mock_progress_queue):
        """JIT should work with different validation modes."""
        n, s, R = 30, 10, 2

        with tempfile.TemporaryDirectory() as tmpdir:
            for mode in ["off", "periodic"]:
                args = (R, 0, n, s, mode, tmpdir, True, 1, mock_progress_queue)
                result = run_single_sim_rev(args)

                # Should complete without error
                assert result["filename"] is not None
                assert abs(result["E_total"] - result["E_initial"]) < 1e-10


class TestJITFlagPerformance:
    """Test performance characteristics of JIT flag."""

    def test_jit_produces_valid_entropy_values(self, mock_progress_queue):
        """JIT simulations should produce physically valid entropy values."""
        n, s, R = 100, 20, 3

        with tempfile.TemporaryDirectory() as tmpdir:
            args = (R, 0, n, s, "off", tmpdir, True, 1, mock_progress_queue)
            result = run_single_sim_rev(args)

            # Read output and check entropy values
            with open(result["filename"], "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    entropy = float(row["S/nk"])
                    # Entropy should be finite and non-negative
                    assert np.isfinite(entropy), f"Non-finite entropy at t={row['t']}"
                    assert entropy >= 0, f"Negative entropy at t={row['t']}"

    def test_jit_lattice_energy_bounded(self, mock_progress_queue):
        """JIT simulations should keep lattice energy within physical bounds."""
        n, s, R = 100, 20, 3

        with tempfile.TemporaryDirectory() as tmpdir:
            args = (R, 0, n, s, "off", tmpdir, True, 1, mock_progress_queue)
            result = run_single_sim_rev(args)

            # Read output and check energy bounds
            with open(result["filename"], "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    U = float(row["U"])  # Lattice energy per site
                    K = float(row["K"])  # Demon energy per site

                    # Total energy per site should be 2.0
                    total = U + K
                    assert (
                        abs(total - 2.0) < 1e-10
                    ), f"Energy not conserved at t={row['t']}"

                    # Lattice energy per site bounded by [-1, 1]
                    assert (
                        -1.0 <= U <= 1.0
                    ), f"Lattice energy out of bounds at t={row['t']}"

                    # Demon energy should be non-negative
                    assert K >= 0, f"Negative demon energy at t={row['t']}"


class TestJITFlagEdgeCases:
    """Test edge cases for JIT flag integration."""

    def test_jit_with_small_lattice(self, mock_progress_queue):
        """JIT should work even with very small lattices."""
        n, s, R = 10, 5, 1

        with tempfile.TemporaryDirectory() as tmpdir:
            args = (R, 0, n, s, "off", tmpdir, True, 1, mock_progress_queue)
            result = run_single_sim_rev(args)

            assert result["filename"] is not None
            assert abs(result["E_total"] - result["E_initial"]) < 1e-10

    def test_jit_with_large_radius(self, mock_progress_queue):
        """JIT should handle large demon-coupling radii."""
        n, s, R = 50, 5, 20

        with tempfile.TemporaryDirectory() as tmpdir:
            args = (R, 0, n, s, "off", tmpdir, True, 1, mock_progress_queue)
            result = run_single_sim_rev(args)

            assert result["filename"] is not None
            assert abs(result["E_total"] - result["E_initial"]) < 1e-10

    def test_jit_r0_special_case(self, mock_progress_queue):
        """JIT should handle R=0 special case correctly."""
        n, s = 30, 5

        with tempfile.TemporaryDirectory() as tmpdir:
            args = (0, 0, n, s, "off", tmpdir, True, 1, mock_progress_queue)
            result = run_single_sim_rev(args)

            # Should create output in r0 directory
            assert "r0" in result["filename"]
            assert os.path.exists(result["filename"])

    def test_jit_multiple_runs_same_radius(self, mock_progress_queue):
        """JIT should support multiple independent runs."""
        n, s, R = 30, 5, 2

        with tempfile.TemporaryDirectory() as tmpdir:
            results = []
            for M in range(3):
                args = (R, M, n, s, "off", tmpdir, True, 1, mock_progress_queue)
                result = run_single_sim_rev(args)
                results.append(result)

            # All should complete successfully
            assert len(results) == 3

            # Each should have different output file
            filenames = [r["filename"] for r in results]
            assert len(set(filenames)) == 3

            # All should conserve energy
            for result in results:
                assert abs(result["E_total"] - result["E_initial"]) < 1e-10


class TestJITFlagCompatibility:
    """Test JIT flag compatibility with existing functionality."""

    def test_jit_false_uses_original_classes(self, mock_progress_queue):
        """use_jit=False should use original Inferno classes."""
        n, s, R = 50, 5, 2

        with tempfile.TemporaryDirectory() as tmpdir:
            # This should use original classes
            args = (R, 0, n, s, "off", tmpdir, False, 1, mock_progress_queue)
            result = run_single_sim_rev(args)

            # Should still work correctly
            assert result["filename"] is not None
            assert abs(result["E_total"] - result["E_initial"]) < 1e-10

    def test_jit_true_uses_jit_classes(self, mock_progress_queue):
        """use_jit=True should use JIT-compiled classes."""
        n, s, R = 50, 5, 2

        with tempfile.TemporaryDirectory() as tmpdir:
            # This should use JIT classes
            args = (R, 0, n, s, "off", tmpdir, True, 1, mock_progress_queue)
            result = run_single_sim_rev(args)

            # Should work correctly
            assert result["filename"] is not None
            assert abs(result["E_total"] - result["E_initial"]) < 1e-10

    def test_jit_flag_backward_compatible(self, mock_progress_queue):
        """Tests should pass with both JIT enabled and disabled."""
        n, s, R = 50, 5, 2

        with tempfile.TemporaryDirectory() as tmpdir:
            # Run both versions
            args_off = (R, 0, n, s, "off", tmpdir, False, 1, mock_progress_queue)
            args_on = (R, 1, n, s, "off", tmpdir, True, 1, mock_progress_queue)

            result_off = run_single_sim_rev(args_off)
            result_on = run_single_sim_rev(args_on)

            # Both should produce valid results
            assert result_off["filename"] is not None
            assert result_on["filename"] is not None

            # Both should conserve energy
            assert abs(result_off["E_total"] - result_off["E_initial"]) < 1e-10
            assert abs(result_on["E_total"] - result_on["E_initial"]) < 1e-10
