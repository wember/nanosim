"""
Copyright (c) 2026 Winry Ember
Licensed under the MIT License.
See LICENSE file in the project root for full license information.

Tests for sim_utils module functions.
"""

import os
import sys

import numpy as np

# Add creutz-sim directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))

import sim_utils


def test_Sk_stable_matches_manual():
    N, K = 10, 5
    # Manual calculation using loggamma
    expected = sim_utils.logg(K + N) - sim_utils.logg(K + 1) - sim_utils.logg(N)
    result = sim_utils.Sk_stable(N, K)
    assert np.isclose(result, expected)


def test_Su_stable_matches_manual():
    N, N0, Nx = 10, 3, 2
    N0_exp = max(N0, 1)
    log_2_term = N0_exp * np.log(2)
    expected = (
        sim_utils.logg(N + 1)
        + log_2_term
        - (
            sim_utils.logg(N - N0 - Nx + 1)
            + sim_utils.logg(N0 + 1)
            + sim_utils.logg(Nx + 1)
        )
    )
    result = sim_utils.Su_stable(N, N0, Nx, N0_exp)
    assert np.isclose(result, expected)


def test_create_argument_parser_defaults():
    parser = sim_utils.create_argument_parser("reversible")
    args = parser.parse_args([])
    assert args.n == 1000000
    assert args.s == 5000
    assert args.r == 11
    assert args.m == 5
    assert args.cores is None
    assert args.validate == "off"
    assert args.no_jit is False


def test_print_simulation_info_and_setup_logging(tmp_path, capsys):
    # First call setup_logging to create the directory
    sim_utils.setup_logging(str(tmp_path), "reversible")
    assert (tmp_path / "logs").exists()

    # Then test print_simulation_info
    sim_utils.print_simulation_info("reversible", 10, 5, 2, 1, 1, "off", True)
    captured = capsys.readouterr()
    assert "Parallel reversible simulation:" in captured.out


def test_create_data_directories(tmp_path):
    sim_utils.create_data_directories(str(tmp_path), 2, "reversible")
    assert (tmp_path / "data" / "r0").exists()
    assert (tmp_path / "data" / "r1").exists()
    sim_utils.create_data_directories(str(tmp_path), 2, "irreversible")
    assert (tmp_path / "data" / "irr" / "r0").exists()
    assert (tmp_path / "data" / "irr" / "r1").exists()


def test_build_simulation_parameters_shape():
    params = sim_utils.build_simulation_parameters(2, 2, 10, 5, "off", "/tmp", True)
    assert len(params) == 4
    for tup in params:
        assert len(tup) == 9


def test_format_time_and_elapsed_time():
    assert sim_utils.format_time(45) == "45 min"
    assert sim_utils.format_time(125) == "2h 5 min"
    assert sim_utils.format_time(1500) == "1d 1h 0 min"
    assert sim_utils.format_elapsed_time(45.2) == "45.2s"
    # Accept either possible rounding for 313.5 seconds
    assert sim_utils.format_elapsed_time(313.5) in ("5.22m (313.5s)", "5.23m (313.5s)")
    assert sim_utils.format_elapsed_time(4000) == "1.11h (66.7m)"


def test_SimulationBase_static_methods():
    N = 10
    order, rev_order = sim_utils.SimulationBase.initialize_order_arrays(N)
    assert set(order) == set(range(N))
    assert np.array_equal(rev_order, order[::-1])
    lattice, bonds, bond_count = sim_utils.SimulationBase.initialize_lattice(N)
    assert len(lattice) == N
    assert len(bonds) == N
    assert bond_count.shape == (3,)
    E_demon, E_demon_sum, d_energy = sim_utils.SimulationBase.initialize_demon_energy(
        N, 20, 10
    )
    assert E_demon.sum() == d_energy == E_demon_sum == 10
    right, left = sim_utils.SimulationBase.setup_neighbor_arrays(N)
    assert np.all(right == (np.arange(1, N + 1) % N))
    assert np.all(left == (np.arange(-1, N - 1) % N))


def test_create_data_directories_already_exists(tmp_path):
    """Test that create_data_directories works when directories already exist."""
    # Create directories first
    os.makedirs(tmp_path / "data" / "r0", exist_ok=True)
    os.makedirs(tmp_path / "data" / "r1", exist_ok=True)

    # Should not raise when directories already exist
    sim_utils.create_data_directories(str(tmp_path), 2, "reversible")
    assert (tmp_path / "data" / "r0").exists()
    assert (tmp_path / "data" / "r1").exists()

    # Same for irreversible
    os.makedirs(tmp_path / "data" / "irr" / "r0", exist_ok=True)
    sim_utils.create_data_directories(str(tmp_path), 2, "irreversible")
    assert (tmp_path / "data" / "irr" / "r0").exists()


def test_build_simulation_parameters_ordering():
    """Test that simulation parameters are ordered correctly (M then R)."""
    params = sim_utils.build_simulation_parameters(2, 3, 10, 5, "off", "/tmp", True)
    # Should have 2 radii * 3 runs = 6 total simulations
    assert len(params) == 6
    # Check ordering: should iterate M first, then R
    # M=0: R=0,1; M=1: R=0,1; M=2: R=0,1
    assert params[0][1] == 0  # M=0
    assert params[0][0] == 0  # R=0
    assert params[1][1] == 0  # M=0
    assert params[1][0] == 1  # R=1
    assert params[2][1] == 1  # M=1
    assert params[2][0] == 0  # R=0
    # Check sim_num increments properly
    assert params[0][7] == 1  # First sim
    assert params[5][7] == 6  # Last sim


def test_calculate_progress_description_early_stage():
    """Test progress description before 1% threshold."""
    active_sims = {
        1: {"R": 0, "M": 0, "phase": "forward", "progress": 5, "start_time": 0}
    }
    sim_times = []
    desc = sim_utils.calculate_progress_description(
        active_sims, sim_times, 1000, 100, 4, 0
    )
    # Should show status without ETA before 1% threshold
    assert "Overall Progress" in desc
    assert "R=0" in desc
    assert "M=0" in desc
    assert "forward" in desc
    assert "5/1000" in desc


def test_calculate_progress_description_with_eta():
    """Test progress description with ETA after 1% threshold."""
    active_sims = {
        1: {
            "R": 0,
            "M": 0,
            "phase": "reverse",
            "progress": 50,
            "start_time": 0,
            "phase_start_time": 0,
        }
    }
    # Add enough sim times to trigger ETA calculation
    sim_times = [10.0] * 5  # 5 completed sims
    desc = sim_utils.calculate_progress_description(
        active_sims, sim_times, 100, 100, 4, 5
    )
    assert "Overall Progress" in desc
    assert "Est:" in desc or "remaining" in desc


def test_handle_progress_message_throttling(monkeypatch):
    """Test that progress updates are throttled to 1 second intervals."""
    import time
    from unittest.mock import MagicMock

    active_sims = {
        1: {
            "R": 0,
            "M": 0,
            "phase": "forward",
            "progress": 5,
            "start_time": 0,
            "phase_start_time": 0,
        }
    }
    pbar = MagicMock()

    msg = {"sim_num": 1, "phase": "forward", "sweep": 10}

    # First call - should update display
    last_refresh = time.time() - 2.0  # 2 seconds ago
    new_refresh = sim_utils.handle_progress_message(
        msg, active_sims, 100, pbar, last_refresh, [], 100, 4, 0
    )
    assert new_refresh > last_refresh
    assert pbar.set_description.called

    # Second call immediately after - should NOT update display (throttled)
    pbar.reset_mock()
    newer_refresh = sim_utils.handle_progress_message(
        msg, active_sims, 100, pbar, new_refresh, [], 100, 4, 0
    )
    assert newer_refresh == new_refresh  # No update
    assert not pbar.set_description.called


def test_handle_progress_message_phase_transition():
    """Test that phase transitions are tracked."""
    import time
    from unittest.mock import MagicMock

    active_sims = {
        1: {
            "R": 0,
            "M": 0,
            "phase": "forward",
            "progress": 50,
            "start_time": 0,
            "phase_start_time": 0,
        }
    }
    pbar = MagicMock()

    # Message indicating phase transition from forward to reverse
    msg = {"sim_num": 1, "phase": "reverse", "sweep": 1}
    last_refresh = time.time() - 2.0

    sim_utils.handle_progress_message(
        msg, active_sims, 100, pbar, last_refresh, [], 100, 4, 0
    )

    # Phase should be updated
    assert active_sims[1]["phase"] == "reverse"
    assert active_sims[1]["progress"] == 1
    assert "phase_start_time" in active_sims[1]


def test_handle_progress_message_unknown_sim():
    """Test handling progress for non-existent simulation."""
    import time
    from unittest.mock import MagicMock

    active_sims = {}
    pbar = MagicMock()
    msg = {"sim_num": 999, "phase": "forward", "sweep": 10}
    last_refresh = time.time()

    # Should return unchanged when sim_num not in active_sims
    result = sim_utils.handle_progress_message(
        msg, active_sims, 100, pbar, last_refresh, [], 100, 4, 0
    )
    assert result == last_refresh
    assert not pbar.set_description.called


def test_setup_logging_creates_directory(tmp_path):
    """Test that setup_logging creates the logs directory."""
    log_dir = sim_utils.setup_logging(str(tmp_path), "reversible")
    assert os.path.exists(log_dir)
    assert log_dir == os.path.join(str(tmp_path), "logs")

    # Should work when called again (directory already exists)
    log_dir2 = sim_utils.setup_logging(str(tmp_path), "irreversible")
    assert os.path.exists(log_dir2)


def test_build_simulation_parameters_with_jit_disabled():
    """Test parameter building with JIT disabled."""
    params = sim_utils.build_simulation_parameters(
        2, 2, 100, 50, "periodic", "/test", False
    )
    assert len(params) == 4
    # Check that use_jit flag is False
    for p in params:
        assert p[6] is False  # use_jit is at index 6
        assert p[4] == "periodic"  # validate_mode at index 4
