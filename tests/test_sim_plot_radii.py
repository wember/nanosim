"""
Tests for refactored sim_plot_radii.py functions.

Tests the extracted testable functions directly.
"""

import os
import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

# Add creutz-sim to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))


class TestFindRadiusCsvFiles:
    """Test find_radius_csv_files function."""

    def test_find_files_for_multiple_radii(self, tmp_path):
        """Test finding CSV files across multiple radius directories."""
        from sim_plot_radii import find_radius_csv_files

        # Create data directories for R=0, R=1, R=2
        for R in [0, 1, 2]:
            data_dir = tmp_path / "data" / f"r{R}"
            data_dir.mkdir(parents=True)
            csv_file = data_dir / "sim_data.csv"
            csv_file.write_text("t,K,U,N0,Nx,n\n1,0.5,0.5,0.1,0.2,100\n")

        file_paths, display_names = find_radius_csv_files(str(tmp_path), max_radius=3)

        assert len(file_paths) == 3
        assert len(display_names) == 3
        assert display_names == ["radius = 0", "radius = 1", "radius = 2"]

    def test_no_files_found(self, tmp_path):
        """Test when no CSV files exist."""
        from sim_plot_radii import find_radius_csv_files

        file_paths, display_names = find_radius_csv_files(str(tmp_path))

        assert len(file_paths) == 0
        assert len(display_names) == 0

    def test_most_recent_file_selected(self, tmp_path):
        """Test that most recent CSV is selected when multiple exist."""
        from sim_plot_radii import find_radius_csv_files

        data_dir = tmp_path / "data" / "r0"
        data_dir.mkdir(parents=True)

        # Create two files with different timestamps
        old_file = data_dir / "sim_data_1.csv"
        old_file.write_text("t,K,U,N0,Nx,n\n1,0.5,0.5,0.1,0.2,100\n")

        import time

        time.sleep(0.01)

        new_file = data_dir / "sim_data_2.csv"
        new_file.write_text("t,K,U,N0,Nx,n\n1,0.5,0.5,0.1,0.2,100\n")

        file_paths, display_names = find_radius_csv_files(str(tmp_path), max_radius=1)

        assert len(file_paths) == 1
        assert file_paths[0] == str(new_file)

    def test_default_max_radius(self, tmp_path):
        """Test default max_radius value searches R=0 to R=10."""
        from sim_plot_radii import find_radius_csv_files

        # Create files for R=0 to R=10
        for R in range(11):
            data_dir = tmp_path / "data" / f"r{R}"
            data_dir.mkdir(parents=True)
            csv_file = data_dir / "sim_data.csv"
            csv_file.write_text("t,K,U,N0,Nx,n\n1,0.5,0.5,0.1,0.2,100\n")

        file_paths, display_names = find_radius_csv_files(str(tmp_path))

        assert len(file_paths) == 11
        assert len(display_names) == 11

    def test_missing_intermediate_radius(self, tmp_path, capsys):
        """Test handling when some radius directories are missing."""
        from sim_plot_radii import find_radius_csv_files

        # Create only R=0 and R=2, skip R=1
        for R in [0, 2]:
            data_dir = tmp_path / "data" / f"r{R}"
            data_dir.mkdir(parents=True)
            csv_file = data_dir / "sim_data.csv"
            csv_file.write_text("t,K,U,N0,Nx,n\n1,0.5,0.5,0.1,0.2,100\n")

        file_paths, display_names = find_radius_csv_files(str(tmp_path), max_radius=3)

        # Should find 2 files and print warning for R=1
        assert len(file_paths) == 2
        captured = capsys.readouterr()
        assert "Warning: No CSV files found" in captured.out


class TestCalculateEntropyWithN0Handling:
    """Test calculate_entropy_with_n0_handling function."""

    def test_basic_entropy_calculation(self):
        """Test entropy calculation with valid data."""
        from sim_plot_radii import calculate_entropy_with_n0_handling

        df = pd.DataFrame(
            {
                "t": [1, 2, 3],
                "K": [50.0, 55.0, 60.0],
                "U": [0.3, 0.4, 0.5],
                "N0": [10.0, 10.0, 10.0],
                "Nx": [20.0, 20.0, 20.0],
                "n": [100, 100, 100],
            }
        )

        entropy = calculate_entropy_with_n0_handling(df)

        assert len(entropy) == 3
        assert all(np.isfinite(entropy))
        assert all(entropy > 0)

    def test_n0_zero_handled(self):
        """Test N0=0 edge case is handled correctly."""
        from sim_plot_radii import calculate_entropy_with_n0_handling

        df = pd.DataFrame(
            {
                "t": [1, 2],
                "K": [50.0, 55.0],
                "U": [0.3, 0.4],
                "N0": [0.0, 10.0],
                "Nx": [20.0, 20.0],
                "n": [100, 100],
            }
        )

        entropy = calculate_entropy_with_n0_handling(df)

        # Should not have NaN or inf
        assert np.isfinite(entropy.iloc[0])
        assert np.isfinite(entropy.iloc[1])
        assert entropy.iloc[0] > 0

    def test_entropy_increases_with_demon_energy(self):
        """Test entropy increases with demon energy."""
        from sim_plot_radii import calculate_entropy_with_n0_handling

        df = pd.DataFrame(
            {
                "t": [1, 2, 3],
                "K": [10.0, 50.0, 90.0],
                "U": [0.3, 0.3, 0.3],
                "N0": [10.0, 10.0, 10.0],
                "Nx": [20.0, 20.0, 20.0],
                "n": [100, 100, 100],
            }
        )

        entropy = calculate_entropy_with_n0_handling(df)

        assert entropy.iloc[1] > entropy.iloc[0]
        assert entropy.iloc[2] > entropy.iloc[1]

    def test_entropy_with_large_lattice(self):
        """Test entropy calculation with large lattice size."""
        from sim_plot_radii import calculate_entropy_with_n0_handling

        df = pd.DataFrame(
            {
                "t": [1],
                "K": [500.0],
                "U": [0.5],
                "N0": [100.0],
                "Nx": [200.0],
                "n": [1000],
            }
        )

        entropy = calculate_entropy_with_n0_handling(df)

        assert np.isfinite(entropy.iloc[0])
        assert entropy.iloc[0] > 0


class TestAddRadiusTraceToSubplot:
    """Test add_radius_trace_to_subplot function."""

    def test_adds_trace_to_figure(self):
        """Test that function adds trace to subplot."""
        from plotly.subplots import make_subplots
        from sim_plot_radii import add_radius_trace_to_subplot

        fig = make_subplots(rows=1, cols=3)
        df = pd.DataFrame(
            {
                "t": [1, 2, 3],
                "K": [0.5, 0.6, 0.7],
                "U": [0.3, 0.4, 0.5],
            }
        )

        add_radius_trace_to_subplot(fig, df, "Test Radius", 1, 1, "K", "Demon Energy")

        # Should have 1 trace added
        assert len(fig.data) >= 1
        assert fig.data[0].name == "Test Radius"

    def test_adds_vline_at_midpoint(self):
        """Test that function adds vertical line at midpoint."""
        from plotly.subplots import make_subplots
        from sim_plot_radii import add_radius_trace_to_subplot

        fig = make_subplots(rows=1, cols=3)
        df = pd.DataFrame(
            {
                "t": list(range(1, 11)),
                "K": [0.5] * 10,
            }
        )

        add_radius_trace_to_subplot(fig, df, "Test", 1, 1, "K", "Energy")

        # Check for vline in layout shapes
        assert len(fig.layout.shapes) >= 1

    def test_multiple_columns(self):
        """Test adding traces to different columns."""
        from plotly.subplots import make_subplots
        from sim_plot_radii import add_radius_trace_to_subplot

        fig = make_subplots(rows=1, cols=3)
        df = pd.DataFrame(
            {
                "t": [1, 2, 3],
                "K": [0.5, 0.6, 0.7],
                "U": [0.3, 0.4, 0.5],
            }
        )

        add_radius_trace_to_subplot(fig, df, "Test", 1, 1, "K", "Demon Energy")
        add_radius_trace_to_subplot(fig, df, "Test", 1, 2, "U", "Lattice Temp")

        assert len(fig.data) >= 2


class TestCreateMultiRadiusPlot:
    """Test create_multi_radius_plot function."""

    def test_creates_figure_with_three_columns(self, tmp_path):
        """Test that function creates 1x3 subplot figure."""
        from sim_plot_radii import create_multi_radius_plot

        # Create test CSV files
        csv_file = tmp_path / "test_data.csv"
        csv_file.write_text(
            "t,K,U,N0,Nx,n\n" "1,0.5,0.3,0.1,0.2,100\n" "2,0.6,0.4,0.2,0.3,100\n"
        )

        file_names = [str(csv_file)]
        data_names = ["radius = 0"]

        fig = create_multi_radius_plot(file_names, data_names)

        # Should have traces for 3 subplots
        assert len(fig.data) == 3

    def test_handles_multiple_radii(self, tmp_path):
        """Test plotting multiple radii."""
        from sim_plot_radii import create_multi_radius_plot

        # Create CSV files for 3 radii
        file_names = []
        for i in range(3):
            csv_file = tmp_path / f"test_r{i}.csv"
            csv_file.write_text(
                "t,K,U,N0,Nx,n\n" "1,0.5,0.3,0.1,0.2,100\n" "2,0.6,0.4,0.2,0.3,100\n"
            )
            file_names.append(str(csv_file))

        data_names = [f"radius = {i}" for i in range(3)]

        fig = create_multi_radius_plot(file_names, data_names)

        # Should have 3 traces per subplot = 9 total
        assert len(fig.data) == 9

    def test_title_includes_lattice_size(self, tmp_path):
        """Test that title includes lattice size."""
        from sim_plot_radii import create_multi_radius_plot

        csv_file = tmp_path / "test_data.csv"
        csv_file.write_text(
            "t,K,U,N0,Nx,n\n"
            "1,0.5,0.3,0.1,0.2,500\n"  # n=500
        )

        fig = create_multi_radius_plot([str(csv_file)], ["radius = 0"])

        assert "500" in fig.layout.title.text
        assert "Multi-Radius Comparison" in fig.layout.title.text

    def test_entropy_calculation_applied(self, tmp_path):
        """Test that entropy is calculated and plotted."""
        from sim_plot_radii import create_multi_radius_plot

        csv_file = tmp_path / "test_data.csv"
        csv_file.write_text(
            "t,K,U,N0,Nx,n\n" "1,50.0,0.3,10.0,20.0,100\n" "2,55.0,0.4,10.0,20.0,100\n"
        )

        fig = create_multi_radius_plot([str(csv_file)], ["radius = 0"])

        # Third trace should be entropy
        entropy_trace = fig.data[2]
        assert len(entropy_trace.y) == 2
        assert all(np.isfinite(entropy_trace.y))


class TestMainFunction:
    """Test the main() function."""

    @patch("sim_plot_radii.sys.argv", ["sim_plot_radii.py"])
    @patch("sim_plot_radii.find_radius_csv_files")
    def test_no_files_found_exits(self, mock_find_files):
        """Test main exits when no files found."""
        from sim_plot_radii import main

        mock_find_files.return_value = ([], [])

        with pytest.raises(SystemExit) as exc_info:
            main()

        assert exc_info.value.code == 1

    @patch("sim_plot_radii.sys.argv", ["sim_plot_radii.py"])
    @patch("sim_plot_radii.find_radius_csv_files")
    @patch("sim_plot_radii.create_multi_radius_plot")
    @patch("plotly.graph_objects.Figure.show")
    def test_successful_execution(
        self, mock_show, mock_create_plot, mock_find_files, tmp_path
    ):
        """Test successful execution flow."""
        from sim_plot_radii import main

        # Mock file discovery
        csv_file = tmp_path / "test.csv"
        csv_file.write_text("t,K,U,N0,Nx,n\n1,0.5,0.3,0.1,0.2,100\n")

        mock_find_files.return_value = ([str(csv_file)], ["radius = 0"])

        # Mock plot creation
        mock_fig = go.Figure()
        mock_create_plot.return_value = mock_fig

        # Should not raise
        main()

        # Verify calls
        mock_find_files.assert_called_once()
        mock_create_plot.assert_called_once_with([str(csv_file)], ["radius = 0"])
        mock_show.assert_called_once()

    @patch("sim_plot_radii.sys.argv", ["sim_plot_radii.py"])
    @patch("sim_plot_radii.find_radius_csv_files")
    @patch("sim_plot_radii.create_multi_radius_plot")
    @patch("plotly.graph_objects.Figure.show")
    def test_prints_file_count(
        self, mock_show, mock_create_plot, mock_find_files, tmp_path, capsys
    ):
        """Test that file count is printed."""
        from sim_plot_radii import main

        csv_file = tmp_path / "test.csv"
        csv_file.write_text("t,K,U,N0,Nx,n\n1,0.5,0.3,0.1,0.2,100\n")

        mock_find_files.return_value = (
            [str(csv_file)] * 3,
            ["radius = 0", "radius = 1", "radius = 2"],
        )
        mock_fig = go.Figure()
        mock_create_plot.return_value = mock_fig

        main()

        captured = capsys.readouterr()
        assert "Found 3 CSV files to plot" in captured.out


class TestMainExecution:
    """Test __main__ block execution."""

    def test_main_function_is_callable(self, tmp_path, monkeypatch, capsys):
        """Test that main() function exists and is callable."""
        from sim_plot_radii import main

        # Mock find_radius_csv_files to return fake files
        with (
            patch("sim_plot_radii.find_radius_csv_files") as mock_find,
            patch("sim_plot_radii.create_multi_radius_plot") as mock_create,
        ):

            # Return some fake files so it doesn't exit early
            mock_find.return_value = (["file1.csv", "file2.csv"], ["R=0", "R=1"])
            mock_fig = MagicMock()
            mock_create.return_value = mock_fig

            # Call main - this tests the function itself
            main()

            # Verify it was called
            mock_find.assert_called_once()
            mock_create.assert_called_once()
            mock_fig.show.assert_called_once()

            # Check output
            captured = capsys.readouterr()
            assert "Found 2 CSV files" in captured.out
