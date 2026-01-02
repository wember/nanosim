"""
Tests for refactored sim_plot.py functions.

Tests the extracted testable functions directly.
"""

import os
import sys
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

# Add creutz-sim to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))


class TestFindMostRecentCsv:
    """Test find_most_recent_csv function."""

    def test_with_files(self, tmp_path):
        """Test finding most recent CSV file."""
        from sim_plot import find_most_recent_csv

        # Create mock data directory structure
        data_dir = tmp_path / "data" / "r1"
        data_dir.mkdir(parents=True)

        # Create multiple CSV files with different timestamps
        old_file = data_dir / "sim_data_1.csv"
        new_file = data_dir / "sim_data_2.csv"

        csv_data = "t,K,U,N0,Nx,S/nk,n\n1,0.5,0.5,0.1,0.2,1.0,100\n"
        old_file.write_text(csv_data)

        import time

        time.sleep(0.01)
        new_file.write_text(csv_data)

        # Test with custom pattern
        patterns = [str(data_dir / "sim_data*.csv")]
        result = find_most_recent_csv(str(tmp_path), patterns)

        assert result == str(new_file)

    def test_no_files(self, tmp_path):
        """Test when no CSV files are found."""
        from sim_plot import find_most_recent_csv

        result = find_most_recent_csv(str(tmp_path))
        assert result is None

    def test_default_patterns(self, tmp_path):
        """Test with default search patterns."""
        from sim_plot import find_most_recent_csv

        # Create file in default location
        data_dir = tmp_path / "data" / "r0"
        data_dir.mkdir(parents=True)

        csv_file = data_dir / "sim_data.csv"
        csv_file.write_text("t,K,U,N0,Nx,S/nk,n\n1,0.5,0.5,0.1,0.2,1.0,100\n")

        result = find_most_recent_csv(str(tmp_path))
        assert result == str(csv_file)

    def test_multiple_directories(self, tmp_path):
        """Test finding most recent across multiple directories."""
        from sim_plot import find_most_recent_csv

        # Create files in different directories
        r1_dir = tmp_path / "data" / "r1"
        r1_dir.mkdir(parents=True)
        r2_dir = tmp_path / "data" / "r2"
        r2_dir.mkdir(parents=True)

        old_file = r1_dir / "sim_data.csv"
        old_file.write_text("t,K,U,N0,Nx,S/nk,n\n1,0.5,0.5,0.1,0.2,1.0,100\n")

        import time

        time.sleep(0.01)

        new_file = r2_dir / "sim_data.csv"
        new_file.write_text("t,K,U,N0,Nx,S/nk,n\n1,0.5,0.5,0.1,0.2,1.0,100\n")

        result = find_most_recent_csv(str(tmp_path))
        assert result == str(new_file)


class TestValidateCsvColumns:
    """Test validate_csv_columns function."""

    def test_all_columns_present(self):
        """Test validation when all columns present."""
        from sim_plot import validate_csv_columns

        df = pd.DataFrame(
            {
                "t": [1, 2],
                "K": [0.5, 0.6],
                "U": [0.3, 0.4],
                "N0": [0.1, 0.2],
                "Nx": [0.2, 0.3],
                "n": [100, 100],
            }
        )

        is_valid, missing = validate_csv_columns(df)
        assert is_valid
        assert len(missing) == 0

    def test_missing_columns(self):
        """Test validation when columns missing."""
        from sim_plot import validate_csv_columns

        df = pd.DataFrame({"t": [1, 2], "K": [0.5, 0.6], "U": [0.3, 0.4]})

        is_valid, missing = validate_csv_columns(df)
        assert not is_valid
        assert "N0" in missing
        assert "Nx" in missing
        assert "n" in missing

    def test_custom_required_columns(self):
        """Test with custom required columns list."""
        from sim_plot import validate_csv_columns

        df = pd.DataFrame({"t": [1, 2], "K": [0.5, 0.6]})

        is_valid, missing = validate_csv_columns(df, required_cols=["t", "K"])
        assert is_valid
        assert len(missing) == 0

        is_valid, missing = validate_csv_columns(df, required_cols=["t", "K", "X"])
        assert not is_valid
        assert "X" in missing

    def test_extra_columns_ok(self):
        """Test that extra columns don't cause validation failure."""
        from sim_plot import validate_csv_columns

        df = pd.DataFrame(
            {
                "t": [1, 2],
                "K": [0.5, 0.6],
                "U": [0.3, 0.4],
                "N0": [0.1, 0.2],
                "Nx": [0.2, 0.3],
                "n": [100, 100],
                "extra": [1, 2],
            }
        )

        is_valid, missing = validate_csv_columns(df)
        assert is_valid


class TestLoadAndValidateCsv:
    """Test load_and_validate_csv function."""

    def test_load_valid_csv(self, tmp_path):
        """Test loading valid CSV file."""
        from sim_plot import load_and_validate_csv

        csv_file = tmp_path / "valid.csv"
        csv_file.write_text(
            "t,K,U,N0,Nx,S/nk,n\n1,0.5,0.5,0.1,0.2,1.0,100\n2,0.6,0.6,0.2,0.3,1.1,100\n"
        )

        df = load_and_validate_csv(str(csv_file))
        assert len(df) == 2
        assert df["n"].iloc[0] == 100
        assert df["t"].iloc[0] == 1
        assert df["t"].iloc[1] == 2

    def test_load_missing_columns(self, tmp_path):
        """Test loading CSV with missing columns raises ValueError."""
        from sim_plot import load_and_validate_csv

        csv_file = tmp_path / "invalid.csv"
        csv_file.write_text("t,K,U\n1,0.5,0.5\n")

        with pytest.raises(ValueError, match="missing required columns"):
            load_and_validate_csv(str(csv_file))

    def test_file_not_found(self):
        """Test loading nonexistent file raises FileNotFoundError."""
        from sim_plot import load_and_validate_csv

        with pytest.raises(FileNotFoundError):
            load_and_validate_csv("/nonexistent/file.csv")

    def test_error_message_includes_missing_columns(self, tmp_path):
        """Test error message includes list of missing columns."""
        from sim_plot import load_and_validate_csv

        csv_file = tmp_path / "partial.csv"
        csv_file.write_text("t,K\n1,0.5\n")

        with pytest.raises(ValueError) as exc_info:
            load_and_validate_csv(str(csv_file))

        assert "N0" in str(exc_info.value)
        assert "Nx" in str(exc_info.value)


class TestCalculateEntropy:
    """Test calculate_entropy function."""

    def test_basic_calculation(self):
        """Test entropy calculation with valid data."""
        from sim_plot import calculate_entropy

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

        entropy = calculate_entropy(df)

        assert len(entropy) == 3
        assert all(np.isfinite(entropy))
        assert all(entropy > 0)

    def test_n0_zero_edge_case(self):
        """Test entropy calculation with N0=0 edge case."""
        from sim_plot import calculate_entropy

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

        entropy = calculate_entropy(df)

        # Should handle N0=0 without errors
        assert np.isfinite(entropy.iloc[0])
        assert np.isfinite(entropy.iloc[1])

    def test_entropy_increases_with_demon_energy(self):
        """Test that entropy generally increases with demon energy."""
        from sim_plot import calculate_entropy

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

        entropy = calculate_entropy(df)

        # Higher demon energy should generally give higher entropy
        assert entropy.iloc[1] > entropy.iloc[0]
        assert entropy.iloc[2] > entropy.iloc[1]

    def test_entropy_with_large_values(self):
        """Test entropy calculation doesn't overflow with large values."""
        from sim_plot import calculate_entropy

        df = pd.DataFrame(
            {
                "t": [1],
                "K": [1000.0],
                "U": [0.5],
                "N0": [100.0],
                "Nx": [200.0],
                "n": [10000],
            }
        )

        entropy = calculate_entropy(df)

        assert np.isfinite(entropy.iloc[0])

    def test_entropy_with_varying_n(self):
        """Test entropy calculation with different lattice sizes."""
        from sim_plot import calculate_entropy

        # Smaller lattice
        df_small = pd.DataFrame(
            {
                "t": [1],
                "K": [50.0],
                "U": [0.3],
                "N0": [10.0],
                "Nx": [20.0],
                "n": [100],
            }
        )

        # Larger lattice
        df_large = pd.DataFrame(
            {
                "t": [1],
                "K": [500.0],
                "U": [0.3],
                "N0": [100.0],
                "Nx": [200.0],
                "n": [1000],
            }
        )

        entropy_small = calculate_entropy(df_small)
        entropy_large = calculate_entropy(df_large)

        # Both should be finite and positive
        assert np.isfinite(entropy_small.iloc[0])
        assert np.isfinite(entropy_large.iloc[0])
        assert entropy_small.iloc[0] > 0
        assert entropy_large.iloc[0] > 0


class TestCreateSimulationPlot:
    """Test create_simulation_plot function."""

    def test_creates_figure(self):
        """Test that function creates a figure object."""
        from sim_plot import create_simulation_plot

        df = pd.DataFrame(
            {
                "t": list(range(1, 11)),
                "K": [0.5 + i * 0.01 for i in range(10)],
                "U": [0.3 + i * 0.01 for i in range(10)],
                "N0": [0.1] * 10,
                "Nx": [0.2] * 10,
                "S/nk": [1.0 + i * 0.01 for i in range(10)],
                "n": [100] * 10,
            }
        )

        fig = create_simulation_plot(df)

        assert fig is not None
        assert hasattr(fig, "data")
        assert hasattr(fig, "layout")

    def test_has_three_traces(self):
        """Test that plot has 3 traces (3 subplots)."""
        from sim_plot import create_simulation_plot

        df = pd.DataFrame(
            {
                "t": [1, 2, 3],
                "K": [0.5, 0.6, 0.7],
                "U": [0.3, 0.4, 0.5],
                "N0": [0.1, 0.2, 0.3],
                "Nx": [0.2, 0.3, 0.4],
                "n": [100, 100, 100],
            }
        )

        fig = create_simulation_plot(df)

        assert len(fig.data) == 3

    def test_default_title(self):
        """Test default title includes lattice size."""
        from sim_plot import create_simulation_plot

        df = pd.DataFrame(
            {
                "t": [1, 2],
                "K": [0.5, 0.6],
                "U": [0.3, 0.4],
                "N0": [0.1, 0.2],
                "Nx": [0.2, 0.3],
                "n": [100, 100],
            }
        )

        fig = create_simulation_plot(df)

        assert "Single Simulation Results" in fig.layout.title.text
        assert "100" in fig.layout.title.text

    def test_custom_title(self):
        """Test creating plot with custom title."""
        from sim_plot import create_simulation_plot

        df = pd.DataFrame(
            {
                "t": [1, 2],
                "K": [0.5, 0.6],
                "U": [0.3, 0.4],
                "N0": [0.1, 0.2],
                "Nx": [0.2, 0.3],
                "n": [100, 100],
            }
        )

        custom_title = "Test Simulation Results"
        fig = create_simulation_plot(df, title=custom_title)

        assert fig.layout.title.text == custom_title

    def test_plot_with_minimal_data(self):
        """Test plot creation with minimal dataset."""
        from sim_plot import create_simulation_plot

        df = pd.DataFrame(
            {
                "t": [1],
                "K": [0.5],
                "U": [0.3],
                "N0": [0.1],
                "Nx": [0.2],
                "n": [100],
            }
        )

        fig = create_simulation_plot(df)

        assert fig is not None
        assert len(fig.data) == 3

    def test_plot_with_long_timeseries(self):
        """Test plot creation with long time series."""
        from sim_plot import create_simulation_plot

        n_points = 1000
        df = pd.DataFrame(
            {
                "t": list(range(1, n_points + 1)),
                "K": [0.5 + i * 0.0001 for i in range(n_points)],
                "U": [0.3 + i * 0.0001 for i in range(n_points)],
                "N0": [0.1] * n_points,
                "Nx": [0.2] * n_points,
                "n": [100] * n_points,
            }
        )

        fig = create_simulation_plot(df)

        assert fig is not None
        # Verify each trace has correct number of points
        assert len(fig.data[0].x) == n_points


class TestMainFunction:
    """Test the main() function."""

    @patch("sim_plot.sys.argv", ["sim_plot.py", "test_data.csv"])
    @patch("sim_plot.load_and_validate_csv")
    @patch("sim_plot.create_simulation_plot")
    @patch("plotly.graph_objects.Figure.show")
    def test_with_file_argument(self, mock_show, mock_create_plot, mock_load_csv):
        """Test main function with file argument."""
        from sim_plot import main

        # Mock return values
        mock_df = pd.DataFrame(
            {"t": [1], "K": [0.5], "U": [0.3], "N0": [0.1], "Nx": [0.2], "n": [100]}
        )
        mock_load_csv.return_value = mock_df

        import plotly.graph_objects as go

        mock_fig = go.Figure()
        mock_create_plot.return_value = mock_fig

        # Run main (should not raise)
        main()

        # Verify functions were called
        mock_load_csv.assert_called_once_with("test_data.csv")
        mock_create_plot.assert_called_once()
        mock_show.assert_called_once()

    @patch("sim_plot.sys.argv", ["sim_plot.py"])
    @patch("sim_plot.find_most_recent_csv")
    def test_no_files_found(self, mock_find_csv):
        """Test main function when no files found."""
        from sim_plot import main

        mock_find_csv.return_value = None

        with pytest.raises(SystemExit) as exc_info:
            main()

        assert exc_info.value.code == 1

    @patch("sim_plot.sys.argv", ["sim_plot.py", "nonexistent.csv"])
    @patch("sim_plot.load_and_validate_csv")
    def test_file_not_found_error(self, mock_load_csv):
        """Test main function when file doesn't exist."""
        from sim_plot import main

        mock_load_csv.side_effect = FileNotFoundError()

        with pytest.raises(SystemExit) as exc_info:
            main()

        assert exc_info.value.code == 1

    @patch("sim_plot.sys.argv", ["sim_plot.py", "invalid.csv"])
    @patch("sim_plot.load_and_validate_csv")
    def test_validation_error(self, mock_load_csv):
        """Test main function when CSV validation fails."""
        from sim_plot import main

        mock_load_csv.side_effect = ValueError("missing required columns")

        with pytest.raises(SystemExit) as exc_info:
            main()

        assert exc_info.value.code == 1

    @patch("sim_plot.sys.argv", ["sim_plot.py"])
    @patch("sim_plot.find_most_recent_csv")
    @patch("sim_plot.load_and_validate_csv")
    @patch("sim_plot.create_simulation_plot")
    @patch("plotly.graph_objects.Figure.show")
    def test_auto_discovery_mode(
        self, mock_show, mock_create_plot, mock_load_csv, mock_find_csv
    ):
        """Test main function in auto-discovery mode."""
        from sim_plot import main

        # Mock file discovery
        mock_find_csv.return_value = "/path/to/recent.csv"

        # Mock data loading
        mock_df = pd.DataFrame(
            {"t": [1], "K": [0.5], "U": [0.3], "N0": [0.1], "Nx": [0.2], "n": [100]}
        )
        mock_load_csv.return_value = mock_df

        import plotly.graph_objects as go

        mock_fig = go.Figure()
        mock_create_plot.return_value = mock_fig

        # Run main
        main()

        # Verify auto-discovery was used
        mock_find_csv.assert_called_once()
        mock_load_csv.assert_called_once_with("/path/to/recent.csv")
        mock_show.assert_called_once()


class TestMainExecution:
    """Test __main__ block execution."""

    def test_main_function_is_callable(self, tmp_path, monkeypatch):
        """Test that main() function exists and is callable."""
        from sim_plot import main

        # Create a fake CSV file
        csv_file = tmp_path / "test.csv"
        csv_file.write_text("t,K,U,N0,Nx,S/nk,n\n1,0.5,0.5,0.1,0.2,1.0,100\n")

        # Mock sys.argv to provide file path
        monkeypatch.setattr("sys.argv", ["sim_plot.py", str(csv_file)])

        # Mock fig.show() to avoid displaying window
        with patch("sim_plot.create_simulation_plot") as mock_create:
            from unittest.mock import MagicMock

            mock_fig = MagicMock()
            mock_create.return_value = mock_fig

            # Call main - this tests the function itself
            main()

            # Verify plot was created and shown
            mock_create.assert_called_once()
            mock_fig.show.assert_called_once()
