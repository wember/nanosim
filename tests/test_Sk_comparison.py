"""
Tests for refactored Sk_comparison.py functions.

Tests the extracted testable functions directly.
"""

import os
import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import plotly.graph_objects as go

# Add creutz-sim to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))


class TestGetColorPalette:
    """Test get_color_palette function."""

    def test_returns_11_colors(self):
        """Test that function returns 11 colors."""
        from Sk_comparison import get_color_palette

        colors = get_color_palette()
        assert len(colors) == 11

    def test_returns_hex_strings(self):
        """Test that all colors are hex strings."""
        from Sk_comparison import get_color_palette

        colors = get_color_palette()
        for color in colors:
            assert isinstance(color, str)
            assert color.startswith("#")
            assert len(color) == 7  # #RRGGBB format

    def test_colors_are_distinct(self):
        """Test that all colors are unique."""
        from Sk_comparison import get_color_palette

        colors = get_color_palette()
        assert len(set(colors)) == 11


class TestLoadAndAverageCsvFiles:
    """Test load_and_average_csv_files function."""

    def test_single_csv_file(self, tmp_path):
        """Test loading single CSV file."""
        from Sk_comparison import load_and_average_csv_files

        csv_file = tmp_path / "sim_data.csv"
        csv_file.write_text(
            "t,K,U,N0,Nx,S/nk,n\n"
            "1,0.5,0.3,0.1,0.2,1.0,100\n"
            "2,0.6,0.4,0.2,0.3,1.1,100\n"
        )

        df = load_and_average_csv_files(str(tmp_path))

        assert len(df) == 2
        assert "t" in df.columns
        assert "S/nk" in df.columns

    def test_multiple_csv_files_averaged(self, tmp_path):
        """Test that multiple CSV files are averaged correctly."""
        from Sk_comparison import load_and_average_csv_files

        # Create two CSV files with different values
        csv1 = tmp_path / "sim_data_1.csv"
        csv1.write_text(
            "t,K,U,N0,Nx,S/nk,n\n"
            "1,0.5,0.3,0.1,0.2,1.0,100\n"
            "2,0.6,0.4,0.2,0.3,1.2,100\n"
        )

        csv2 = tmp_path / "sim_data_2.csv"
        csv2.write_text(
            "t,K,U,N0,Nx,S/nk,n\n"
            "1,0.7,0.5,0.3,0.4,2.0,100\n"
            "2,0.8,0.6,0.4,0.5,2.2,100\n"
        )

        df = load_and_average_csv_files(str(tmp_path))

        # Check averaged values (should be mean of the two files)
        assert np.isclose(df["S/nk"].iloc[0], 1.5)  # (1.0 + 2.0) / 2
        assert np.isclose(df["S/nk"].iloc[1], 1.7)  # (1.2 + 2.2) / 2

    def test_empty_folder(self, tmp_path):
        """Test handling of empty folder."""
        from Sk_comparison import load_and_average_csv_files

        df = load_and_average_csv_files(str(tmp_path))

        assert df.empty

    def test_nonexistent_folder(self):
        """Test handling of nonexistent folder."""
        from Sk_comparison import load_and_average_csv_files

        df = load_and_average_csv_files("/nonexistent/path")

        assert df.empty


class TestCalculateZoomIndices:
    """Test calculate_zoom_indices function."""

    def test_basic_calculation(self):
        """Test zoom indices for typical length."""
        from Sk_comparison import calculate_zoom_indices

        start, end = calculate_zoom_indices(100)

        # Should be middle 25 elements (25% of 100)
        assert end - start == 25
        # Should be centered around index 50
        assert 37 <= start <= 38
        assert 62 <= end <= 63

    def test_small_length(self):
        """Test zoom indices for small length."""
        from Sk_comparison import calculate_zoom_indices

        start, end = calculate_zoom_indices(20)

        # Should be middle 5 elements (25% of 20)
        assert end - start == 5

    def test_large_length(self):
        """Test zoom indices for large length."""
        from Sk_comparison import calculate_zoom_indices

        start, end = calculate_zoom_indices(1000)

        # Should be middle 250 elements (25% of 1000)
        assert end - start == 250

    def test_indices_are_integers(self):
        """Test that returned indices are integers."""
        from Sk_comparison import calculate_zoom_indices

        start, end = calculate_zoom_indices(100)

        assert isinstance(start, int)
        assert isinstance(end, int)


class TestAddEntropyTraces:
    """Test add_entropy_traces function."""

    def test_adds_four_traces(self):
        """Test that function adds 4 traces (one per subplot)."""
        from plotly.subplots import make_subplots
        from Sk_comparison import add_entropy_traces

        fig = make_subplots(rows=2, cols=2)
        df = pd.DataFrame(
            {
                "t": list(range(1, 101)),
                "S/nk": [1.0 + i * 0.01 for i in range(100)],
            }
        )

        initial_traces = len(fig.data)
        add_entropy_traces(fig, df, radius=0, color="#800020", label_prefix="irr ")

        # Should add 4 traces
        assert len(fig.data) == initial_traces + 4

    def test_trace_names(self):
        """Test that trace names include radius."""
        from plotly.subplots import make_subplots
        from Sk_comparison import add_entropy_traces

        fig = make_subplots(rows=2, cols=2)
        df = pd.DataFrame(
            {
                "t": list(range(1, 101)),
                "S/nk": [1.0] * 100,
            }
        )

        add_entropy_traces(fig, df, radius=5, color="#800020", label_prefix="irr ")

        # Check that traces reference radius
        assert any("5" in trace.name for trace in fig.data)

    def test_returns_zoom_indices(self):
        """Test that function returns zoom indices."""
        from plotly.subplots import make_subplots
        from Sk_comparison import add_entropy_traces

        fig = make_subplots(rows=2, cols=2)
        df = pd.DataFrame(
            {
                "t": list(range(1, 101)),
                "S/nk": [1.0] * 100,
            }
        )

        start, end = add_entropy_traces(
            fig, df, radius=0, color="#800020", label_prefix=""
        )

        assert isinstance(start, int)
        assert isinstance(end, int)
        assert end > start

    def test_color_applied(self):
        """Test that specified color is applied to traces."""
        from plotly.subplots import make_subplots
        from Sk_comparison import add_entropy_traces

        fig = make_subplots(rows=2, cols=2)
        df = pd.DataFrame(
            {
                "t": list(range(1, 101)),
                "S/nk": [1.0] * 100,
            }
        )

        test_color = "#FF0000"
        add_entropy_traces(fig, df, radius=0, color=test_color, label_prefix="")

        # At least one trace should have the specified color
        assert any(
            trace.line.color == test_color
            for trace in fig.data
            if hasattr(trace, "line")
        )


class TestCreateEntropyComparisonPlot:
    """Test create_entropy_comparison_plot function."""

    def test_creates_figure(self, tmp_path):
        """Test that function creates a figure."""
        from Sk_comparison import create_entropy_comparison_plot

        # Create minimal test data
        for sim_type in ["", "irr/"]:
            data_dir = tmp_path / "data" / sim_type / "r0"
            data_dir.mkdir(parents=True)
            csv_file = data_dir / "sim_data.csv"
            csv_file.write_text(
                "t,K,U,N0,Nx,S/nk,n\n"
                + "\n".join(
                    [f"{i},0.5,0.3,0.1,0.2,{1.0 + i*0.01},100" for i in range(1, 101)]
                )
            )

        fig = create_entropy_comparison_plot(str(tmp_path), max_radius=1)

        assert fig is not None
        assert hasattr(fig, "data")
        assert hasattr(fig, "layout")

    def test_handles_empty_directories(self, tmp_path):
        """Test handling when no data directories exist."""
        from Sk_comparison import create_entropy_comparison_plot

        fig = create_entropy_comparison_plot(str(tmp_path), max_radius=1)

        assert fig is not None
        # Should create empty figure without errors

    def test_title_includes_lattice_size(self, tmp_path):
        """Test that title includes lattice size when data exists."""
        from Sk_comparison import create_entropy_comparison_plot

        data_dir = tmp_path / "data" / "r0"
        data_dir.mkdir(parents=True)
        csv_file = data_dir / "sim_data.csv"
        csv_file.write_text(
            "t,K,U,N0,Nx,S/nk,n\n"
            + "\n".join([f"{i},0.5,0.3,0.1,0.2,1.0,500" for i in range(1, 11)])
        )

        fig = create_entropy_comparison_plot(str(tmp_path), max_radius=1)

        assert "500" in fig.layout.title.text

    def test_custom_bin_size(self, tmp_path):
        """Test using custom bin_size parameter."""
        from Sk_comparison import create_entropy_comparison_plot

        data_dir = tmp_path / "data" / "r0"
        data_dir.mkdir(parents=True)
        csv_file = data_dir / "sim_data.csv"
        csv_file.write_text(
            "t,K,U,N0,Nx,S/nk,n\n"
            + "\n".join([f"{i},0.5,0.3,0.1,0.2,1.0,100" for i in range(1, 101)])
        )

        # Should not raise with custom bin_size
        fig = create_entropy_comparison_plot(str(tmp_path), max_radius=1, bin_size=20)

        assert fig is not None

    def test_processes_multiple_radii(self, tmp_path):
        """Test processing multiple radius directories."""
        from Sk_comparison import create_entropy_comparison_plot

        # Create data for R=0, R=1, R=2
        for R in range(3):
            data_dir = tmp_path / "data" / f"r{R}"
            data_dir.mkdir(parents=True)
            csv_file = data_dir / "sim_data.csv"
            csv_file.write_text(
                "t,K,U,N0,Nx,S/nk,n\n"
                + "\n".join([f"{i},0.5,0.3,0.1,0.2,1.0,100" for i in range(1, 101)])
            )

        fig = create_entropy_comparison_plot(str(tmp_path), max_radius=3)

        # Should have traces from all 3 radii (4 traces per radius)
        assert len(fig.data) >= 12

    def test_adds_midpoint_lines(self, tmp_path):
        """Test that midpoint vertical lines are added."""
        from Sk_comparison import create_entropy_comparison_plot

        data_dir = tmp_path / "data" / "r0"
        data_dir.mkdir(parents=True)
        csv_file = data_dir / "sim_data.csv"
        csv_file.write_text(
            "t,K,U,N0,Nx,S/nk,n\n"
            + "\n".join([f"{i},0.5,0.3,0.1,0.2,1.0,100" for i in range(1, 101)])
        )

        fig = create_entropy_comparison_plot(str(tmp_path), max_radius=1)

        # Should have vlines in layout shapes
        assert len(fig.layout.shapes) >= 2


class TestMainFunction:
    """Test the main() function."""

    @patch("Sk_comparison.create_entropy_comparison_plot")
    @patch("plotly.graph_objects.Figure.show")
    def test_main_execution(self, mock_show, mock_create_plot):
        """Test main function execution flow."""
        from Sk_comparison import main

        mock_fig = go.Figure()
        mock_create_plot.return_value = mock_fig

        # Should not raise
        main()

        # Verify calls
        mock_create_plot.assert_called_once()
        mock_show.assert_called_once()


class TestMainExecution:
    """Test __main__ block execution."""

    def test_main_function_is_callable(self):
        """Test that main() function exists and is callable."""
        from Sk_comparison import main

        # Mock fig.show() to avoid displaying window
        with patch("Sk_comparison.create_entropy_comparison_plot") as mock_create:
            mock_fig = MagicMock()
            mock_create.return_value = mock_fig

            # Call main - this tests the function itself
            main()

            # Verify it was called
            mock_create.assert_called_once()
            mock_fig.show.assert_called_once()
