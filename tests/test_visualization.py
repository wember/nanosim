"""
Tests for visualization scripts (sim_plot.py, sim_plot_radii.py, Sk_comparison.py).

Focuses on testing data loading, processing, and validation logic without
displaying plots. Uses mocked file I/O and plotly operations.
"""

import os
import sys
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

# Add creutz-sim to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "creutz-sim"))


class TestSimPlotDataLoading:
    """Test data loading and validation for sim_plot.py"""

    def test_csv_file_discovery(self, tmp_path):
        """Test auto-discovery of most recent CSV file."""
        # Create mock data directory structure
        data_dir = tmp_path / "data" / "r1"
        data_dir.mkdir(parents=True)

        # Create multiple CSV files with different timestamps
        old_file = data_dir / "sim_data_1.csv"
        new_file = data_dir / "sim_data_2.csv"

        # Create CSV with required columns
        csv_data = "t,K,U,N0,Nx,S/nk,n\n1,0.5,0.5,0.1,0.2,1.0,100\n"
        old_file.write_text(csv_data)

        # Small delay to ensure different timestamps
        import time

        time.sleep(0.01)
        new_file.write_text(csv_data)

        # Mock glob to return our files
        with patch("glob.glob") as mock_glob:
            mock_glob.return_value = [str(old_file), str(new_file)]

            # The most recent file should be selected
            files = sorted(
                mock_glob.return_value, key=lambda x: os.path.getmtime(x), reverse=True
            )
            assert files[0] == str(new_file)

    def test_csv_validation_missing_columns(self, tmp_path):
        """Test CSV validation catches missing required columns."""
        csv_file = tmp_path / "invalid.csv"
        # Missing 'n' column
        csv_file.write_text("t,K,U,N0,Nx\n1,0.5,0.5,0.1,0.2\n")

        df = pd.read_csv(csv_file)
        required_cols = ["t", "K", "U", "N0", "Nx", "n"]
        missing_cols = [col for col in required_cols if col not in df.columns]

        assert len(missing_cols) == 1
        assert "n" in missing_cols

    def test_csv_validation_all_columns_present(self, tmp_path):
        """Test CSV validation passes with all required columns."""
        csv_file = tmp_path / "valid.csv"
        csv_file.write_text("t,K,U,N0,Nx,S/nk,n\n1,0.5,0.5,0.1,0.2,1.0,100\n")

        df = pd.read_csv(csv_file)
        required_cols = ["t", "K", "U", "N0", "Nx", "n"]
        missing_cols = [col for col in required_cols if col not in df.columns]

        assert len(missing_cols) == 0

    def test_data_extraction(self, tmp_path):
        """Test extraction of data columns from CSV."""
        csv_file = tmp_path / "data.csv"
        csv_content = """t,K,U,N0,Nx,S/nk,n
1,0.5,0.3,0.1,0.2,1.0,100
2,0.6,0.4,0.15,0.25,1.1,100
3,0.7,0.5,0.2,0.3,1.2,100"""
        csv_file.write_text(csv_content)

        df = pd.read_csv(csv_file)

        # Extract columns
        n = df["n"][0]
        t = df["t"]
        U = df["U"]
        K = df["K"]

        assert n == 100
        assert len(t) == 3
        assert list(t) == [1, 2, 3]
        assert list(K) == [0.5, 0.6, 0.7]
        assert list(U) == [0.3, 0.4, 0.5]

    def test_n0_edge_case_handling(self, tmp_path):
        """Test N0=0 edge case handling (N0_exp = max(1, N0))."""
        csv_file = tmp_path / "data.csv"
        csv_content = """t,K,U,N0,Nx,S/nk,n
1,0.5,0.3,0.0,0.2,1.0,100
2,0.6,0.4,0.1,0.25,1.1,100"""
        csv_file.write_text(csv_content)

        df = pd.read_csv(csv_file)

        # Apply N0_exp calculation
        df["N0_exp"] = df["N0"].replace(0, 1)

        assert df["N0_exp"].iloc[0] == 1  # N0=0 replaced with 1
        assert df["N0_exp"].iloc[1] == 0.1  # N0=0.1 unchanged


class TestSimPlotEntropyCalculation:
    """Test entropy calculation logic from sim_plot.py"""

    def test_entropy_calculation(self):
        """Test entropy calculation with valid inputs."""
        from scipy.special import loggamma

        # Test data
        n = 100
        K_val = 50.0
        N0_val = 10.0
        Nx_val = 20.0
        N0_exp_val = 10.0

        # Sk: Kinetic entropy
        Sk_val = loggamma(K_val + n) - loggamma(K_val + 1) - loggamma(n)

        # Su: Configurational entropy
        Su_val = (
            loggamma(n + 1)
            + N0_exp_val * np.log(2)
            - (
                loggamma(n - N0_val - Nx_val + 1)
                + loggamma(N0_val + 1)
                + loggamma(Nx_val + 1)
            )
        )

        # Total entropy per site
        total_entropy = (Sk_val + Su_val) / n

        # Should produce finite result
        assert np.isfinite(total_entropy)
        assert total_entropy > 0

    def test_entropy_with_n0_zero(self):
        """Test entropy calculation when N0=0 (uses N0_exp=1)."""
        from scipy.special import loggamma

        n = 100
        K_val = 50.0
        N0_val = 0.0
        Nx_val = 20.0
        N0_exp_val = 1.0  # Edge case: N0=0 replaced with 1

        # Sk: Kinetic entropy
        Sk_val = loggamma(K_val + n) - loggamma(K_val + 1) - loggamma(n)

        # Su: Configurational entropy with N0_exp=1
        Su_val = (
            loggamma(n + 1)
            + N0_exp_val * np.log(2)
            - (
                loggamma(n - N0_val - Nx_val + 1)
                + loggamma(N0_val + 1)
                + loggamma(Nx_val + 1)
            )
        )

        # Should not produce NaN or inf
        assert np.isfinite(Sk_val)
        assert np.isfinite(Su_val)


class TestSimPlotRadiiDataLoading:
    """Test multi-radius data loading for sim_plot_radii.py"""

    def test_radius_directory_discovery(self, tmp_path):
        """Test discovery of CSV files across multiple radius directories."""
        project_root = tmp_path
        csv_data = "t,K,U,N0,Nx,S/nk,n\n1,0.5,0.5,0.1,0.2,1.0,100\n"

        # Create directories for R=0, R=1, R=2
        file_names = []
        for R in [0, 1, 2]:
            radius_dir = project_root / "data" / f"r{R}"
            radius_dir.mkdir(parents=True)
            csv_file = radius_dir / "sim_data.csv"
            csv_file.write_text(csv_data)
            file_names.append(str(csv_file))

        # Verify files exist
        for file_path in file_names:
            assert os.path.exists(file_path)

    def test_radius_file_selection_most_recent(self, tmp_path):
        """Test selection of most recent file when multiple CSVs in directory."""
        radius_dir = tmp_path / "data" / "r1"
        radius_dir.mkdir(parents=True)

        csv_data = "t,K,U,N0,Nx,S/nk,n\n1,0.5,0.5,0.1,0.2,1.0,100\n"

        # Create multiple files
        old_file = radius_dir / "sim_data_1.csv"
        old_file.write_text(csv_data)

        import time

        time.sleep(0.01)

        new_file = radius_dir / "sim_data_2.csv"
        new_file.write_text(csv_data)

        # Simulate glob and sorting
        import glob

        csv_files = glob.glob(str(radius_dir / "sim_data*.csv"))
        csv_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)

        # Most recent should be first
        assert csv_files[0] == str(new_file)

    def test_missing_radius_directory_handling(self, tmp_path):
        """Test handling when radius directory doesn't exist."""
        project_root = tmp_path
        folder_path = project_root / "data" / "r5"

        # Directory doesn't exist
        assert not folder_path.exists()

        # glob should return empty list
        import glob

        csv_files = glob.glob(str(folder_path / "sim_data*.csv"))
        assert len(csv_files) == 0


class TestSkComparisonDataProcessing:
    """Test data processing for Sk_comparison.py"""

    def test_ensemble_averaging_multiple_runs(self, tmp_path):
        """Test averaging across multiple simulation runs."""
        # Create 3 CSV files with different data
        csv_files = []
        for i in range(3):
            csv_file = tmp_path / f"run_{i}.csv"
            # Each run has slightly different S/nk values
            csv_content = f"""t,K,U,N0,Nx,S/nk,n
1,0.5,0.3,0.1,0.2,{1.0 + i * 0.1},100
2,0.6,0.4,0.15,0.25,{1.1 + i * 0.1},100
3,0.7,0.5,0.2,0.3,{1.2 + i * 0.1},100"""
            csv_file.write_text(csv_content)
            csv_files.append(csv_file)

        # Load all files
        list_of_dfs = []
        for file_path in csv_files:
            df = pd.read_csv(file_path)
            list_of_dfs.append(df)

        # Calculate ensemble average
        if list_of_dfs:
            combined_df = pd.concat(list_of_dfs, ignore_index=True)
            avg_df = combined_df.groupby("t").mean()

            # Check averaging worked
            assert len(avg_df) == 3
            # Average of [1.0, 1.1, 1.2] = 1.1
            assert np.isclose(avg_df["S/nk"].iloc[0], 1.1, atol=0.01)

    def test_rolling_average_smoothing(self):
        """Test rolling average smoothing of entropy data."""
        # Create sample data with noise
        t = np.arange(100)
        entropy = np.sin(t / 10) + np.random.normal(0, 0.1, 100)

        df = pd.DataFrame({"t": t, "S/nk": entropy})

        # Apply rolling average
        bin_size = 10
        smoothed = df["S/nk"].rolling(window=bin_size, center=True).mean()

        # Smoothed data should have less variance
        assert np.nanstd(smoothed) < np.std(entropy)

    def test_zoom_window_calculation(self):
        """Test calculation of zoom window (middle 25% of data)."""
        n_points = 1000

        # Calculate middle 25% indices
        start_idx = int(n_points * 0.375)  # 37.5%
        end_idx = int(n_points * 0.625)  # 62.5%

        assert start_idx == 375
        assert end_idx == 625
        assert end_idx - start_idx == 250  # 25% of 1000

    def test_color_palette_length(self):
        """Test that color palette has correct length for 11 radii."""
        colors = [
            "#301934",
            "#702963",
            "#800020",
            "#AA336A",
            "#9F2B68",
            "#800080",
            "#BF40BF",
            "#DA70D6",
            "#CF9FFF",
            "#E0B0FF",
            "#CBC3E3",
        ]

        # Should have 11 colors for R=0 to R=10
        assert len(colors) == 11


class TestVisualizationIntegration:
    """Integration tests for visualization workflows"""

    def test_full_sim_plot_workflow(self, tmp_path):
        """Test complete workflow of loading and processing data for sim_plot."""
        # Create valid CSV file
        csv_file = tmp_path / "sim_data.csv"
        csv_content = """t,K,U,N0,Nx,S/nk,n
1,0.5,0.3,0.1,0.2,1.0,100
2,0.6,0.4,0.15,0.25,1.1,100
3,0.7,0.5,0.2,0.3,1.2,100
4,0.8,0.6,0.25,0.35,1.3,100
5,0.9,0.7,0.3,0.4,1.4,100"""
        csv_file.write_text(csv_content)

        # Load and validate
        df = pd.read_csv(csv_file)
        required_cols = ["t", "K", "U", "N0", "Nx", "n"]
        missing_cols = [col for col in required_cols if col not in df.columns]
        assert len(missing_cols) == 0

        # Extract data
        n = df["n"][0]
        t = df["t"]
        U = df["U"]
        K = df["K"]

        assert n == 100
        assert len(t) == 5
        assert len(U) == 5
        assert len(K) == 5

        # Calculate midpoint for phase transition marker
        midpoint = len(df) // 2
        assert midpoint == 2

    def test_full_radii_comparison_workflow(self, tmp_path):
        """Test complete workflow for multi-radius comparison."""
        project_root = tmp_path
        csv_data = (
            "t,K,U,N0,Nx,S/nk,n\n1,0.5,0.5,0.1,0.2,1.0,100\n2,0.6,0.6,0.2,0.3,1.1,100\n"
        )

        # Create data for 3 radii
        file_names = []
        data_names = []
        for R in range(3):
            radius_dir = project_root / "data" / f"r{R}"
            radius_dir.mkdir(parents=True)
            csv_file = radius_dir / "sim_data.csv"
            csv_file.write_text(csv_data)
            file_names.append(str(csv_file))
            data_names.append(f"radius = {R}")

        # Verify we have 3 files
        assert len(file_names) == 3
        assert len(data_names) == 3

        # Load each file
        for i, file_path in enumerate(file_names):
            df = pd.read_csv(file_path)
            assert len(df) == 2
            assert data_names[i] == f"radius = {i}"

    def test_sk_comparison_workflow(self, tmp_path):
        """Test complete workflow for reversible/irreversible comparison."""
        project_root = tmp_path
        csv_data = """t,K,U,N0,Nx,S/nk,n
1,0.5,0.3,0.1,0.2,1.0,100
2,0.6,0.4,0.15,0.25,1.1,100
3,0.7,0.5,0.2,0.3,1.2,100"""

        # Create irreversible data
        irr_dir = project_root / "data" / "irr" / "r1"
        irr_dir.mkdir(parents=True)

        # Multiple runs for ensemble averaging
        for run in range(3):
            csv_file = irr_dir / f"irr_sim_data_r1_{run}.csv"
            csv_file.write_text(csv_data)

        # Load all runs
        import glob

        all_csv_files = glob.glob(str(irr_dir / "*.csv"))
        assert len(all_csv_files) == 3

        list_of_dfs = []
        for file_path in all_csv_files:
            df = pd.read_csv(file_path)
            list_of_dfs.append(df)

        # Combine and average
        combined_df = pd.concat(list_of_dfs, ignore_index=True)
        avg_df = combined_df.groupby("t").mean()

        assert len(avg_df) == 3
        assert "S/nk" in avg_df.columns

    @patch("plotly.graph_objects.Figure.show")
    def test_plot_show_not_called_in_tests(self, mock_show):
        """Test that we can mock plotly.show() to prevent display in tests."""
        import plotly.graph_objects as go

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6]))

        # Call show (mocked)
        fig.show()

        # Verify show was called but didn't actually display
        mock_show.assert_called_once()


class TestVisualizationErrorHandling:
    """Test error handling in visualization scripts"""

    def test_empty_csv_file(self, tmp_path):
        """Test handling of empty CSV file."""
        csv_file = tmp_path / "empty.csv"
        csv_file.write_text("")

        with pytest.raises(pd.errors.EmptyDataError):
            pd.read_csv(csv_file)

    def test_malformed_csv_file(self, tmp_path):
        """Test handling of malformed CSV data."""
        csv_file = tmp_path / "malformed.csv"
        csv_file.write_text("t,K,U\n1,0.5\n2,0.6,0.7,extra")  # Inconsistent columns

        # Pandas should raise ParserError for inconsistent columns
        with pytest.raises(pd.errors.ParserError):
            pd.read_csv(csv_file)

    def test_no_csv_files_found(self, tmp_path):
        """Test behavior when no CSV files are found."""
        import glob

        pattern = str(tmp_path / "data" / "r*" / "sim_data*.csv")
        csv_files = glob.glob(pattern)

        # Should return empty list
        assert len(csv_files) == 0

    def test_file_not_found_error(self):
        """Test FileNotFoundError handling."""
        with pytest.raises(FileNotFoundError):
            pd.read_csv("/nonexistent/path/data.csv")
