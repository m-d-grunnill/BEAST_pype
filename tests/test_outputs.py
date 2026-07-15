"""Unit tests for beast_pype.outputs module."""
import os
import pytest
import tempfile
import shutil
import numpy as np
import pandas as pd
from unittest.mock import patch
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for testing
import matplotlib.pyplot as plt

from beast_pype.outputs import (
    read_log_file,
    read_xml_set_logs_for_plotting,
    percentile_5th,
    percentile_95th,
    hdi_pivot,
    hdi_columns_starting_with,
    plot_box_violin,
    plot_hist_kde,
    plot_origin_or_tmrca,
    plot_comparative_box_violin,
    summary_stats_and_plot,
    _gridded_skyline,
)


# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def tmp_dir():
    """Create a temporary directory and clean up after test."""
    d = tempfile.mkdtemp()
    yield d
    shutil.rmtree(d)


@pytest.fixture
def synthetic_log_file(tmp_dir):
    """Create a synthetic .log file."""
    filepath = os.path.join(tmp_dir, "test.log")
    n = 100
    data = {
        "Sample": np.arange(0, n * 1000, 1000),
        "posterior": np.random.randn(n).cumsum(),
        "likelihood": np.random.randn(n).cumsum(),
        "treeHeight": np.abs(np.random.randn(n)) + 5,
    }
    df = pd.DataFrame(data)
    with open(filepath, "w") as f:
        f.write("# comment line\n")
    df.to_csv(filepath, sep="\t", index=False, mode="a")
    return filepath


@pytest.fixture
def synthetic_log_file_csv(tmp_dir):
    """Create a synthetic .csv log file."""
    filepath = os.path.join(tmp_dir, "test.csv")
    n = 100
    data = {
        "Sample": np.arange(0, n * 1000, 1000),
        "posterior": np.random.randn(n).cumsum(),
        "likelihood": np.random.randn(n).cumsum(),
    }
    df = pd.DataFrame(data)
    df.to_csv(filepath, index=False)
    return filepath


@pytest.fixture
def synthetic_log_file_tsv(tmp_dir):
    """Create a synthetic .tsv log file."""
    filepath = os.path.join(tmp_dir, "test.tsv")
    n = 100
    data = {
        "Sample": np.arange(0, n * 1000, 1000),
        "posterior": np.random.randn(n).cumsum(),
        "likelihood": np.random.randn(n).cumsum(),
    }
    df = pd.DataFrame(data)
    df.to_csv(filepath, sep="\t", index=False)
    return filepath


@pytest.fixture
def synthetic_log_with_origin(tmp_dir):
    """Create a synthetic .log file with origin column."""
    filepath = os.path.join(tmp_dir, "test_origin.log")
    n = 100
    youngest_tip = 2024.5
    data = {
        "Sample": np.arange(0, n * 1000, 1000),
        "posterior": np.random.randn(n).cumsum(),
        "likelihood": np.random.randn(n).cumsum(),
        "treeHeight": np.abs(np.random.randn(n)) + 5,
        "origin_BDSKY_Serial": np.abs(np.random.randn(n)) + 10,
    }
    df = pd.DataFrame(data)
    with open(filepath, "w") as f:
        f.write("# comment line\n")
    df.to_csv(filepath, sep="\t", index=False, mode="a")
    return filepath, youngest_tip


@pytest.fixture
def synthetic_log_with_sampling_prop(tmp_dir):
    """Create a synthetic .log file with samplingProportion columns."""
    filepath = os.path.join(tmp_dir, "test_sampling.log")
    n = 50
    data = {
        "Sample": np.arange(0, n * 1000, 1000),
        "posterior": np.random.randn(n).cumsum(),
        "samplingProportion.1": np.zeros(n),
        "samplingProportion.2": np.random.rand(n),
    }
    df = pd.DataFrame(data)
    with open(filepath, "w") as f:
        f.write("# comment\n")
    df.to_csv(filepath, sep="\t", index=False, mode="a")
    return filepath


@pytest.fixture
def synthetic_log_with_growth_rate(tmp_dir):
    """Create a synthetic .log file with growthRate column."""
    filepath = os.path.join(tmp_dir, "test_growth.log")
    n = 50
    data = {
        "Sample": np.arange(0, n * 1000, 1000),
        "posterior": np.random.randn(n).cumsum(),
        "treeHeight": np.abs(np.random.randn(n)) + 5,
        "growthRate": np.abs(np.random.randn(n)) + 0.5,
    }
    df = pd.DataFrame(data)
    with open(filepath, "w") as f:
        f.write("# comment\n")
    df.to_csv(filepath, sep="\t", index=False, mode="a")
    return filepath


@pytest.fixture
def trace_df():
    """Create a synthetic trace DataFrame for plotting tests."""
    np.random.seed(42)
    n = 200
    return pd.DataFrame({
        "posterior": np.random.randn(n).cumsum(),
        "likelihood": np.random.randn(n).cumsum(),
        "treeHeight": np.abs(np.random.randn(n)) + 5,
        "reproductiveNumber.1": np.abs(np.random.randn(n)) + 1,
        "reproductiveNumber.2": np.abs(np.random.randn(n)) + 0.8,
        "Origin": np.random.uniform(2020.0, 2022.0, n),
        "TMRCA": np.random.uniform(2019.0, 2021.0, n),
    })


@pytest.fixture
def xml_set_trace_df(trace_df):
    """Create a trace DataFrame with xml_set column."""
    df1 = trace_df.copy()
    df1["xml_set"] = "strain_A"
    df2 = trace_df.copy()
    df2["xml_set"] = "strain_B"
    return pd.concat([df1, df2], ignore_index=True)


@pytest.fixture
def melted_df(xml_set_trace_df):
    """Create a melted DataFrame for violin/box plots."""
    id_vars = ["xml_set"]
    value_vars = [c for c in xml_set_trace_df.columns if c != "xml_set"]
    return xml_set_trace_df.melt(
        id_vars=id_vars, value_vars=value_vars, value_name="Estimate"
    )


# ============================================================================
# Tests: read_log_file
# ============================================================================


class TestReadLogFile:
    """Tests for read_log_file function."""

    def test_reads_log_format(self, synthetic_log_file):
        df = read_log_file(synthetic_log_file)
        assert isinstance(df, pd.DataFrame)
        assert "Sample" in df.columns
        assert "posterior" in df.columns
        assert len(df) == 100

    def test_reads_csv_format(self, synthetic_log_file_csv):
        df = read_log_file(synthetic_log_file_csv)
        assert isinstance(df, pd.DataFrame)
        assert "Sample" in df.columns
        assert len(df) == 100

    def test_reads_tsv_format(self, synthetic_log_file_tsv):
        df = read_log_file(synthetic_log_file_tsv)
        assert isinstance(df, pd.DataFrame)
        assert "Sample" in df.columns
        assert len(df) == 100

    def test_invalid_extension_raises(self, tmp_dir):
        filepath = os.path.join(tmp_dir, "test.xyz")
        with open(filepath, "w") as f:
            f.write("data")
        with pytest.raises(ValueError, match="file_path should end with"):
            read_log_file(filepath)

    def test_removes_zero_sampling_prop(self, synthetic_log_with_sampling_prop):
        df = read_log_file(synthetic_log_with_sampling_prop)
        # First samplingProportion column (all zeros) should be removed
        assert "samplingProportion.1" not in df.columns
        assert "samplingProportion.2" in df.columns

    def test_keeps_sampling_prop_when_disabled(self, synthetic_log_with_sampling_prop):
        df = read_log_file(
            synthetic_log_with_sampling_prop,
            remove_0_first_sampling_prop=False,
        )
        assert "samplingProportion.1" in df.columns

    def test_origin_conversion(self, synthetic_log_with_origin):
        filepath, youngest_tip = synthetic_log_with_origin
        df = read_log_file(filepath, youngest_tip_date=youngest_tip)
        assert "Origin" in df.columns
        assert "TMRCA" in df.columns

    def test_growth_rate_conversion(self, synthetic_log_with_growth_rate):
        filepath = synthetic_log_with_growth_rate
        df = read_log_file(filepath, youngest_tip_date=2024.0)
        assert "Growth Rate per day" in df.columns
        assert "Doubling Time in days" in df.columns


# ============================================================================
# Tests: read_xml_set_logs_for_plotting
# ============================================================================


class TestReadXmlSetLogsForPlotting:
    """Tests for read_xml_set_logs_for_plotting."""

    def test_returns_df_and_melted(self, synthetic_log_file, tmp_dir):
        # Create a second log file
        filepath2 = os.path.join(tmp_dir, "test2.log")
        shutil.copy(synthetic_log_file, filepath2)
        file_dict = {"strain_A": synthetic_log_file, "strain_B": filepath2}
        df, df_melt = read_xml_set_logs_for_plotting(file_dict)
        assert isinstance(df, pd.DataFrame)
        assert isinstance(df_melt, pd.DataFrame)
        assert "xml_set" in df.columns
        assert "Estimate" in df_melt.columns

    def test_xml_set_label_custom(self, synthetic_log_file, tmp_dir):
        filepath2 = os.path.join(tmp_dir, "test2.log")
        shutil.copy(synthetic_log_file, filepath2)
        file_dict = {"A": synthetic_log_file, "B": filepath2}
        df, df_melt = read_xml_set_logs_for_plotting(
            file_dict, xml_set_label="lineage"
        )
        assert "lineage" in df.columns


# ============================================================================
# Tests: Utility functions
# ============================================================================


class TestPercentileFunctions:
    """Tests for percentile helper functions."""

    def test_percentile_5th(self):
        data = np.arange(100)
        assert percentile_5th(data) == pytest.approx(np.percentile(data, 5))

    def test_percentile_95th(self):
        data = np.arange(100)
        assert percentile_95th(data) == pytest.approx(np.percentile(data, 95))


class TestHdiPivot:
    """Tests for hdi_pivot."""

    def test_returns_dataframe(self, xml_set_trace_df):
        result = hdi_pivot(xml_set_trace_df, "treeHeight")
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2  # two xml_sets

    def test_contains_expected_columns(self, xml_set_trace_df):
        result = hdi_pivot(xml_set_trace_df, "treeHeight", hdi_prob=0.9)
        assert "Median" in result.columns
        assert "Lower 0.9 HDI" in result.columns
        assert "Upper 0.9 HDI" in result.columns


class TestHdiColumnsStartingWith:
    """Tests for hdi_columns_starting_with."""

    def test_returns_dataframe(self, trace_df):
        result = hdi_columns_starting_with(trace_df, "reproductiveNumber")
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2  # reproductiveNumber.1 and .2

    def test_has_expected_columns(self, trace_df):
        result = hdi_columns_starting_with(trace_df, "reproductiveNumber")
        assert "Parameter" in result.columns
        assert "Median" in result.columns


# ============================================================================
# Tests: Plotting functions
# ============================================================================


class TestPlotBoxViolin:
    """Tests for plot_box_violin."""

    def test_returns_axes(self, xml_set_trace_df):
        fig = plot_box_violin(xml_set_trace_df, x="xml_set", y="treeHeight")
        assert fig is not None
        plt.close("all")

    def test_with_grid(self, xml_set_trace_df):
        fig = plot_box_violin(
            xml_set_trace_df, x="xml_set", y="treeHeight", include_grid=True
        )
        assert fig is not None
        plt.close("all")

    def test_without_grid(self, xml_set_trace_df):
        fig = plot_box_violin(
            xml_set_trace_df, x="xml_set", y="treeHeight", include_grid=False
        )
        assert fig is not None
        plt.close("all")


class TestPlotHistKde:
    """Tests for plot_hist_kde."""

    def test_returns_fig_ax_hdi(self, trace_df):
        fig, ax, hdi_est = plot_hist_kde(trace_df, "treeHeight")
        assert fig is not None
        assert ax is not None
        assert "Median" in hdi_est
        assert "Lower 0.95 HDI" in hdi_est
        assert "Upper 0.95 HDI" in hdi_est
        plt.close("all")

    def test_no_hdi(self, trace_df):
        result = plot_hist_kde(trace_df, "treeHeight", hdi_prob=None)
        fig, ax = result
        assert fig is not None
        assert ax is not None
        plt.close("all")

    def test_custom_hdi_prob(self, trace_df):
        fig, ax, hdi_est = plot_hist_kde(trace_df, "treeHeight", hdi_prob=0.9)
        assert "Lower 0.9 HDI" in hdi_est
        assert "Upper 0.9 HDI" in hdi_est
        plt.close("all")

    def test_custom_xlabel(self, trace_df):
        fig, ax, _ = plot_hist_kde(trace_df, "treeHeight", x_label="My Label")
        assert ax.get_xlabel() == "My Label"
        plt.close("all")

    def test_default_xlabel_is_parameter(self, trace_df):
        fig, ax, _ = plot_hist_kde(trace_df, "treeHeight")
        assert ax.get_xlabel() == "treeHeight"
        plt.close("all")


class TestPlotOriginOrTmrca:
    """Tests for plot_origin_or_tmrca."""

    def test_origin_plot(self, trace_df):
        fig, ax, hdi_est = plot_origin_or_tmrca(trace_df, "Origin", hdi_prob=0.95)
        assert fig is not None
        plt.close("all")

    def test_tmrca_plot(self, trace_df):
        fig, ax, hdi_est = plot_origin_or_tmrca(trace_df, "TMRCA", hdi_prob=0.95)
        assert fig is not None
        plt.close("all")

    def test_invalid_parameter_raises(self, trace_df):
        with pytest.raises(ValueError, match='parameter must be "Origin" or "TMRCA"'):
            plot_origin_or_tmrca(trace_df, "posterior")

    def test_no_hdi(self, trace_df):
        result = plot_origin_or_tmrca(trace_df, "Origin", hdi_prob=None)
        fig, ax = result
        assert fig is not None
        plt.close("all")


class TestPlotComparativeBoxViolin:
    """Tests for plot_comparative_box_violin."""

    def test_basic_plot(self, melted_df):
        fig = plot_comparative_box_violin(melted_df, "treeHeight")
        assert fig is not None
        plt.close("all")

    def test_with_prior_draws(self, melted_df):
        prior = np.random.randn(100)
        fig = plot_comparative_box_violin(melted_df, "treeHeight", prior_draws=prior)
        assert fig is not None
        plt.close("all")


# ============================================================================
# Tests: _gridded_skyline
# ============================================================================


class TestGriddedSkyline:
    """Tests for _gridded_skyline internal function."""

    def test_returns_dataframe(self):
        n = 100
        data = {
            "reproductiveNumber.1": np.abs(np.random.randn(n)) + 1,
            "reproductiveNumber.2": np.abs(np.random.randn(n)) + 0.8,
            "reproductiveNumber.3": np.abs(np.random.randn(n)) + 0.6,
            "origin_BDSKY_Serial": np.abs(np.random.randn(n)) + 10,
            "Origin": np.random.uniform(2020, 2022, n),
        }
        df = pd.DataFrame(data)
        result = _gridded_skyline(df, youngest_tip=2024.0)
        assert isinstance(result, pd.DataFrame)
        assert "lower" in result.columns
        assert "median" in result.columns
        assert "upper" in result.columns
        assert "year_decimal" in result.columns

    def test_grid_size(self):
        n = 100
        data = {
            "reproductiveNumber.1": np.abs(np.random.randn(n)) + 1,
            "reproductiveNumber.2": np.abs(np.random.randn(n)) + 0.8,
            "origin_BDSKY_Serial": np.abs(np.random.randn(n)) + 10,
            "Origin": np.random.uniform(2020, 2022, n),
        }
        df = pd.DataFrame(data)
        result = _gridded_skyline(df, youngest_tip=2024.0, grid_size=50)
        assert len(result) == 50


# ============================================================================
# Tests: summary_stats_and_plot
# ============================================================================


class TestSummaryStatsAndPlot:
    """Tests for summary_stats_and_plot."""

    def test_returns_fig_df_stats(self, xml_set_trace_df):
        fig, y_df, stats = summary_stats_and_plot(
            xml_set_trace_df, x="xml_set", y="treeHeight"
        )
        assert fig is not None
        assert isinstance(y_df, pd.DataFrame)
        assert isinstance(stats, pd.DataFrame)
        plt.close("all")
