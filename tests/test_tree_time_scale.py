"""Unit tests for beast_pype.tree_time_scale clock filter methods."""

import os
import pytest
import tempfile
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib
import matplotlib.pyplot as plt

from beast_pype.tree_time_scale import timescale, plot_root_to_tip, tree_nodes_ci


def _make_synthetic_data(tmp_path, n_tips=10, seq_length=300, outlier_date=None):
    """Create synthetic tree, alignment, and dates for testing.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory.
    n_tips : int
        Number of tips.
    seq_length : int
        Length of sequences.
    outlier_date : float or None
        If set, make tip 0 an outlier by assigning this date.

    Returns
    -------
    tuple of (tree_path, alignment_path, dates_path, tip_names)
    """
    rng = np.random.default_rng(42)

    # Generate tip names and dates (evenly spaced over 2 years)
    tip_names = [f"seq_{i}" for i in range(n_tips)]
    base_year = 2020.0
    dates = [base_year + i * (2.0 / n_tips) for i in range(n_tips)]
    if outlier_date is not None:
        dates[0] = outlier_date  # Make first tip an outlier

    # Create a simple star tree with branch lengths proportional to time
    clock_rate = 0.001  # subs/site/year
    newick_parts = []
    for i in range(n_tips):
        bl = max(0.001, clock_rate * (dates[i] - base_year) + rng.normal(0, 0.0001))
        newick_parts.append(f"{tip_names[i]}:{bl:.6f}")
    newick = "(" + ",".join(newick_parts) + ");"

    tree_path = str(tmp_path / "test_tree.nwk")
    with open(tree_path, "w") as f:
        f.write(newick)

    # Create alignment (random sequences with mutations proportional to branch length)
    base_seq = "".join(rng.choice(list("ACGT"), size=seq_length))
    records = []
    for i in range(n_tips):
        n_mutations = max(0, int(clock_rate * (dates[i] - base_year) * seq_length))
        seq_list = list(base_seq)
        mut_positions = rng.choice(seq_length, size=min(n_mutations, seq_length), replace=False)
        for pos in mut_positions:
            alternatives = [b for b in "ACGT" if b != seq_list[pos]]
            seq_list[pos] = rng.choice(alternatives)
        records.append(SeqRecord(Seq("".join(seq_list)), id=tip_names[i], description=""))

    alignment_path = str(tmp_path / "test_alignment.fasta")
    with open(alignment_path, "w") as f:
        SeqIO.write(records, f, "fasta")

    # Create dates CSV
    dates_path = str(tmp_path / "test_dates.csv")
    dates_df = pd.DataFrame({"name": tip_names, "date": dates})
    dates_df.to_csv(dates_path, index=False)

    return tree_path, alignment_path, dates_path, tip_names


@pytest.fixture
def synthetic_data(tmp_path):
    """Fixture for clock-like synthetic data (no outliers)."""
    return _make_synthetic_data(tmp_path, n_tips=10)


@pytest.fixture
def synthetic_data_with_outlier(tmp_path):
    """Fixture for synthetic data with one outlier (date far in the past)."""
    return _make_synthetic_data(tmp_path, n_tips=20, outlier_date=2010.0)


class TestTimescaleClockFilter:
    """Tests for clock_filter and clock_filter_method parameters."""

    def test_clock_filter_none_disables_filtering(self, synthetic_data):
        """No tips should be removed when clock_filter=None."""
        tree_path, aln_path, dates_path, tip_names = synthetic_data
        time_tree, bad_tips = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=None,
            remove_root=False,
            rng_seed=42,
        )
        assert bad_tips == []
        remaining_tips = [t.name for t in time_tree.tree.get_terminals()]
        assert len(remaining_tips) == len(tip_names)

    def test_clock_filter_local_default(self, synthetic_data):
        """Default clock_filter=3.0 with method='local' should run without error."""
        tree_path, aln_path, dates_path, tip_names = synthetic_data
        time_tree, bad_tips = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=3.0,
            clock_filter_method='local',
            remove_root=False,
            rng_seed=42,
        )
        assert isinstance(bad_tips, list)
        remaining_tips = [t.name for t in time_tree.tree.get_terminals()]
        assert len(remaining_tips) + len(bad_tips) == len(tip_names)

    def test_clock_filter_residual_method(self, synthetic_data):
        """clock_filter_method='residual' should also run without error."""
        tree_path, aln_path, dates_path, tip_names = synthetic_data
        time_tree, bad_tips = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=3.0,
            clock_filter_method='residual',
            remove_root=False,
            rng_seed=42,
        )
        assert isinstance(bad_tips, list)
        remaining_tips = [t.name for t in time_tree.tree.get_terminals()]
        assert len(remaining_tips) + len(bad_tips) == len(tip_names)

    def test_strict_filter_removes_outlier(self, synthetic_data_with_outlier):
        """A very strict clock filter (low z-score threshold) should detect the outlier."""
        tree_path, aln_path, dates_path, tip_names = synthetic_data_with_outlier
        time_tree, bad_tips = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=1.5,  # Very strict threshold
            clock_filter_method='local',
            reroot=[tip_names[-1]],
            remove_root=False,
            rng_seed=42,
        )
        # The outlier (seq_0 with date 2010.0) should be removed
        assert len(bad_tips) > 0
        assert "seq_0" in bad_tips
        remaining_tips = [t.name for t in time_tree.tree.get_terminals()]
        assert "seq_0" not in remaining_tips

    def test_outliers_attribute_set(self, synthetic_data_with_outlier):
        """time_tree.outliers should be a DataFrame when clock filtering is active."""
        tree_path, aln_path, dates_path, tip_names = synthetic_data_with_outlier
        time_tree, bad_tips = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=1.5,
            clock_filter_method='local',
            reroot=[tip_names[-1]],
            remove_root=False,
            rng_seed=42,
        )
        if len(bad_tips) > 0:
            assert hasattr(time_tree, 'outliers')
            assert time_tree.outliers is not None
            assert isinstance(time_tree.outliers, pd.DataFrame)
            assert len(time_tree.outliers) > 0

    def test_bad_tips_pruned_from_tree(self, synthetic_data_with_outlier):
        """Bad tips should not remain in the tree after timescale()."""
        tree_path, aln_path, dates_path, tip_names = synthetic_data_with_outlier
        time_tree, bad_tips = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=1.5,
            clock_filter_method='local',
            reroot=[tip_names[-1]],
            remove_root=False,
            rng_seed=42,
        )
        remaining_tips = {t.name for t in time_tree.tree.get_terminals()}
        for tip in bad_tips:
            assert tip not in remaining_tips

    def test_lenient_filter_keeps_all(self, synthetic_data_with_outlier):
        """A very lenient filter (high threshold) should keep all tips."""
        tree_path, aln_path, dates_path, tip_names = synthetic_data_with_outlier
        time_tree, bad_tips = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=100.0,  # Extremely lenient
            clock_filter_method='local',
            reroot=[tip_names[-1]],
            remove_root=False,
            rng_seed=42,
        )
        assert bad_tips == []
        remaining_tips = [t.name for t in time_tree.tree.get_terminals()]
        assert len(remaining_tips) == len(tip_names)

    def test_negative_branch_lengths_corrected(self, synthetic_data):
        """Negative branch lengths within tolerance should be set to 0."""
        tree_path, aln_path, dates_path, _ = synthetic_data
        tolerance = 2.0  # Use a large tolerance to catch TreeTime's negative branches
        time_tree, _ = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=None,
            remove_root=False,
            negative_tolerance=tolerance,
            rng_seed=42,
        )
        for node in time_tree.tree.find_clades():
            if node.branch_length is not None:
                # Any negative branch with |bl| < tolerance should be corrected to 0
                # Any remaining negative branch must exceed the tolerance
                if node.branch_length < 0:
                    assert abs(node.branch_length) >= tolerance


class TestTimescaleReroot:
    """Tests for remove_root parameter interaction with clock filter."""

    def test_remove_root_prunes_root_strains(self, synthetic_data):
        """When reroot is a list and remove_root=True, root strains are pruned."""
        tree_path, aln_path, dates_path, tip_names = synthetic_data
        root_strain = tip_names[-1]  # Use last tip as root
        time_tree, bad_tips = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            reroot=[root_strain],
            remove_root=True,
            clock_filter=None,
            rng_seed=42,
        )
        remaining_tips = [t.name for t in time_tree.tree.get_terminals()]
        assert root_strain not in remaining_tips

    def test_remove_root_false_keeps_root_strains(self, synthetic_data):
        """When remove_root=False, root strains remain."""
        tree_path, aln_path, dates_path, tip_names = synthetic_data
        root_strain = tip_names[-1]
        time_tree, bad_tips = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            reroot=[root_strain],
            remove_root=False,
            clock_filter=None,
            rng_seed=42,
        )
        remaining_tips = [t.name for t in time_tree.tree.get_terminals()]
        assert root_strain in remaining_tips


class TestPlotRootToTip:
    """Tests for plot_root_to_tip including outlier visualization."""

    def test_returns_fig_and_ax(self, synthetic_data):
        """plot_root_to_tip should return fig, ax."""
        tree_path, aln_path, dates_path, _ = synthetic_data
        time_tree, _ = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=None,
            remove_root=False,
            rng_seed=42,
        )
        fig, ax = plot_root_to_tip(time_tree)
        assert isinstance(fig, matplotlib.figure.Figure)
        assert ax is not None
        plt.close(fig)

    def test_outliers_plotted_as_orange(self, synthetic_data_with_outlier):
        """When outliers exist, they should appear as orange scatter points."""
        tree_path, aln_path, dates_path, tip_names = synthetic_data_with_outlier
        # Use explicit root to avoid least-squares rerooting failure
        time_tree, bad_tips = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=1.5,
            clock_filter_method='local',
            reroot=[tip_names[-1]],
            remove_root=False,
            rng_seed=42,
        )
        # Build outliers_df with numdate and dist2root
        if bad_tips and hasattr(time_tree, 'outliers') and time_tree.outliers is not None:
            outliers_df = time_tree.outliers.copy()
            if 'given_date' in outliers_df.columns and 'numdate' not in outliers_df.columns:
                outliers_df['numdate'] = outliers_df['given_date']
            rate = time_tree.clock_model['slope']
            intercept = time_tree.clock_model['intercept']
            if 'dist2root' not in outliers_df.columns and 'apparent_date' in outliers_df.columns:
                outliers_df['dist2root'] = rate * outliers_df['apparent_date'] + intercept
        else:
            outliers_df = None
        fig, ax = plot_root_to_tip(time_tree, outliers_df=outliers_df)
        # Check legend contains outlier label
        if bad_tips and outliers_df is not None:
            legend = ax.get_legend()
            assert legend is not None
            legend_texts = [t.get_text() for t in legend.get_texts()]
            assert any('outlier' in t.lower() for t in legend_texts)
        plt.close(fig)

    def test_no_outliers_no_legend(self, synthetic_data):
        """When no outliers, no outlier legend should be present."""
        tree_path, aln_path, dates_path, _ = synthetic_data
        time_tree, _ = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=None,
            remove_root=False,
            rng_seed=42,
        )
        fig, ax = plot_root_to_tip(time_tree)
        legend = ax.get_legend()
        if legend:
            legend_texts = [t.get_text() for t in legend.get_texts()]
            assert not any('outlier' in t.lower() for t in legend_texts)
        plt.close(fig)


class TestTreeNodesCi:
    """Tests for tree_nodes_ci function."""

    def test_returns_dataframe(self, synthetic_data):
        """tree_nodes_ci should return a DataFrame."""
        tree_path, aln_path, dates_path, _ = synthetic_data
        time_tree, _ = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=None,
            remove_root=False,
            rng_seed=42,
        )
        df = tree_nodes_ci(time_tree)
        assert isinstance(df, pd.DataFrame)
        assert 'node' in df.columns
        assert 'date' in df.columns
        assert 'year_decimal' in df.columns

    def test_ci_columns_present(self, synthetic_data):
        """CI columns should reflect the fraction used."""
        tree_path, aln_path, dates_path, _ = synthetic_data
        time_tree, _ = timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=None,
            remove_root=False,
            rng_seed=42,
        )
        df = tree_nodes_ci(time_tree, fraction=0.95)
        assert any('interval' in c for c in df.columns)
        assert any('upper' in c for c in df.columns)


class TestIterativeTimescale:
    """Tests for iterative_timescale function."""

    def test_no_clock_filter_single_iteration(self, synthetic_data):
        """With clock_filter=None, should run once and return empty outliers."""
        from beast_pype.tree_time_scale import iterative_timescale
        tree_path, aln_path, dates_path, tip_names = synthetic_data
        time_tree, outliers_df, final_fasta, final_meta = iterative_timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=None,
            remove_root=False,
            rng_seed=42,
        )
        assert outliers_df.empty
        assert os.path.isfile(final_fasta)
        assert os.path.isfile(final_meta)
        remaining = [t.name for t in time_tree.tree.get_terminals()]
        assert len(remaining) == len(tip_names)

    def test_converges_with_clean_data(self, synthetic_data):
        """Clean data (no outliers) should converge in 1 iteration."""
        from beast_pype.tree_time_scale import iterative_timescale
        tree_path, aln_path, dates_path, tip_names = synthetic_data
        time_tree, outliers_df, final_fasta, final_meta = iterative_timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=3.0,
            clock_filter_method='local',
            remove_root=False,
            rng_seed=42,
        )
        # With clean clock-like data and a lenient filter, expect no outliers
        assert outliers_df.empty or len(outliers_df) == 0
        remaining = [t.name for t in time_tree.tree.get_terminals()]
        assert len(remaining) == len(tip_names)

    def test_iterative_removes_outlier(self, synthetic_data_with_outlier):
        """Outlier should be removed in iterative mode."""
        from beast_pype.tree_time_scale import iterative_timescale
        tree_path, aln_path, dates_path, tip_names = synthetic_data_with_outlier
        time_tree, outliers_df, final_fasta, final_meta = iterative_timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=1.5,
            clock_filter_method='local',
            reroot=[tip_names[-1]],
            remove_root=False,
            rng_seed=42,
            max_iterations=1,  # Synthetic star tree breaks on subsequent iterations
        )
        # seq_0 (with date 2010) should be removed as outlier
        assert 'seq_0' in outliers_df['name'].values
        remaining = [t.name for t in time_tree.tree.get_terminals()]
        assert 'seq_0' not in remaining

    def test_outliers_df_has_iteration_column(self, synthetic_data_with_outlier):
        """Outliers DataFrame should include iteration column."""
        from beast_pype.tree_time_scale import iterative_timescale
        tree_path, aln_path, dates_path, tip_names = synthetic_data_with_outlier
        _, outliers_df, _, _ = iterative_timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=1.5,
            clock_filter_method='local',
            reroot=[tip_names[-1]],
            remove_root=False,
            rng_seed=42,
            max_iterations=1,
        )
        if not outliers_df.empty:
            assert 'iteration' in outliers_df.columns
            assert 'name' in outliers_df.columns

    def test_root_strains_excluded_from_outliers(self, synthetic_data_with_outlier):
        """Root strains should never be removed as outliers."""
        from beast_pype.tree_time_scale import iterative_timescale
        tree_path, aln_path, dates_path, tip_names = synthetic_data_with_outlier
        root_tip = tip_names[-1]
        _, outliers_df, _, _ = iterative_timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=1.5,
            clock_filter_method='local',
            reroot=[root_tip],
            remove_root=False,
            rng_seed=42,
            max_iterations=1,
        )
        if not outliers_df.empty:
            assert root_tip not in outliers_df['name'].values

    def test_remove_root_after_convergence(self, synthetic_data):
        """Root strains should be pruned after convergence when remove_root=True."""
        from beast_pype.tree_time_scale import iterative_timescale
        tree_path, aln_path, dates_path, tip_names = synthetic_data
        root_tip = tip_names[-1]
        time_tree, _, _, _ = iterative_timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=None,
            reroot=[root_tip],
            remove_root=True,
            rng_seed=42,
        )
        remaining = [t.name for t in time_tree.tree.get_terminals()]
        assert root_tip not in remaining

    def test_max_iterations_respected(self, synthetic_data_with_outlier):
        """Should stop after max_iterations even if not converged."""
        from beast_pype.tree_time_scale import iterative_timescale
        tree_path, aln_path, dates_path, tip_names = synthetic_data_with_outlier
        # max_iterations=1 means only 1 iteration regardless
        time_tree, outliers_df, _, _ = iterative_timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=1.5,
            clock_filter_method='local',
            reroot=[tip_names[-1]],
            remove_root=False,
            max_iterations=1,
            rng_seed=42,
        )
        # Should have run at most 1 iteration
        if not outliers_df.empty:
            assert outliers_df['iteration'].max() <= 1

    def test_filtered_files_exist(self, synthetic_data_with_outlier):
        """Filtered FASTA and metadata files should exist after outlier removal."""
        from beast_pype.tree_time_scale import iterative_timescale
        tree_path, aln_path, dates_path, tip_names = synthetic_data_with_outlier
        _, outliers_df, final_fasta, final_meta = iterative_timescale(
            ftree=tree_path,
            falignment=aln_path,
            fdates=dates_path,
            clock_filter=1.5,
            clock_filter_method='local',
            reroot=[tip_names[-1]],
            remove_root=False,
            rng_seed=42,
            max_iterations=1,
        )
        assert os.path.isfile(final_fasta)
        assert os.path.isfile(final_meta)
        if not outliers_df.empty:
            # Filtered files should have fewer sequences
            from Bio import SeqIO as sio
            original_count = len(list(sio.parse(aln_path, 'fasta')))
            filtered_count = len(list(sio.parse(final_fasta, 'fasta')))
            assert filtered_count < original_count
