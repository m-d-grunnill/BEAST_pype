"""Unit tests for beast_pype.diagnostics.mcmc module."""
import os
import pytest
import tempfile
import shutil
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

from beast_pype.diagnostics.mcmc import (
    read_log_file_as_dataframe,
    read_log_files_as_xarraydata,
    read_log_files_as_posterior,
    burn_posterior,
    select_chains_from_posterior,
    subset_and_merge_trees,
    merge_logs_to_csv,
    BEASTDiag,
)


# ============================================================================
# Fixtures
# ============================================================================

EXAMPLE_DIR = os.path.join(
    os.path.dirname(__file__),
    "..",
    "example_beast2_outputs_to_diagnose",
)


@pytest.fixture
def example_log_files():
    """Return sorted list of example .log file paths."""
    log_dir = os.path.abspath(EXAMPLE_DIR)
    files = sorted(
        os.path.join(log_dir, f)
        for f in os.listdir(log_dir)
        if f.endswith(".log")
    )
    assert len(files) > 0, "No .log files found in example directory"
    return files


@pytest.fixture
def example_trees_files():
    """Return sorted list of example .trees file paths."""
    trees_dir = os.path.abspath(EXAMPLE_DIR)
    files = sorted(
        os.path.join(trees_dir, f)
        for f in os.listdir(trees_dir)
        if f.endswith(".trees")
    )
    assert len(files) > 0, "No .trees files found in example directory"
    return files


@pytest.fixture
def tmp_dir():
    """Create a temporary directory and clean up after test."""
    d = tempfile.mkdtemp()
    yield d
    shutil.rmtree(d)


@pytest.fixture
def synthetic_log_files(tmp_dir):
    """Create synthetic .log files for testing."""
    n_samples = 100
    step_size = 1000
    files = []
    for seed in [111, 222, 333]:
        filepath = os.path.join(tmp_dir, f"run-with-seed-{seed}-BEAST.log")
        samples = np.arange(0, n_samples * step_size, step_size)
        data = {
            "Sample": samples,
            "posterior": np.random.randn(n_samples).cumsum(),
            "likelihood": np.random.randn(n_samples).cumsum(),
            "prior": np.random.randn(n_samples).cumsum(),
            "treeLikelihood": np.random.randn(n_samples).cumsum(),
            "clockRate": np.abs(np.random.randn(n_samples)) * 1e-3,
        }
        df = pd.DataFrame(data)
        # Write with comment header like BEAST 2 format
        with open(filepath, "w") as f:
            f.write("# comment line 1\n")
            f.write("# comment line 2\n")
            df.to_csv(f, sep="\t", index=False)
        files.append(filepath)
    return files


@pytest.fixture
def synthetic_trees_files(tmp_dir):
    """Create synthetic .trees (NEXUS) files for testing."""
    header = (
        "#NEXUS\n"
        "\n"
        "Begin taxa;\n"
        "\tDimensions ntax=3;\n"
        "\t\tTaxlabels\n"
        "\t\t\ttip_A\n"
        "\t\t\ttip_B\n"
        "\t\t\ttip_C\n"
        "\t\t\t;\n"
        "End;\n"
        "Begin trees;\n"
        "\tTranslate\n"
        "\t\t1 tip_A,\n"
        "\t\t2 tip_B,\n"
        "\t\t3 tip_C\n"
        "\t\t;\n"
    )
    footer = "End;\n"
    files = []
    for seed in [111, 222]:
        filepath = os.path.join(tmp_dir, f"run-with-seed-{seed}.trees")
        with open(filepath, "w") as f:
            f.write(header)
            for state in range(0, 100000, 10000):
                f.write(
                    f"tree STATE_{state} = "
                    f"((1:{0.1 + state * 0.00001:.5f},2:{0.1 + state * 0.00001:.5f}):"
                    f"{0.2 + state * 0.00001:.5f},3:{0.3 + state * 0.00001:.5f});\n"
                )
            f.write(footer)
        files.append(filepath)
    return files


@pytest.fixture
def synthetic_posterior(synthetic_log_files):
    """Create a posterior DataTree from synthetic log files."""
    return read_log_files_as_posterior(synthetic_log_files)


# ============================================================================
# Tests for read_log_file_as_dataframe
# ============================================================================

class TestReadLogFileAsDataframe:
    def test_returns_dataframe(self, example_log_files):
        df = read_log_file_as_dataframe(example_log_files[0])
        assert isinstance(df, pd.DataFrame)

    def test_has_sample_column(self, example_log_files):
        df = read_log_file_as_dataframe(example_log_files[0])
        assert "Sample" in df.columns

    def test_cut_to_first(self, example_log_files):
        df_full = read_log_file_as_dataframe(example_log_files[0])
        cut_value = df_full["Sample"].iloc[5]
        df_cut = read_log_file_as_dataframe(example_log_files[0], cut_to_first=cut_value)
        assert len(df_cut) <= len(df_full)
        assert df_cut["Sample"].max() <= cut_value

    def test_cut_to_first_none(self, example_log_files):
        df = read_log_file_as_dataframe(example_log_files[0], cut_to_first=None)
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

    def test_synthetic_log(self, synthetic_log_files):
        df = read_log_file_as_dataframe(synthetic_log_files[0])
        assert isinstance(df, pd.DataFrame)
        assert "Sample" in df.columns
        assert len(df) == 100


# ============================================================================
# Tests for read_log_files_as_xarraydata
# ============================================================================

class TestReadLogFilesAsXarrayData:
    def test_returns_dataset(self, synthetic_log_files):
        ds = read_log_files_as_xarraydata(synthetic_log_files)
        assert isinstance(ds, xr.Dataset)

    def test_has_chain_and_draw_dims(self, synthetic_log_files):
        ds = read_log_files_as_xarraydata(synthetic_log_files)
        assert "chain" in ds.dims
        assert "draw" in ds.dims

    def test_chain_count_from_list(self, synthetic_log_files):
        ds = read_log_files_as_xarraydata(synthetic_log_files)
        assert ds.sizes["chain"] == 3

    def test_chain_count_from_dict(self, synthetic_log_files):
        log_dict = {f"chain_{i}": f for i, f in enumerate(synthetic_log_files)}
        ds = read_log_files_as_xarraydata(log_dict)
        assert ds.sizes["chain"] == 3

    def test_start_index(self, synthetic_log_files):
        ds = read_log_files_as_xarraydata(synthetic_log_files, start=5)
        chains = ds.coords["chain"].values
        assert chains[0] == 5

    def test_invalid_type_raises(self):
        with pytest.raises(TypeError):
            read_log_files_as_xarraydata("not_a_list_or_dict")

    def test_dict_keys_as_chain_labels(self, synthetic_log_files):
        log_dict = {"alpha": synthetic_log_files[0], "beta": synthetic_log_files[1]}
        ds = read_log_files_as_xarraydata(log_dict)
        chains = list(ds.coords["chain"].values)
        assert "alpha" in chains
        assert "beta" in chains


# ============================================================================
# Tests for read_log_files_as_posterior
# ============================================================================

class TestReadLogFilesAsPosterior:
    def test_returns_datatree(self, synthetic_log_files):
        posterior = read_log_files_as_posterior(synthetic_log_files)
        assert isinstance(posterior, xr.DataTree)

    def test_has_posterior_group(self, synthetic_log_files):
        posterior = read_log_files_as_posterior(synthetic_log_files)
        assert hasattr(posterior, "posterior")

    def test_posterior_has_data_vars(self, synthetic_log_files):
        posterior = read_log_files_as_posterior(synthetic_log_files)
        data_vars = list(posterior.posterior.data_vars)
        assert "posterior" in data_vars
        assert "likelihood" in data_vars

    def test_with_dict_input(self, synthetic_log_files):
        log_dict = {"chain_a": synthetic_log_files[0], "chain_b": synthetic_log_files[1]}
        posterior = read_log_files_as_posterior(log_dict)
        assert isinstance(posterior, xr.DataTree)
        chains = list(posterior.posterior.coords["chain"].values)
        assert "chain_a" in chains
        assert "chain_b" in chains


# ============================================================================
# Tests for burn_posterior
# ============================================================================

class TestBurnPosterior:
    def test_burn_by_percentage(self, synthetic_posterior):
        burned = burn_posterior(synthetic_posterior, in_percentage=10)
        original_draws = len(synthetic_posterior.posterior["draw"])
        burned_draws = len(burned.posterior["draw"])
        expected = original_draws - round(0.10 * original_draws)
        assert burned_draws == expected

    def test_burn_by_proportion(self, synthetic_posterior):
        burned = burn_posterior(synthetic_posterior, in_proportion=0.2)
        original_draws = len(synthetic_posterior.posterior["draw"])
        burned_draws = len(burned.posterior["draw"])
        expected = original_draws - round(0.2 * original_draws)
        assert burned_draws == expected

    def test_front_by_percentage(self, synthetic_posterior):
        burned = burn_posterior(synthetic_posterior, front_percentage=50)
        original_draws = len(synthetic_posterior.posterior["draw"])
        burned_draws = len(burned.posterior["draw"])
        expected = round(0.50 * original_draws)
        assert burned_draws == expected

    def test_burn_in_and_front(self, synthetic_posterior):
        burned = burn_posterior(synthetic_posterior, in_percentage=10, front_percentage=80)
        original_draws = len(synthetic_posterior.posterior["draw"])
        in_num = round(0.10 * original_draws)
        front_num = round(0.80 * original_draws)
        expected = front_num - in_num
        burned_draws = len(burned.posterior["draw"])
        assert burned_draws == expected

    def test_no_burn(self, synthetic_posterior):
        burned = burn_posterior(synthetic_posterior)
        original_draws = len(synthetic_posterior.posterior["draw"])
        burned_draws = len(burned.posterior["draw"])
        assert burned_draws == original_draws

    def test_multiple_in_params_raises(self, synthetic_posterior):
        with pytest.raises(ValueError):
            burn_posterior(synthetic_posterior, in_proportion=0.1, in_percentage=10)

    def test_multiple_front_params_raises(self, synthetic_posterior):
        with pytest.raises(ValueError):
            burn_posterior(synthetic_posterior, front_proportion=0.5, front_percentage=50)

    def test_invalid_proportion_type_raises(self, synthetic_posterior):
        with pytest.raises(TypeError):
            burn_posterior(synthetic_posterior, in_proportion=10)  # should be float

    def test_invalid_percentage_type_raises(self, synthetic_posterior):
        with pytest.raises(TypeError):
            burn_posterior(synthetic_posterior, in_percentage="ten")


# ============================================================================
# Tests for select_chains_from_posterior
# ============================================================================

class TestSelectChainsFromPosterior:
    def test_select_single_chain(self, synthetic_posterior):
        selected = select_chains_from_posterior(synthetic_posterior, [0])
        chains = list(selected.posterior.coords["chain"].values)
        assert len(chains) == 1
        assert chains[0] == 0

    def test_select_multiple_chains(self, synthetic_posterior):
        selected = select_chains_from_posterior(synthetic_posterior, [0, 2])
        chains = list(selected.posterior.coords["chain"].values)
        assert len(chains) == 2
        assert 0 in chains
        assert 2 in chains

    def test_select_with_string_labels(self, synthetic_log_files):
        log_dict = {"alpha": synthetic_log_files[0], "beta": synthetic_log_files[1]}
        posterior = read_log_files_as_posterior(log_dict)
        selected = select_chains_from_posterior(posterior, ["alpha"])
        chains = list(selected.posterior.coords["chain"].values)
        assert chains == ["alpha"]


# ============================================================================
# Tests for subset_and_merge_trees
# ============================================================================

class TestSubsetAndMergeTrees:
    def test_basic_subset(self, synthetic_trees_files, tmp_dir):
        output_file = os.path.join(tmp_dir, "merged.trees")
        subset_and_merge_trees(
            file_list=synthetic_trees_files,
            in_number=20000,
            front_number=80000,
            output_file=output_file,
        )
        assert os.path.isfile(output_file)

    def test_output_contains_correct_states(self, synthetic_trees_files, tmp_dir):
        output_file = os.path.join(tmp_dir, "merged.trees")
        subset_and_merge_trees(
            file_list=synthetic_trees_files,
            in_number=20000,
            front_number=80000,
            output_file=output_file,
        )
        with open(output_file) as f:
            content = f.read()
        # Should contain states 20000-80000 from both files
        assert "STATE_" in content

    def test_relabeling(self, synthetic_trees_files, tmp_dir):
        output_file = os.path.join(tmp_dir, "merged.trees")
        subset_and_merge_trees(
            file_list=synthetic_trees_files,
            in_number=20000,
            front_number=80000,
            output_file=output_file,
            relabel_start=0,
        )
        with open(output_file) as f:
            content = f.read()
        assert "STATE_0" in content

    def test_in_number_greater_than_front_raises(self, synthetic_trees_files, tmp_dir):
        output_file = os.path.join(tmp_dir, "merged.trees")
        with pytest.raises(ValueError, match="in_number must be <= front_number"):
            subset_and_merge_trees(
                file_list=synthetic_trees_files,
                in_number=90000,
                front_number=10000,
                output_file=output_file,
            )

    def test_file_not_found_raises(self, tmp_dir):
        output_file = os.path.join(tmp_dir, "merged.trees")
        with pytest.raises(FileNotFoundError):
            subset_and_merge_trees(
                file_list=["nonexistent.trees"],
                in_number=0,
                front_number=100000,
                output_file=output_file,
            )

    def test_wrong_extension_raises(self, tmp_dir):
        bad_file = os.path.join(tmp_dir, "bad.txt")
        with open(bad_file, "w") as f:
            f.write("not a trees file")
        output_file = os.path.join(tmp_dir, "merged.trees")
        with pytest.raises(ValueError, match="does not have .trees or .tree extension"):
            subset_and_merge_trees(
                file_list=[bad_file],
                in_number=0,
                front_number=100000,
                output_file=output_file,
            )

    def test_header_mismatch_raises(self, tmp_dir):
        """Test that files with different headers raise ValueError."""
        header_a = "#NEXUS\nBegin taxa;\n\tDimensions ntax=3;\nEnd;\nBegin trees;\n"
        header_b = "#NEXUS\nBegin taxa;\n\tDimensions ntax=4;\nEnd;\nBegin trees;\n"
        footer = "End;\n"

        file_a = os.path.join(tmp_dir, "a.trees")
        file_b = os.path.join(tmp_dir, "b.trees")

        with open(file_a, "w") as f:
            f.write(header_a)
            f.write("tree STATE_0 = ((1:0.1,2:0.1):0.2,3:0.3);\n")
            f.write(footer)
        with open(file_b, "w") as f:
            f.write(header_b)
            f.write("tree STATE_0 = ((1:0.1,2:0.1):0.2,3:0.3);\n")
            f.write(footer)

        output_file = os.path.join(tmp_dir, "merged.trees")
        with pytest.raises(ValueError, match="NEXUS header mismatch"):
            subset_and_merge_trees(
                file_list=[file_a, file_b],
                in_number=0,
                front_number=100000,
                output_file=output_file,
            )

    def test_nexus_structure_preserved(self, synthetic_trees_files, tmp_dir):
        """Output file should start with #NEXUS and end with End;"""
        output_file = os.path.join(tmp_dir, "merged.trees")
        subset_and_merge_trees(
            file_list=synthetic_trees_files,
            in_number=0,
            front_number=90000,
            output_file=output_file,
        )
        with open(output_file) as f:
            lines = f.readlines()
        assert lines[0].strip() == "#NEXUS"
        assert lines[-1].strip() == "End;"

    def test_tree_count(self, synthetic_trees_files, tmp_dir):
        """Check correct number of trees selected."""
        output_file = os.path.join(tmp_dir, "merged.trees")
        # Each file has states 0, 10000, 20000, ..., 90000 (10 trees)
        # Selecting 20000-80000 gives 7 trees per file, 14 total from 2 files
        subset_and_merge_trees(
            file_list=synthetic_trees_files,
            in_number=20000,
            front_number=80000,
            output_file=output_file,
        )
        with open(output_file) as f:
            content = f.read()
        import re
        tree_count = len(re.findall(r"tree STATE_\d+", content))
        assert tree_count == 14  # 7 per file * 2 files


# ============================================================================
# Tests for merge_logs_to_csv
# ============================================================================

class TestMergeLogsToCsv:
    def test_creates_csv(self, synthetic_posterior, tmp_dir):
        output_file = os.path.join(tmp_dir, "merged.csv")
        merge_logs_to_csv(synthetic_posterior, output_file=output_file)
        assert os.path.isfile(output_file)

    def test_csv_content(self, synthetic_posterior, tmp_dir):
        output_file = os.path.join(tmp_dir, "merged.csv")
        merge_logs_to_csv(synthetic_posterior, output_file=output_file)
        df = pd.read_csv(output_file)
        assert "Sample" in df.columns
        assert "posterior" in df.columns

    def test_like_logcombiner_false(self, synthetic_posterior, tmp_dir):
        output_file = os.path.join(tmp_dir, "merged.csv")
        merge_logs_to_csv(
            synthetic_posterior, output_file=output_file, like_logcombiner=False
        )
        df = pd.read_csv(output_file)
        assert os.path.isfile(output_file)
        assert isinstance(df, pd.DataFrame)


# ============================================================================
# Tests for BEASTDiag class
# ============================================================================

class TestBEASTDiag:
    def test_init(self):
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        assert diag.directory == os.path.abspath(EXAMPLE_DIR)

    def test_init_missing_dir_raises(self, tmp_dir):
        with pytest.raises(FileNotFoundError):
            BEASTDiag(os.path.join(tmp_dir, "nonexistent_dir"))

    def test_has_parameters(self):
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        assert len(diag.parameters) > 0

    def test_has_chains(self):
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        assert len(diag.original_chains) > 0

    def test_chain_names_cleaned(self):
        """Chain names should not end with -BEAST or .log."""
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        for chain in diag.original_chains:
            assert not chain.endswith("-BEAST")
            assert not chain.endswith(".log")

    def test_set_burnin(self):
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        diag.set_burnin(in_percentage=20, front_percentage=80)
        assert diag.burinin_percentage == 20
        assert diag.keep_front_percentage == 80
        assert diag.selected_posterior is not None

    def test_set_burnin_invalid_percentage_raises(self):
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        with pytest.raises(TypeError):
            diag.set_burnin(in_percentage=-5)

    def test_set_burnin_front_less_than_in_raises(self):
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        with pytest.raises(ValueError):
            diag.set_burnin(in_percentage=50, front_percentage=30)

    def test_select_chains_with_list(self):
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        chains_to_select = diag.original_chains[:3]
        diag.select_chains(chains=chains_to_select)
        assert diag.selected_chains == chains_to_select
        assert diag.selected_posterior is not None

    def test_select_chains_with_kwargs(self):
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        kwargs = {chain: (i < 3) for i, chain in enumerate(diag.original_chains)}
        diag.select_chains(**kwargs)
        assert len(diag.selected_chains) == 3

    def test_select_chains_invalid_chain_raises(self):
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        with pytest.raises(ValueError):
            diag.select_chains(chains=["nonexistent_chain"])

    def test_select_chains_no_args_raises(self):
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        with pytest.raises(ValueError):
            diag.select_chains()

    def test_select_chains_both_args_raises(self):
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        kwargs = {chain: True for chain in diag.original_chains}
        with pytest.raises(ValueError):
            diag.select_chains(chains=diag.original_chains[:2], **kwargs)

    def test_output_prefix(self, tmp_dir):
        """Test that output_prefix is set correctly."""
        src = os.path.abspath(EXAMPLE_DIR)
        diag = BEASTDiag(src, output_prefix=os.path.join(tmp_dir, "test_"))
        assert diag.output_prefix == os.path.join(tmp_dir, "test_")

    def test_merge_logs_to_csv(self, tmp_dir):
        diag = BEASTDiag(os.path.abspath(EXAMPLE_DIR))
        diag.set_burnin(in_percentage=10)
        output_file = os.path.join(tmp_dir, "merged.csv")
        diag.merge_logs_to_csv(output_file=output_file)
        assert os.path.isfile(output_file)
        df = pd.read_csv(output_file)
        assert len(df) > 0
