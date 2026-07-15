"""Unit tests for beast_pype.diagnostics.static_diag module."""
import os
import pytest
import tempfile
import shutil
import numpy as np
import nbformat as nbf
from unittest.mock import patch

from beast_pype.diagnostics.static_diag import gen_static_diagnostic_nb


# ============================================================================
# Fixtures
# ============================================================================

EXAMPLE_DIR = os.path.join(
    os.path.dirname(__file__),
    "..",
    "example_beast2_outputs_to_diagnose",
)


@pytest.fixture
def tmp_dir():
    """Create a temporary directory and clean up after test."""
    d = tempfile.mkdtemp()
    yield d
    shutil.rmtree(d)


@pytest.fixture
def synthetic_beast_dir(tmp_dir):
    """Create a directory with synthetic .log and .trees files."""
    n_samples = 50
    step_size = 1000

    for seed in [111, 222]:
        # Create log files
        filepath = os.path.join(tmp_dir, f"run-with-seed-{seed}-BEAST.log")
        samples = np.arange(0, n_samples * step_size, step_size)
        lines = ["# comment line\n"]
        header = "Sample\tposterior\tlikelihood\ttreeHeight\n"
        lines.append(header)
        for s in samples:
            vals = [str(s),
                    str(np.random.randn()),
                    str(np.random.randn()),
                    str(abs(np.random.randn()) + 1)]
            lines.append("\t".join(vals) + "\n")
        with open(filepath, "w") as f:
            f.writelines(lines)

        # Create .trees files
        trees_path = os.path.join(tmp_dir, f"run-with-seed-{seed}-BEAST.trees")
        tree_lines = [
            "#NEXUS\n",
            "Begin taxa;\n",
            "    Dimensions ntax=3;\n",
            "    Taxlabels\n",
            "        A\n",
            "        B\n",
            "        C\n",
            "    ;\n",
            "End;\n",
            "Begin trees;\n",
        ]
        for s in samples:
            tree_lines.append(
                f"    tree STATE_{s} = ((A:1,B:1):1,C:2);\n"
            )
        tree_lines.append("End;\n")
        with open(trees_path, "w") as f:
            f.writelines(tree_lines)

    return tmp_dir


# ============================================================================
# Tests: Input validation
# ============================================================================


class TestGenStaticDiagnosticNbValidation:
    """Tests for input validation in gen_static_diagnostic_nb."""

    def test_directory_not_found(self, tmp_dir):
        fake_dir = os.path.join(tmp_dir, "nonexistent")
        with pytest.raises(FileNotFoundError, match="Directory not found"):
            gen_static_diagnostic_nb(fake_dir)

    def test_burnin_negative(self, synthetic_beast_dir):
        with pytest.raises(ValueError, match="burnin must be between 0 and 100"):
            gen_static_diagnostic_nb(synthetic_beast_dir, burnin=-5)

    def test_burnin_100(self, synthetic_beast_dir):
        with pytest.raises(ValueError, match="burnin must be between 0 and 100"):
            gen_static_diagnostic_nb(synthetic_beast_dir, burnin=100)

    def test_no_log_files(self, tmp_dir):
        with pytest.raises(FileNotFoundError, match="No .log files found"):
            gen_static_diagnostic_nb(tmp_dir)


# ============================================================================
# Tests: Notebook generation (mocking execution)
# ============================================================================


class TestGenStaticDiagnosticNbOutput:
    """Tests for gen_static_diagnostic_nb output (mocking execution and HTML export)."""

    @patch("beast_pype.diagnostics.static_diag.execute_notebook")
    def test_returns_dict_with_expected_keys(self, mock_execute, synthetic_beast_dir, tmp_dir):
        output_prefix = os.path.join(tmp_dir, "test_output_")
        mock_execute.return_value = None

        result = gen_static_diagnostic_nb(
            synthetic_beast_dir,
            burnin=10,
            output_prefix=output_prefix,
        )
        assert "notebook" in result
        assert "notebook_html" in result
        assert "merged_log" in result
        assert "merged_trees" in result

    @patch("beast_pype.diagnostics.static_diag.execute_notebook")
    def test_notebook_file_created(self, mock_execute, synthetic_beast_dir, tmp_dir):
        output_prefix = os.path.join(tmp_dir, "test_output_")
        mock_execute.return_value = None

        result = gen_static_diagnostic_nb(
            synthetic_beast_dir,
            burnin=10,
            output_prefix=output_prefix,
        )
        assert os.path.isfile(result["notebook"])

    @patch("beast_pype.diagnostics.static_diag.execute_notebook")
    def test_merged_log_csv_created(self, mock_execute, synthetic_beast_dir, tmp_dir):
        output_prefix = os.path.join(tmp_dir, "test_output_")
        mock_execute.return_value = None

        result = gen_static_diagnostic_nb(
            synthetic_beast_dir,
            burnin=10,
            output_prefix=output_prefix,
        )
        assert os.path.isfile(result["merged_log"])
        assert result["merged_log"].endswith(".csv")

    @patch("beast_pype.diagnostics.static_diag.execute_notebook")
    def test_merged_trees_created(self, mock_execute, synthetic_beast_dir, tmp_dir):
        output_prefix = os.path.join(tmp_dir, "test_output_")
        mock_execute.return_value = None

        result = gen_static_diagnostic_nb(
            synthetic_beast_dir,
            burnin=10,
            output_prefix=output_prefix,
        )
        assert os.path.isfile(result["merged_trees"])
        assert result["merged_trees"].endswith(".trees")

    @patch("beast_pype.diagnostics.static_diag.execute_notebook")
    def test_notebook_has_correct_structure(self, mock_execute, synthetic_beast_dir, tmp_dir):
        output_prefix = os.path.join(tmp_dir, "test_output_")
        mock_execute.return_value = None

        result = gen_static_diagnostic_nb(
            synthetic_beast_dir,
            burnin=10,
            output_prefix=output_prefix,
            parameters_per_section=2,
        )
        nb = nbf.read(result["notebook"], as_version=4)
        cell_types = [c["cell_type"] for c in nb["cells"]]
        assert "markdown" in cell_types
        assert "code" in cell_types

    @patch("beast_pype.diagnostics.static_diag.execute_notebook")
    def test_notebook_contains_parameter_sections(self, mock_execute, synthetic_beast_dir, tmp_dir):
        output_prefix = os.path.join(tmp_dir, "test_output_")
        mock_execute.return_value = None

        result = gen_static_diagnostic_nb(
            synthetic_beast_dir,
            burnin=10,
            output_prefix=output_prefix,
            parameters_per_section=1,
        )
        nb = nbf.read(result["notebook"], as_version=4)
        markdown_cells = [c for c in nb["cells"] if c["cell_type"] == "markdown"]
        param_headings = [c for c in markdown_cells if "## Parameters:" in c["source"]]
        assert len(param_headings) >= 1

    @patch("beast_pype.diagnostics.static_diag.execute_notebook")
    def test_notebook_has_kernelspec(self, mock_execute, synthetic_beast_dir, tmp_dir):
        output_prefix = os.path.join(tmp_dir, "test_output_")
        mock_execute.return_value = None

        result = gen_static_diagnostic_nb(
            synthetic_beast_dir,
            burnin=10,
            output_prefix=output_prefix,
            kernel_name="test_kernel",
        )
        nb = nbf.read(result["notebook"], as_version=4)
        assert "kernelspec" in nb["metadata"]
        assert nb["metadata"]["kernelspec"]["name"] == "test_kernel"

    @patch("beast_pype.diagnostics.static_diag.execute_notebook")
    def test_default_output_prefix(self, mock_execute, synthetic_beast_dir):
        mock_execute.return_value = None

        result = gen_static_diagnostic_nb(
            synthetic_beast_dir,
            burnin=10,
        )
        expected_prefix = os.path.join(synthetic_beast_dir, "static_diag_")
        assert result["notebook"].startswith(expected_prefix)

    @patch("beast_pype.diagnostics.static_diag.execute_notebook")
    def test_burnin_zero(self, mock_execute, synthetic_beast_dir, tmp_dir):
        """burnin=0 should work (no samples removed)."""
        output_prefix = os.path.join(tmp_dir, "test_output_")
        mock_execute.return_value = None

        result = gen_static_diagnostic_nb(
            synthetic_beast_dir,
            burnin=0,
            output_prefix=output_prefix,
        )
        assert os.path.isfile(result["notebook"])

    @patch("beast_pype.diagnostics.static_diag.execute_notebook")
    def test_no_trees_files(self, mock_execute, tmp_dir):
        """If no .trees files, merged_trees should be None."""
        n_samples = 20
        step_size = 1000
        filepath = os.path.join(tmp_dir, "run1-BEAST.log")
        lines = ["# comment\n", "Sample\tposterior\tlikelihood\n"]
        for s in range(0, n_samples * step_size, step_size):
            lines.append(f"{s}\t{np.random.randn()}\t{np.random.randn()}\n")
        with open(filepath, "w") as f:
            f.writelines(lines)

        output_prefix = os.path.join(tmp_dir, "output_")
        mock_execute.return_value = None

        result = gen_static_diagnostic_nb(
            tmp_dir,
            burnin=10,
            output_prefix=output_prefix,
        )
        assert result["merged_trees"] is None

    @patch("beast_pype.diagnostics.static_diag.execute_notebook")
    def test_execute_notebook_called_with_kernel_name(self, mock_execute, synthetic_beast_dir, tmp_dir):
        output_prefix = os.path.join(tmp_dir, "test_output_")
        mock_execute.return_value = None

        gen_static_diagnostic_nb(
            synthetic_beast_dir,
            burnin=10,
            output_prefix=output_prefix,
            kernel_name="my_custom_kernel",
        )
        mock_execute.assert_called_once()
        call_kwargs = mock_execute.call_args[1]
        assert call_kwargs["kernel_name"] == "my_custom_kernel"
