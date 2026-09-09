"""Unit tests for beast_pype.report_gen.gen_metadata_report."""
import os
import pytest
import tempfile
import shutil
import pandas as pd
import nbformat as nbf

from beast_pype.report_gen import gen_metadata_report


@pytest.fixture
def tmp_dir():
    d = tempfile.mkdtemp()
    yield d
    shutil.rmtree(d)


@pytest.fixture
def metadata_a(tmp_dir):
    path = os.path.join(tmp_dir, "A_metadata.csv")
    pd.DataFrame({
        "strain": [f"a{i}" for i in range(6)],
        "collection date": [
            "2023-01-05", "2023-02-10", "2023-03-15",
            "2023-04-12", "2023-05-01", "2023-06-10",
        ],
    }).to_csv(path, index=False)
    return path


@pytest.fixture
def metadata_b(tmp_dir):
    path = os.path.join(tmp_dir, "B_metadata.tsv")
    pd.DataFrame({
        "strain": [f"b{i}" for i in range(5)],
        "collection date": [
            "2023-03-01", "2023-04-05", "2023-05-14",
            "2023-06-02", "2023-07-09",
        ],
    }).to_csv(path, sep="\t", index=False)
    return path


class TestGenCollectionDateReportSimple:

    def test_creates_notebook(self, tmp_dir, metadata_a):
        out = os.path.join(tmp_dir, "report.ipynb")
        result = gen_metadata_report(
            out, metadata_paths=metadata_a,
            collection_date_field="collection date",
            xml_set_comparisons=False, kernel_name="beast_pype")
        assert result == out
        assert os.path.exists(out)

    def test_notebook_is_valid_and_has_kernelspec(self, tmp_dir, metadata_a):
        out = os.path.join(tmp_dir, "report.ipynb")
        gen_metadata_report(
            out, metadata_paths=metadata_a,
            collection_date_field="collection date",
            xml_set_comparisons=False, kernel_name="my_kernel")
        nb = nbf.read(out, as_version=4)
        nbf.validate(nb)
        assert nb.metadata.kernelspec.name == "my_kernel"

    def test_references_describe_and_histogram(self, tmp_dir, metadata_a):
        out = os.path.join(tmp_dir, "report.ipynb")
        gen_metadata_report(
            out, metadata_paths=metadata_a,
            collection_date_field="collection date",
            xml_set_comparisons=False)
        nb = nbf.read(out, as_version=4)
        source = "\n".join(c.source for c in nb.cells)
        assert "describe_collection_dates(" in source
        assert "plot_collection_date_histogram(" in source
        assert "describe_collection_dates_by_xml_set" not in source
        assert "plot_stacked_collection_date_histogram" not in source

    def test_dict_input_raises(self, tmp_dir, metadata_a, metadata_b):
        out = os.path.join(tmp_dir, "report.ipynb")
        with pytest.raises(ValueError):
            gen_metadata_report(
                out, metadata_paths={"A": metadata_a, "B": metadata_b},
                collection_date_field="collection date",
                xml_set_comparisons=False)


class TestGenCollectionDateReportComparative:

    def test_creates_notebook(self, tmp_dir, metadata_a, metadata_b):
        out = os.path.join(tmp_dir, "report.ipynb")
        result = gen_metadata_report(
            out, metadata_paths={"A": metadata_a, "B": metadata_b},
            collection_date_field="collection date",
            xml_set_comparisons=True, xml_set_label="xml set")
        assert result == out
        assert os.path.exists(out)

    def test_has_all_sequences_and_comparison_sections(
            self, tmp_dir, metadata_a, metadata_b):
        out = os.path.join(tmp_dir, "report.ipynb")
        gen_metadata_report(
            out, metadata_paths={"A": metadata_a, "B": metadata_b},
            collection_date_field="collection date",
            xml_set_comparisons=True, xml_set_label="xml set")
        nb = nbf.read(out, as_version=4)
        nbf.validate(nb)
        source = "\n".join(c.source for c in nb.cells)
        assert "describe_collection_dates(" in source
        assert "plot_collection_date_histogram(" in source
        assert "describe_collection_dates_by_xml_set(" in source
        assert "plot_stacked_collection_date_histogram(" in source
        markdown = "\n".join(
            c.source for c in nb.cells if c.cell_type == "markdown")
        assert "All Sequences" in markdown
        assert "Comparison of XML Sets" in markdown

    def test_string_input_raises(self, tmp_dir, metadata_a):
        out = os.path.join(tmp_dir, "report.ipynb")
        with pytest.raises(ValueError):
            gen_metadata_report(
                out, metadata_paths=metadata_a,
                collection_date_field="collection date",
                xml_set_comparisons=True)
