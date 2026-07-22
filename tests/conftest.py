import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--kernel-name",
        action="store",
        default=None,
        help="Jupyter kernel name to use in workflow tests (overrides WorkflowVariationTest.kernel_name)."
    )
