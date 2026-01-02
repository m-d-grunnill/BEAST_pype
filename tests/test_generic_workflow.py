from click.testing import CliRunner
from beast_pype.cli import beast_pype
import yaml
from papermill.iorw import read_yaml_file
from pathlib import Path
from papermill import execute_notebook
import os
from _test_parts import WorkFlowVariationTest

class TestGenericFull(WorkFlowVariationTest):
    parameter_path = 'parameters/locally_run_examples/Generic_full.yml'
    workflow = 'Generic'
    variation = 'full'