from click.testing import CliRunner
from beast_pype.cli import beast_pype
import yaml
from papermill.iorw import read_yaml_file
from pathlib import Path
from papermill import execute_notebook
import os
from _test_parts import _add_tmp_path_to_params_yml_save, _test_running_of_workflow



def test_full_workflow(tmp_path, subtests):
    new_parameters_path = _add_tmp_path_to_params_yml_save(
        tmp_path,
        'parameters/locally_run_examples/Generic_full.yml')
    worklflow = 'Generic'
    _test_running_of_workflow(tmp_path, subtests, worklflow, new_parameters_path)