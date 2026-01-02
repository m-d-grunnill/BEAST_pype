from click.testing import CliRunner
from beast_pype.cli import beast_pype
import yaml
from papermill.iorw import read_yaml_file
from pathlib import Path
from papermill import execute_notebook
import os

def _add_tmp_path_to_params_yml_save(tmp_path, parameters_path):
    parameters = read_yaml_file(parameters_path)
    parameters['overall_save_dir'] = f"{tmp_path}/{parameters['overall_save_dir']}"
    new_parameters_path = f"{tmp_path}/parameters.yml"
    with open(new_parameters_path, 'w') as file:
        yaml.safe_dump(parameters, file)
    file.close()
    return new_parameters_path

def _test_running_of_workflow(tmp_path, subtests, workflow, new_parameters_path):
    runner = CliRunner()
    result = runner.invoke(beast_pype, ['run-workflow', workflow, new_parameters_path])
    with subtests.test("Check for error generation."):
        assert result.exit_code == 0
    possible_phase_5s = [path for path in Path(tmp_path).rglob('Phase-5-Diagnosing-Outputs-and-Generate-Report.ipynb')]
    with subtests.test("Check Diagnostic notebook generated."):
        assert len(possible_phase_5s) == 1
    phase_5_path = possible_phase_5s[0]
    os.chdir(str(phase_5_path.parent))
    execute_notebook(
        input_path='Phase-5-Diagnosing-Outputs-and-Generate-Report.ipynb',
        output_path='Phase-5-Diagnosing-Outputs-and-Generate-Report.ipynb')
    with subtests.test("Check Report generated."):
        assert os.path.exists('outputs_and_reports/BEASTPype-Report.ipynb')
    with subtests.test("Check mcc tree was generated."):
        assert os.path.exists('outputs_and_reports/mcc_tree.nexus')