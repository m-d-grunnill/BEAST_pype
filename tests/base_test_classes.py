from click.testing import CliRunner
from beast_pype.cli import beast_pype
import yaml
from papermill.iorw import read_yaml_file
from pathlib import Path
from papermill import execute_notebook
import os


def _check_file_was_generated(
        subtests,
        directory,
        filename,
        return_path=False):
    possible_notebooks = [path for path in Path(directory).rglob(filename)]
    with subtests.test(f"Check {filename} notebook was generated."):
        assert len(possible_notebooks) == 1
    notebook_path = possible_notebooks[0]
    if return_path:
        return notebook_path

def _check_file_was_not_generated(
        subtests,
        directory,
        filename):
    possible_notebooks = [path for path in Path(directory).rglob(filename)]
    with subtests.test(f"Check {filename} notebook was generated."):
        assert len(possible_notebooks) == 0

class WorkflowVariationTest:

    parameters_path = None
    workflow = None
    variation  = None
    xml_generation_notebook = None
    diagnostic_notebook= None

    def test_running_of_workflow(self, subtests, tmp_path):
        parameters = read_yaml_file(self.parameters_path)
        parameters['overall_save_dir'] = f"{tmp_path}/{parameters['overall_save_dir']}"
        tmp_parameters_path = f"{tmp_path}/parameters.yml"
        with open(tmp_parameters_path, 'w') as file:
            yaml.safe_dump(parameters, file)
        file.close()
        runner = CliRunner()
        result = runner.invoke(beast_pype, ['run-workflow', self.workflow, tmp_parameters_path])
        with subtests.test("Check for error generation."):
            assert result.exit_code == 0
        should_not_be_generated = []
        should_be_generated = ['Phase-4-GNU-Parallel-Running-BEAST.ipynb']
        self.adding_notebooks_to_lists(should_be_generated, should_not_be_generated,
                                       subtests, tmp_path)
        for notebook in should_be_generated:
            _check_file_was_generated(
                subtests=subtests,
                directory=tmp_path,
                filename=notebook)
        for notebook in should_not_be_generated:
            _check_file_was_not_generated(
                subtests=subtests,
                directory=tmp_path,
                filename=notebook)
        phase_5_path = _check_file_was_generated(
            subtests=subtests,
            directory=tmp_path,
            filename=self.diagnostic_notebook,
            return_path=True
        )
        start_working_dir = os.getcwd()
        os.chdir(phase_5_path.parent)
        execute_notebook(
            input_path=self.diagnostic_notebook,
            output_path=self.diagnostic_notebook)
        with subtests.test("Check Report generated."):
            assert os.path.exists(f'outputs_and_reports/BEASTPype-Report.ipynb')

        os.chdir(start_working_dir)


class SimpleWorkflowVariationTest(WorkflowVariationTest):
    diagnostic_notebook = 'Phase-5-Diagnosing-Outputs-and-Generate-Report.ipynb'

    def adding_notebooks_to_lists(self, should_be_generated, should_not_be_generated, subtests, tmp_path):
        if self.variation == 'full':
            should_be_generated += [
                'Phase-2i-IQTree-Building.ipynb',
                'Phase-2i-IQTree-Correction.ipynb',
                'Phase-2ii-TreeTime-and-Down-Sampling.ipynb'
            ]
        else:
            should_not_be_generated += [
                'Phase-2i-IQTree-Building.ipynb',
                'Phase-2i-IQTree-Correction.ipynb',
                'Phase-2ii-TreeTime-and-Down-Sampling.ipynb'
            ]
        if self.variation in ['full', 'no initial tree']:
            should_be_generated += [self.xml_generation_notebook]
        else:
            should_not_be_generated += [self.xml_generation_notebook]

class ComparativeWorkflowVariationTest(WorkflowVariationTest):
    diagnostic_notebook =  'Phase-5-Diagnosing-XML-sets-and-Generate-Report.ipynb'

    def adding_notebooks_to_lists(self, should_be_generated, should_not_be_generated, subtests, tmp_path):
        phase_1_path = _check_file_was_generated(
            subtests=subtests,
            directory=tmp_path,
            filename='Phase-1-Metadata-and-Sequence-Separation.ipynb',
            return_path=True)
        xml_set_directories = [f.path for f in os.scandir(phase_1_path.parent) if f.is_dir() and not f.name.startswith('.')]
        if self.variation == 'full':
            should_be_generated += [
                'Phase-2i-IQTree-Building.ipynb',
                'Phase-2i-IQTree-Correction.ipynb',
            ] + [f"{xml_set_directory}/Phase-2ii-TreeTime-and-Down-Sampling.ipynb"
                 for xml_set_directory in xml_set_directories]
        else:
            should_not_be_generated += [
                                       'Phase-2i-IQTree-Building.ipynb',
                                       'Phase-2i-IQTree-Correction.ipynb',
                                   ] + [f"{xml_set_directory}/Phase-2ii-TreeTime-and-Down-Sampling.ipynb"
                                        for xml_set_directory in xml_set_directories]
        if self.variation in ['full', 'no initial tree']:
                should_be_generated += [f"{xml_set_directory}/{self.xml_generation_notebook}"
                                        for xml_set_directory in xml_set_directories]
        else:
                should_not_be_generated += [f"{xml_set_directory}/{self.xml_generation_notebook}"
                                        for xml_set_directory in xml_set_directories]
