from click.testing import CliRunner
from beast_pype.cli import beast_pype
import yaml
from papermill.iorw import read_yaml_file
from pathlib import Path
from papermill import execute_notebook
import os


def _check_file_was_generated(
        subtests,
        filename,
        directory=None,
        return_path=False):
    if directory is None:
        directory = os.getcwd()
    possible_notebooks = [path for path in Path(directory).rglob(filename)]
    with subtests.test(f"Check {filename} notebook was generated."):
        assert len(possible_notebooks) == 1
    notebook_path = possible_notebooks[0]
    if return_path:
        return notebook_path

def _check_file_was_not_generated(
        subtests,
        filename,
        directory=None):
    if directory is None:
        directory = os.getcwd()
    possible_notebooks = [path for path in Path(directory).rglob(filename)]
    with subtests.test(f"Check {filename} notebook was NOT generated."):
        assert len(possible_notebooks) == 0

class WorkflowVariationTest:

    parameters_path = None
    workflow = None
    variation  = None
    xml_generation_notebook = None
    diagnostic_notebook= None
    kernel_name = 'beast_pype'


    def test_running_of_workflow(self, subtests, tmp_path):
        self.start_working_dir = os.getcwd()
        parameters = read_yaml_file(self.parameters_path)
        os.chdir(tmp_path)
        for param in ['fasta_path', 'metadata_path', 'template_xml_path', 'ready_to_go_xml']:
            if param in parameters:
                parameters[param] = f"{self.start_working_dir}/{parameters[param]}"
        parameters['max_threads'] = 1  # This allows all tests to be run in parallel as it stops beast and IQ tree running parallel.
        self.parameters = parameters
        tmp_parameters_path = "parameters.yml"
        with open(tmp_parameters_path, 'w') as file:
            yaml.safe_dump(self.parameters, file)
        file.close()
        runner = CliRunner()
        result = runner.invoke(beast_pype, ['run-workflow', self.workflow, tmp_parameters_path])
        with subtests.test(f"Check for error generation: {result.exc_info}"):
            assert result.exit_code == 0
        should_be_generated, should_not_be_generated = self.adding_notebooks_to_lists()
        should_be_generated.append('Phase-4-GNU-Parallel-Running-BEAST.ipynb') # Appended here so it is checked last.
        for notebook in should_be_generated:
            _check_file_was_generated(
                subtests=subtests,
                filename=notebook)
        for notebook in should_not_be_generated:
            _check_file_was_not_generated(
                subtests=subtests,
                filename=notebook)
        phase_5_path = _check_file_was_generated(
            subtests=subtests,
            filename=self.diagnostic_notebook,
            return_path=True
        )
        ### Unfortunately when running this section of certain tests from command
        ### line instead of pycharm pytest seems to get stuck running the notebook
        ### (some point after creating outputs_and_reports).
        ### Hence this is commented out.
        # os.chdir(phase_5_path.parent)
        # execute_notebook(
        #     input_path=self.diagnostic_notebook,
        #     output_path=self.diagnostic_notebook)
        # with subtests.test("Check Report generated."):
        #     assert os.path.exists(f'outputs_and_reports/BEAST_pype-Report.ipynb')
        ###
        os.chdir(self.start_working_dir)


class SimpleWorkflowVariationTest(WorkflowVariationTest):
    diagnostic_notebook = 'Phase-5-Diagnosing-Outputs-and-Generate-Report.ipynb'

    def adding_notebooks_to_lists(self):
        should_be_generated = []
        should_not_be_generated = []
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
        return should_be_generated, should_not_be_generated

class ComparativeWorkflowVariationTest(WorkflowVariationTest):
    diagnostic_notebook =  'Phase-5-Diagnosing-XML-sets-and-Generate-Report.ipynb'

    def adding_notebooks_to_lists(self):
        should_be_generated = ['Phase-1-Metadata-and-Sequence-Separation.ipynb']
        should_not_be_generated = []
        xml_set_labels = list(self.parameters['xml_set_definitions'].keys())
        if self.variation == 'full':
            should_be_generated += [
                                       'Phase-2i-IQTree-Building.ipynb',
                                       'Phase-2i-IQTree-Correction.ipynb',
                                   ] + [f"{xml_set_directory}/Phase-2ii-TreeTime-and-Down-Sampling.ipynb"
                                        for xml_set_directory in xml_set_labels]
        else:
            should_not_be_generated += [
                                           'Phase-2i-IQTree-Building.ipynb',
                                           'Phase-2i-IQTree-Correction.ipynb',
                                       ] + [f"{xml_set_directory}/Phase-2ii-TreeTime-and-Down-Sampling.ipynb"
                                            for xml_set_directory in xml_set_labels]
        if self.variation in ['full', 'no initial tree']:
            should_be_generated += [f"{xml_set_directory}/{self.xml_generation_notebook}"
                                    for xml_set_directory in xml_set_labels]
        else:
            should_not_be_generated += [f"{xml_set_directory}/{self.xml_generation_notebook}"
                                        for xml_set_directory in xml_set_labels]
        return should_be_generated, should_not_be_generated