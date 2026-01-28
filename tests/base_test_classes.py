from click.testing import CliRunner
from beast_pype.cli import beast_pype
import yaml
from papermill.iorw import read_yaml_file
from pathlib import Path
from papermill import execute_notebook
import os
import shutil
from copy import deepcopy


def _check_file_was_generated(
        subtests,
        filename,
        directory,
        test_ran_list,
        return_path=False):
    possible_notebooks = [path for path in Path(directory).rglob(filename)]
    with subtests.test(f"Check {filename} notebook was generated."):
        test_ran_list.append(len(possible_notebooks) == 1)
        assert test_ran_list[-1]
    notebook_path = possible_notebooks[0]
    if return_path:
        return notebook_path

def _check_file_was_not_generated(
        subtests,
        filename,
        directory,
        test_ran_list):
    possible_notebooks = [path for path in Path(directory).rglob(filename)]
    with subtests.test(f"Check {filename} notebook was NOT generated."):
        test_ran_list.append(len(possible_notebooks) == 0)
        assert test_ran_list[-1]

class WorkflowVariationTest:
    name = None
    parameters_path = None
    workflow = None
    variation  = None
    xml_generation_notebooks = None
    xml_generation_notebook = None
    diagnostic_notebook= None
    tree_building_notebooks = None
    downsampled_tree_building_notebooks = None
    beast_running_notebook = 'Phase-4-GNU-Parallel-Running-BEAST.ipynb'
    kernel_name = 'dev_beast_pype'


    def test_running_of_workflow(self, subtests, tmp_path):
        self.start_working_dir = os.getcwd()
        parameters = read_yaml_file(self.parameters_path)
        parameters['overall_save_dir'] = self.name
        overall_save_dir = parameters['overall_save_dir']
        os.makedirs(overall_save_dir, exist_ok=True)
        for param in ['fasta_path', 'metadata_path', 'template_xml_path', 'ready_to_go_xml']:
            if param in parameters:
                parameters[param] = f"{self.start_working_dir}/{parameters[param]}"
        parameters['max_threads'] = 1  # This allows all tests to be run in parallel as it stops beast and IQ tree running parallel.
        self.parameters = parameters
        tmp_parameters_path = f"{overall_save_dir}/parameters.yml"
        with open(tmp_parameters_path, 'w') as file:
            yaml.safe_dump(self.parameters, file)
        file.close()
        test_ran_ok_list = []
        runner = CliRunner()
        result = runner.invoke(beast_pype, ['run-workflow','-k', self.kernel_name, self.workflow, tmp_parameters_path])
        with subtests.test(f"Check for error generation: {result.exc_info}"):
            test_ran_ok_list.append(result.exit_code == 0)
            assert test_ran_ok_list[-1]
        should_be_generated = []
        should_not_be_generated = []
        if self.variation in ['initial tree building', 'downsampling'] :
            should_be_generated += self.tree_building_notebooks
        else:
            should_not_be_generated += self.tree_building_notebooks
        if self.variation == 'downsampling':
            should_be_generated += self.downsampled_tree_building_notebooks
        else:
            should_not_be_generated += self.downsampled_tree_building_notebooks
        if self.variation != 'xml ready-to-go':
            should_be_generated += self.xml_generation_notebooks
        else:
            should_not_be_generated += self.xml_generation_notebooks
        should_be_generated.append(self.beast_running_notebook) # Appended here so it is checked last.
        for notebook in should_be_generated:
            _check_file_was_generated(
                subtests=subtests,
                filename=notebook,
                directory=overall_save_dir,
                test_ran_list=test_ran_ok_list,)
        for notebook in should_not_be_generated:
            _check_file_was_not_generated(
                subtests=subtests,
                filename=notebook,
                directory=overall_save_dir,
                test_ran_list=test_ran_ok_list)
        phase_5_path = _check_file_was_generated(
            subtests=subtests,
            filename=self.diagnostic_notebook,
            directory=overall_save_dir,
            test_ran_list=test_ran_ok_list,
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
        # os.chdir(self.start_working_dir)
        ###
        if all(test_ran_ok_list):
            shutil.rmtree(overall_save_dir)


class SimpleWorkflowVariationTest(WorkflowVariationTest):
    diagnostic_notebook = 'Phase-5-Diagnosing-Outputs-and-Generate-Report.ipynb'
    tree_building_notebooks = [
        'Phase-2i-IQTree-Building.ipynb', 'Phase-2i-IQTree-Correction.ipynb', 'Phase-2ii-TreeTime-and-Downsampling.ipynb',
                               ]
    downsampled_tree_building_notebooks = [
        'Phase-2iii-IQTree-Building.ipynb', 'Phase-2iii-IQTree-Correction.ipynb', 'Phase-2iv-TreeTime-with-Downsampled-Data.ipynb'
                               ]

    @property
    def xml_generation_notebooks(self):
        return [self.xml_generation_notebook]



class ComparativeWorkflowVariationTest(WorkflowVariationTest):
    diagnostic_notebook =  'Phase-5-Diagnosing-XML-sets-and-Generate-Report.ipynb'

    @property
    def xml_set_labels(self):
        return list(self.parameters['xml_set_definitions'].keys())

    @property
    def tree_building_notebooks(self):
        notebook_list = [
                            'Phase-2i-IQTree-Building.ipynb', 'Phase-2i-IQTree-Correction.ipynb']
        notebook_list += [f"{xml_set_directory}/Phase-2iv-TreeTime-with-Downsampled-Data.ipynb"
                          for xml_set_directory in self.xml_set_labels]
        return notebook_list

    @property
    def downsampled_tree_building_notebooks(self):
        notebook_list = ['Phase-2iii-IQTree-Building.ipynb', 'Phase-2iii-IQTree-Correction.ipynb']
        notebook_list += [f"{xml_set_directory}/Phase-2iv-TreeTime-with-Downsampled-Data.ipynb"
                          for xml_set_directory in self.xml_set_labels]
        return notebook_list

    @property
    def xml_generation_notebooks(self):
        return [f"{xml_set_directory}/{self.xml_generation_notebook}"
                for xml_set_directory in self.xml_set_labels]