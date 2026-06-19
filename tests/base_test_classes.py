from click.testing import CliRunner
from beast_pype.cli import beast_pype
import yaml
from papermill.iorw import read_yaml_file
from pathlib import Path
import os
import shutil
from datetime import datetime
# from papermill import execute_notebook

def _check_file_was_generated(
        subtests,
        filename,
        directory,
        test_ran_list):
    with subtests.test(f"Check {filename} notebook was generated."):
        file_path = Path(f"{directory}/{filename}")
        test_ran_list.append(file_path.is_file())
        assert test_ran_list[-1]

def _check_file_was_not_generated(
        subtests,
        filename,
        directory,
        test_ran_list):
    with subtests.test(f"Check {filename} notebook was NOT generated."):
        file_path = Path(f"{directory}/{filename}")
        test_ran_list.append(not file_path.is_file())
        assert test_ran_list[-1]

class WorkflowVariationTest:
    name = None
    parameters_path = None
    workflow = None
    xml_generation_notebook = None
    diagnostic_notebook= None
    tree_building_notebooks = None
    downsampled_tree_building_notebooks = None
    beast_running_notebook = 'Phase-4-GNU-Parallel-Running-BEAST.ipynb'
    kernel_name = 'dev_beast_pype'
    initial_tree_building = False
    downsampling = False
    xml_generation = False
    notebooks_in_xml_set_dirs = []
    tree_building_notebooks = [
        'Phase-2i-IQTree-Building.ipynb',
        'Phase-2i-IQTree-Correction.ipynb',
        'Phase-2ii-TreeTime-and-Downsampling.ipynb',
                               ]
    downsampled_tree_building_notebooks = [
        'Phase-2iii-IQTree-Building-with-Downsampled-Data.ipynb',
        'Phase-2iii-IQTree-Correction-with-Downsampled-Data.ipynb',
        'Phase-2iv-TreeTime-with-Downsampled-Data.ipynb'
                               ]
    xml_set_labels = None
    start_timeout = 600


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
        parameters['specific_run_save_dir'] = '2026-02-20_15-28-29' # datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        self.parameters = parameters
        tmp_parameters_path = f"{overall_save_dir}/parameters.yml"
        with open(tmp_parameters_path, 'w') as file:
            yaml.safe_dump(self.parameters, file)
        file.close()
        self.save_dir = f"{parameters['overall_save_dir']}/{parameters['specific_run_save_dir']}"
        test_ran_ok_list = []
        runner = CliRunner()
        run_args = ['run-workflow',
            '-k', self.kernel_name]
        if self.start_timeout is not None:
            run_args += [
                '--start_timeout', self.start_timeout,
                self.workflow,
                tmp_parameters_path]
        else:
            run_args += [
                self.workflow,
                tmp_parameters_path]
        result = runner.invoke(beast_pype, args=run_args)
        with subtests.test(f"Check for error generation: {result.exc_info}"):
            test_ran_ok_list.append(result.exit_code == 0)
            assert test_ran_ok_list[-1]
        should_be_generated = []
        should_not_be_generated = []
        if self.initial_tree_building:
            should_be_generated += self.tree_building_notebooks
        else:
            should_not_be_generated += self.tree_building_notebooks
        if self.downsampling:
            should_be_generated += self.downsampled_tree_building_notebooks
        else:
            should_not_be_generated += self.downsampled_tree_building_notebooks
        if self.xml_generation:
            should_be_generated.append(self.xml_generation_notebook)
        else:
            should_not_be_generated.append(self.xml_generation_notebook)
        should_be_generated += [self.beast_running_notebook, self.diagnostic_notebook]
        for notebook in should_be_generated:
            if notebook in self.notebooks_in_xml_set_dirs:
                for xml_set_directory in self.xml_set_labels:
                    _check_file_was_generated(
                        subtests=subtests,
                        filename=notebook,
                        directory=f'{self.save_dir}/{xml_set_directory}',
                        test_ran_list=test_ran_ok_list)
            else:
                _check_file_was_generated(
                    subtests=subtests,
                    filename=notebook,
                    directory=self.save_dir,
                    test_ran_list=test_ran_ok_list)
        for notebook in should_not_be_generated:
            if notebook in self.notebooks_in_xml_set_dirs:
                for xml_set_directory in self.xml_set_labels:
                    _check_file_was_not_generated(
                        subtests=subtests,
                        filename=notebook,
                        directory=f'{self.save_dir}/{xml_set_directory}',
                        test_ran_list=test_ran_ok_list)
            else:
                _check_file_was_not_generated(
                    subtests=subtests,
                    filename=notebook,
                    directory=self.save_dir,
                    test_ran_list=test_ran_ok_list)
        ### Unfortunately when running this section of certain tests from command
        ### line instead of pycharm pytest seems to get stuck running the notebook
        ### (some point after creating outputs_and_reports).
        ### Hence this is commented out.
        # phase_5_path = f'{self.name}/{self.diagnostic_notebook}"
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




class ComparativeWorkflowVariationTest(WorkflowVariationTest):
    diagnostic_notebook =  'Phase-5-Diagnosing-XML-sets-and-Generate-Report.ipynb'

    notebooks_in_xml_set_dirs = ['Phase-2i-IQTree-Correction.ipynb',
                                 'Phase-2ii-TreeTime-and-Downsampling.ipynb',
                                 'Phase-2iii-IQTree-Correction-with-Downsampled-Data.ipynb',
                                 'Phase-2iv-TreeTime-with-Downsampled-Data.ipynb',
                                 'Phase-3-Gen-BDSKY-Serial-xml.ipynb',
                                 'Phase-3-Gen-Generic-xml.ipynb']

    @property
    def xml_set_labels(self):
        return list(self.parameters['xml_set_definitions'].keys())
