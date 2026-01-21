import click
import importlib.resources as importlib_resources
import os
from beast_pype.nb_utils import execute_notebook
from .diagnostics import gen_beast_diagnostic_nb
from datetime import datetime
from papermill.iorw import read_yaml_file

workflows_path = importlib_resources.path('beast_pype', 'workflows')
available_workflows = [file for file in os.listdir(workflows_path) if file.endswith('.ipynb')]
default_workflow_names = [
    file.replace('.ipynb', '')
    for file in available_workflows]

reports_to_exclude = ['COVID-Strain-Surveillance.ipynb']
reports_path = importlib_resources.path('beast_pype', 'report_templates')
available_reports = [file for file in os.listdir(reports_path)
                     if file.endswith('.ipynb') and file not in reports_to_exclude]
default_report_names = [
    file.replace('.ipynb', '')
    for file in available_reports]

diag_valid_params = [
    'metadata_path',
    'rt_partitions',
    'sampling_prop_partition_freq',
    'collection_date_field'
]

def _is_int(value):
    """Use casting to check if value can convert to an `int`."""
    try:
        int(value)
    except ValueError:
        return False
    else:
        return True


def _is_float(value):
    """Use casting to check if value can convert to a `float`."""
    try:
        float(value)
    except ValueError:
        return False
    else:
        return True

def _resolve_type(value):
    if value == "True":
        return True
    elif value == "False":
        return False
    elif value in ["None", "null"]:
        return None
    elif _is_int(value):
        return int(value)
    elif _is_float(value):
        return float(value)
    else:
        return value



@click.group(context_settings=dict(help_option_names=['-h', '--help']),
             epilog='See https://github.com/m-d-grunnill/BEAST_pype/wiki for further documentation.')
def beast_pype():
    pass


@beast_pype.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('workflow', required=1, type=click.Choice(default_workflow_names))
@click.argument('parameters', required=1, type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True))
@click.option('--kernel_name', '-k', default='beast_pype', type=str,
              help='Name of Jupyter python kernel_name to use when running workflow.\n' +
                   'This is also the name of the conda environment to use in phases 4 & ' +
                   'phase 2ii (as these Jupyter notebooks use the `bash` kernel_name).\n' +
                   'If not given "beast_pype" is used.'
              )
def run_workflow(workflow,
                 parameters, kernel_name):
    """
    WORKFLOW: Workflow you wish to execute, or a path to jupyter notebook to be run  as a workflow.\n
    PARAMETERS: Path to YAML file containing parameters.
    """
    parameters = read_yaml_file(parameters)
    parameters['kernel_name'] = kernel_name
    workflow_save_name = f'{workflow}.ipynb'
    workflow = f'{workflows_path}/{workflow_save_name}'

    os.makedirs(parameters['overall_save_dir'], exist_ok=True)
    if 'specific_run_save_dir' not in parameters:
        parameters['specific_run_save_dir'] = datetime.now().strftime(
            "%Y-%m-%d_%H-%M-%S")
    save_dir = f"{parameters['overall_save_dir']}/{parameters['specific_run_save_dir']}"
    os.makedirs(save_dir)
    execute_notebook(
        input_path=workflow,
        output_path=f"{save_dir}/{workflow_save_name}",
        parameters=parameters,
        kernel_name=kernel_name,
        progress_bar=True,
        nest_asyncio=True
    )

@beast_pype.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('report_template', required=1,
                type=click.Choice(default_report_names))
@click.argument('beast_outputs',
                 required=1,
                 type=click.Path(exists=True, dir_okay=True, file_okay=False, readable=True, writable=True))
@click.option('--parameters', '-p', nargs=2, multiple=True, help='Parameters to be passed to workflow.')
@click.option('--parameters_file', '-f', multiple=True, help='Path to YAML file containing parameters.')
@click.option('--kernel_name', '-k', default='beast_pype', type=str,
              help='Name of Jupyter python kernel_name to use when running diagnostic & report template notebooks.\n' +
                   'If not given "beast_pype" is used.'
              )
def diagnose_results(report_template,
                     beast_outputs,
                     parameters,
                     parameters_file,
                     kernel_name):
    """
    REPORT_TEMPLATE: Report template to use after diagnosing BEAST 2 outputs.\n
    BEAST_OUTPUTS: Path to directory containing BEAST 2 outputs to diagnose.
    """
    # Read in Parameters
    parameters_final = {}
    for name, value in parameters or []:
        parameters_final[name] = _resolve_type(value)
    for files in parameters_file or []:
        parameters_final.update(read_yaml_file(files) or {})

    if 'beast_xml_path' in parameters_final:
        beast_xml_path = parameters_final.pop('beast_xml_path')
    else:
        beast_xml_path = None
    for param, value in parameters_final.items():
        if param not in diag_valid_params:
            raise ValueError(f"Parameter {param} is not a valid parameter for use in the diagnostic workflow.")
    if report_template is None:
        raise ValueError('A value for --report_template is required for the diagnostic workflow.')
    if beast_outputs is None:
        raise ValueError('A value for --beast_outputs is required for the diagnostic workflow.')

    gen_beast_diagnostic_nb(
        beast_outputs=beast_outputs,
        report_template=report_template,
        beast_xml_path=beast_xml_path,
        kernel_name=kernel_name,
        **parameters_final)
