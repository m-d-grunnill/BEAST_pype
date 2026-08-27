import click
import os
import nbformat
from nbconvert import HTMLExporter
from beast_pype.nb_utils import execute_notebook
from beast_pype.path_utils import path_to_workflows, path_to_report_templates
from beast_pype.diagnostics import gen_beast_diagnostic_nb, gen_static_diagnostic_nb
from beast_pype.report_gen import gen_summary_tree_report, gen_parameters_report, gen_metadata_report
from datetime import datetime
from papermill.iorw import read_yaml_file


workflows_path = path_to_workflows()
available_workflows = [file.name for file in workflows_path.glob('*.ipynb')]
default_workflow_names = [
    file.replace('.ipynb', '')
    for file in available_workflows]

reports_to_exclude = ['COVID-Strain-Surveillance.ipynb']
reports_path = path_to_report_templates()
available_reports = [file.name for file in reports_path.glob('*.ipynb')
                     if not file.name in reports_to_exclude]
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
@click.option('--start_timeout', '-t', default=60, type=int,
              help='Time in seconds to wait for the kernel to start before raising an error when launching a jupyter notebook.\n' +
                   'If not given 60 seconds (1 minutes) is used.'
              )
def run_workflow(workflow,
                 parameters,
                 kernel_name,
                 start_timeout):
    """
    WORKFLOW: Workflow you wish to execute, or a path to jupyter notebook to be run  as a workflow.\n
    PARAMETERS: Path to YAML file containing parameters.
    """
    parameters = read_yaml_file(parameters)
    parameters['kernel_name'] = kernel_name
    workflow_save_name = f'{workflow}.ipynb'
    workflow = f'{workflows_path}/{workflow_save_name}'
    if 'specific_run_save_dir' not in parameters:
        parameters['specific_run_save_dir'] = datetime.now().strftime(
            "%Y-%m-%d_%H-%M-%S")
    os.makedirs(parameters['overall_save_dir'], exist_ok=True)
    save_dir = f"{parameters['overall_save_dir']}/{parameters['specific_run_save_dir']}"
    os.makedirs(save_dir)
    parameters['start_timeout'] = start_timeout
    execute_notebook(
        input_path=workflow,
        output_path=f"{save_dir}/{workflow_save_name}",
        parameters=parameters,
        kernel_name=kernel_name,
        progress_bar=True,
        nest_asyncio=True, start_timeout=start_timeout
    )

@beast_pype.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('beast_outputs',
                 required=1,
                 type=click.Path(exists=True, dir_okay=True, file_okay=False, readable=True, writable=True))
@click.option('--burnin', '-b', default=10, type=float,
              help='Burn-in percentage (0-100) to remove from the start of each chain.\n' +
                   'If not given 10%% is used.'
              )
@click.option('--output_prefix', '-o', default=None, type=str,
              help='Prefix (including path) for output files.\n' +
                   'If not given, defaults to <beast_outputs>/static_diag_.'
              )
@click.option('--parameters_per_section', '-n', default=1, type=int,
              help='Number of parameters to display per notebook section.\n' +
                   'If not given 1 is used.'
              )
@click.option('--kernel_name', '-k', default='beast_pype', type=str,
              help='Name of the Jupyter kernel to use when executing the diagnostic notebook.\n' +
                   'If not given "beast_pype" is used.'
              )
def static_diagnose_and_merge(beast_outputs,
                              burnin,
                              output_prefix,
                              parameters_per_section,
                              kernel_name):
    """
    BEAST_OUTPUTS: Path to directory containing BEAST 2 outputs to statically diagnose and merge.
    """
    results = gen_static_diagnostic_nb(
        directory=beast_outputs,
        burnin=burnin,
        output_prefix=output_prefix,
        parameters_per_section=parameters_per_section,
        kernel_name=kernel_name,
    )
    click.echo(f"Notebook: {results['notebook']}")
    click.echo(f"Notebook HTML: {results['notebook_html']}")
    click.echo(f"Merged log: {results['merged_log']}")
    if results['merged_trees']:
        click.echo(f"Merged trees: {results['merged_trees']}")
    else:
        click.echo("No .trees files found; merged trees not generated.")


@beast_pype.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--summary_tree', '-s', multiple=True, required=True,
              type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
              help='Path to a summary tree file (nexus format). Can be specified multiple times.')
@click.option('--beast_xml', '-x', multiple=True, required=True,
              type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
              help='Path to a BEAST 2 XML file corresponding to a summary tree. '
                   'Must be specified the same number of times as --summary_tree, in the same order.')
@click.option('--output', '-o', required=True, type=str,
              help='Path to save the output report notebook (.ipynb).')
@click.option('--collection_date_field', '-c', default='collection_date', type=str,
              help='Name of field in metadata containing collection dates (YYYY-MM-DD). '
                   'Default: "collection_date".')
@click.option('--plot_width', default=12, type=float,
              help='Width of tree plots in inches. Default: 12.')
@click.option('--plot_height', default=6, type=float,
              help='Height of tree plots in inches. Default: 6.')
@click.option('--plot_res', default=120, type=int,
              help='PPI for rasterization. Default: 120.')
@click.option('--kernel_name', '-k', default='beast_pype_R', type=str,
              help='Name of the Jupyter kernel to use when executing the notebook. '
                   'Default: "beast_pype_R".')
def summary_tree_report(summary_tree,
                        beast_xml,
                        output,
                        collection_date_field,
                        plot_width,
                        plot_height,
                        plot_res,
                        kernel_name):
    """Generate a report notebook summarizing BEAST 2 summary trees."""
    summary_trees = list(summary_tree) if len(summary_tree) > 1 else summary_tree[0]
    beast_xml_paths = list(beast_xml) if len(beast_xml) > 1 else beast_xml[0]
    results = gen_summary_tree_report(
        output_report_path=output,
        summary_trees=summary_trees,
        beast_xml_paths=beast_xml_paths,
        collection_date_field=collection_date_field,
        plot_width=plot_width,
        plot_height=plot_height,
        plot_res=plot_res,
        kernel_name=kernel_name,
    )
    click.echo(f"Notebook: {results['notebook']}")
    click.echo(f"Notebook HTML: {results['notebook_html']}")



@beast_pype.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--metadata', '-m', multiple=True, required=True,
              type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
              help='Path to a metadata file (.csv or .tsv). Can be specified multiple times. '
                   'When more than one is provided, comparative (xml-set) visualisations are used.')
@click.option('--xml_set_name', '-n', multiple=True, type=str,
              help='Name for an xml set, paired with each --metadata in order. '
                   'Only used (and required) when more than one --metadata is provided.')
@click.option('--output', '-o', required=True, type=str,
              help='Path to save the output report notebook (.ipynb).')
@click.option('--collection_date_field', '-c', default='collection_date', type=str,
              help='Name of field in metadata containing collection dates (YYYY-MM-DD). '
                   'Default: "collection_date".')
@click.option('--xml_set_label', default='xml set', type=str,
              help='Label for the xml set grouping variable in comparative reports. '
                   'Default: "xml set".')
@click.option('--kernel_name', '-k', default='beast_pype', type=str,
              help='Name of the Jupyter kernel to use when executing the notebook. '
                   'Default: "beast_pype".')
def metadata_report(metadata,
                    xml_set_name,
                    output,
                    collection_date_field,
                    xml_set_label,
                    kernel_name):
    """Generate a report notebook summarising sample metadata.

    Summarises sample collection date metadata via a describe table and a
    histogram. When more than one --metadata is given, an additional section
    compares each xml set.
    """
    xml_set_comparisons = len(metadata) > 1
    if xml_set_comparisons:
        if len(xml_set_name) != len(metadata):
            raise click.UsageError(
                'When more than one --metadata is provided, a matching '
                '--xml_set_name must be given for each (same number, same order).')
        metadata_paths = dict(zip(xml_set_name, metadata))
    else:
        metadata_paths = metadata[0]

    notebook_path = gen_metadata_report(
        output_report_path=output,
        metadata_paths=metadata_paths,
        collection_date_field=collection_date_field,
        xml_set_comparisons=xml_set_comparisons,
        xml_set_label=xml_set_label,
        kernel_name=kernel_name,
    )
    execute_notebook(
        input_path=notebook_path,
        output_path=notebook_path,
        kernel_name=kernel_name,
        progress_bar=True,
    )
    executed_nb = nbformat.read(notebook_path, as_version=4)
    html_exporter = HTMLExporter(exclude_input=True)
    html_body, _ = html_exporter.from_notebook_node(executed_nb)
    html_path = notebook_path.replace('.ipynb', '.html')
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_body)
    click.echo(f"Notebook: {notebook_path}")
    click.echo(f"Notebook HTML: {html_path}")



@beast_pype.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('report_template', required=1,
                type=click.Choice(default_report_names))
@click.option('--output', '-o', required=True, type=str,
              help='Path to save the output report notebook (.ipynb).')
@click.option('--merged_log', '-l', multiple=True, required=True,
              type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
              help='Path to a merged log file (.csv, .tsv or .log). Can be specified multiple times. '
                   'When more than one is provided, comparative (xml-set) visualisations are used.')
@click.option('--beast_xml', '-x', multiple=True, required=True,
              type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
              help='Path to a BEAST 2 XML file. Can be specified multiple times.')
@click.option('--kernel_name', '-k', default='beast_pype', type=str,
              help='Name of the Jupyter kernel to use when executing the notebook.\n' +
                   'If not given "beast_pype" is used.')
@click.option('--xml_set_label', default=None, type=str,
              help='Label for xml set grouping variable. Passed as a parameter '
                   'to the executed notebook.')
def parameters_report(report_template,
                      output,
                      merged_log,
                      beast_xml,
                      kernel_name,
                      xml_set_label):
    """
    REPORT_TEMPLATE: Report template to use for generating the parameters report.

    Generates a parameters report notebook from a report template,
    executes it, and exports it to HTML.
    """
    merged_log_paths = list(merged_log) if len(merged_log) > 1 else merged_log[0]
    beast_xml_paths = list(beast_xml) if len(beast_xml) > 1 else beast_xml[0]
    results = gen_parameters_report(
        report_template=report_template,
        output_path=output,
        merged_log_paths=merged_log_paths,
        beast_xml_paths=beast_xml_paths,
        kernel_name=kernel_name,
        xml_set_label=xml_set_label,
    )
    click.echo(f"Notebook: {results['notebook']}")
    click.echo(f"Notebook HTML: {results['notebook_html']}")



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
@click.option('--topology', '-t', default='CCD0', type=str,
              help='Topology to use for merged BEAST tree summarization and plotting, e.g. "MCC" or "CCD0".\n' +
                   'See BEAST2\'s TreeAnnotator documentation or https://www.beast2.org/2024/06/24/what-is-new-in-v2.7.7.html for more details on tree summarization methods.\n' +
                   'If not given "CCD0" is used.'
              )
@click.option('--beast-xml-path', '-r', default=None, type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
              help='Path to BEAST XML file to use for diagnostics.'
              )
def diagnose_results(report_template,
                     beast_outputs,
                     parameters,
                     parameters_file,
                     kernel_name,
                     topology,
                     beast_xml_path):
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

    for param, value in parameters_final.items():
        if param not in diag_valid_params:
            raise ValueError(f"Parameter {param} is not a valid parameter for use in the diagnostic workflow.")
    if report_template is None:
        raise ValueError('A value for --report_template is required for the diagnostic workflow.')
    if beast_outputs is None:
        raise ValueError('A value for --beast_outputs is required for the diagnostic workflow.')

    gen_beast_diagnostic_nb(
        beast_outputs=beast_outputs,
        parameters_report_template=report_template,
        beast_xml_path=beast_xml_path,
        kernel_name=kernel_name,
        topology=topology,
        **parameters_final)
