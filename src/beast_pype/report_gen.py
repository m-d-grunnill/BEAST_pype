from PIL.Image import init
import pandas as pd
from beast_pype.outputs import read_log_file
import nbformat as nbf
import re
import os
from beast_pype.path_utils import path_to_workflow_modules, path_to_report_templates
from copy import deepcopy


def add_unreported_outputs(notebook_template_path,
                           merged_logs_path,
                           output_path,
                           xml_set_comparisons=False,
                           as_version=4):
    """
    Add sections for unreported outputs to a notebook.

    Parameters
    ----------
    notebook_template_path: str
        Path to the notebook template. Metadata of the notebook will be searched for
        the entry "BEAST outputs reported". Any parameter not listed but occurring
        in the file at merged_log_path will then be added to the notebook.
    merged_logs_path: str
        Path to directory containing merged log files you wish to report on.
    output_path: str
        Path to save new report notebook file
    xml_set_comparisons: bool, default False
        Are the sections to be added for xml_set comparisons.
    as_version: int
        Ipnyb notebook version.

    Returns
    -------
    None
    """
    report = nbf.read(notebook_template_path, as_version=as_version)
    columns_already_reported = ['xml set', 'Sample']
    columns_already_reported += report.metadata["BEAST outputs reported"]['parameters']
    merged_log_paths = [file for file in os.listdir(merged_logs_path) if file.endswith('merged.log')]
    merged_log_columns = []
    for merged_log_path in merged_log_paths:
        if merged_log_path.endswith('.log'):
            merged_log_df = read_log_file(f'{merged_logs_path}/{merged_log_path}')
        elif merged_log_path.endswith('.csv'):
            merged_log_df = pd.read_csv(f'{merged_logs_path}/{merged_log_path}')
        elif merged_log_path.endswith('.tsv'):
            merged_log_df = pd.read_csv(f'{merged_logs_path}/{merged_log_path}', sep='\t')
        else:
            raise ValueError(f'Only csv, tsv or log files are supported for merged_log_path ({merged_logs_path}/{merged_log_path}).')
        columns_already_reported += [
            col for col in merged_log_df.columns
            if col.startswith(('Unnamed',
                *report.metadata["BEAST outputs reported"]["parameters starting with"]))
        ]
        merged_log_columns += merged_log_df.columns.to_list()
    columns_to_report = set(merged_log_columns) - set(columns_already_reported)
    columns_to_report = sorted(columns_to_report)
    for column in columns_to_report:
        report['cells'].append(nbf.v4.new_markdown_cell(f"## {column}"))
        variable_pre_fix = re.sub('[^a-zA-Z0-9]', '_', column)
        if xml_set_comparisons:
            report['cells'].append(
                nbf.v4.new_code_cell(
                    f"{variable_pre_fix}_ax = plot_comparative_box_violin(df_melted_for_seaborn, '{column}', xml_set_label=xml_set_label)\n" +
                    f"{variable_pre_fix}_hdi_df = hdi_pivot(df, '{column}', xml_set_label=xml_set_label)\n" +
                    f"display({variable_pre_fix}_hdi_df)")
            )
        else:
            report['cells'].append(
                nbf.v4.new_code_cell(
                    f"{variable_pre_fix}_fig, ax, {variable_pre_fix}_hdi_est ="+
                    f"plot_hist_kde(trace_df=trace_df, parameter='{column}', hdi_prob=0.95)\n"+
                    f"display({variable_pre_fix}_hdi_est)")
            )
    with open(output_path, 'w') as f:
        nbf.write(report, f)


def _update_notebook_metadata(notebook_path, metadata_update, output_path = None, as_version=4):
    notebook = nbf.read(notebook_path, as_version=as_version)
    notebook.metadata.update(metadata_update)
    if output_path is None:
        output_path = notebook_path
    with open(output_path , 'w') as f:
        nbf.write(notebook, f)

def gen_mcc_notebook(merged_trees_path, output_path, kernel_name='beast_pype', as_version=4):
    workflow_modules_path = path_to_workflow_modules()
    notebook = nbf.read(f"{workflow_modules_path}/Phase-5ii-Gen-MCC-Trees.ipynb", as_version=as_version)
    notebook['cells'].append(nbf.v4.new_code_cell(
        f'source activate {kernel_name}'
    ))
    merged_tree_paths = [file for file in os.listdir(merged_trees_path) if file.endswith('merged.trees')]
    for merged_tree_path in merged_tree_paths:
        prefix = merged_tree_path.replace('merged.trees', '')
        if prefix != '':
            notebook['cells'].append(nbf.v4.new_markdown_cell(
                f"## Producing {prefix.replace('_','')} MCC tree."
            ))
        notebook['cells'].append(nbf.v4.new_code_cell(
            f"treeannotator -burnin 0 {merged_trees_path}/{prefix}merged.trees {merged_trees_path}/{prefix}mcc_tree.nexus"
        ))
    with open(output_path , 'w') as f:
        nbf.write(notebook, f)

def gen_tree_report(output_report_path, xml_set_comparisons=False, as_version=4):
    """Generate a notebook report for tree outputs of BEAST_pype.

    Parameters
    ----------
    output_report_path : str
        Path to save the report notebook to.
    xml_set_comparisons : bool, optional
        Wheather xml_set comparisons were made, by default False.
    as_version : int, optional
        Jupyter notebook version to use, by default 4
    """
    if xml_set_comparisons:
        mcc_files = [f for f in os.listdir("outputs_and_reports") if f.endswith("mcc_tree.nexus")]
        xml_set_dict = {f"## {file.replace('_mcc_tree.nexus', '')}\n" :
                        {'directory': file.replace('_mcc_tree.nexus', ''),
                        'mcc_file':  file}
                        for file in mcc_files}
        heading_hashes = '###'
    else:
        xml_set_dict = {'' : {'directory': os.getcwd(), 'mcc_file': "mcc_tree.nexus"}}
        heading_hashes = '##'

    
    start_nb = nbf.read(f"{path_to_report_templates()}/Trees_Report_Start.ipynb", as_version=as_version)
    tree_report_nb = deepcopy(start_nb)
    mcc_nb = nbf.read(f"{path_to_report_templates()}/MCC_Report_Start.ipynb", as_version=as_version)
    initial_tree_nb = nbf.read(f"{path_to_report_templates()}/Initial_Trees_Report.ipynb", as_version=as_version)
    downsampled_nb = nbf.read(f"{path_to_report_templates()}/Downsampled_Trees_Report.ipynb", as_version=as_version)
    for xml_set_label, xml_set_sub_dict in xml_set_dict.items():
        if os.file.exists(f"{xml_set_sub_dict['directory']}/downsampled_metadata.csv"):
            mcc_metadata = f"{xml_set_sub_dict['directory']}/downsampled_metadata.csv"
        else:
            mcc_metadata = f"{xml_set_sub_dict['directory']}/metadata.csv"
        tree_report_nb['cells'] += [
            nbf.v4.new_markdown_cell(f"{xml_set_label}{heading_hashes} MCC Summary Tree Plot of BEAST Runs"),
            nbf.v4.new_code_cell(
                f"mcc_tree_path <- 'outputs_and_reports/{xml_set_sub_dict['mcc_file']}'\n" +
                f"mcc_metadata_path <- '{mcc_metadata}'"            
            )] + mcc_nb['cells']

        if os.path.exists(f"{xml_set_sub_dict['directory']}/downsampled_tree"):
            tree_report_nb['cells'] += [
                nbf.v4.new_markdown_cell(f"{heading_hashes} Downsampled Tree Plot"),
                nbf.v4.new_code_cell(
                    f"downsampled_tree_node_ci_path <- '{xml_set_sub_dict['directory']}/downsampled_tree/treetime_node_confidence.csv'\n" +
                    f"downsampled_metadata_path <- '{xml_set_sub_dict['directory']}/metadata.csv'\n" +
                    f"downsampled_tree_path <- '{xml_set_sub_dict['directory']}/downsampled_tree/treetime.nwk'"
                )]
            tree_report_nb['cells'] += downsampled_nb['cells']
            tree_report_nb['cells'] +=  [
                nbf.v4.new_markdown_cell(
                    f"{heading_hashes}# Root-to-tip Plot & Stats.\n\n" +
                    f"![]({xml_set_sub_dict['directory']}/downsampled_initial_tree/treetime_root_to_tip.png)"
                ),
                nbf.v4.new_code_cell(
                    f"cat(readLines('{xml_set_sub_dict['directory']}/downsampled_initial_tree/treetime_clock_model_stats.yml'), sep = '\\n')"
                )]

        if os.path.exists(f"{xml_set_sub_dict['directory']}/initial_tree"):
            tree_report_nb['cells'] += [
                nbf.v4.new_markdown_cell(f"{heading_hashes} Initial Tree Plot"),
                nbf.v4.new_code_cell(
                    f"initial_tree_node_ci_path <- '{xml_set_sub_dict['directory']}/initial_tree/treetime_node_confidence.csv'\n" +
                    f"initial_metadata_path <- '{xml_set_sub_dict['directory']}/metadata.csv'\n" +
                    f"initial_tree_path <- '{xml_set_sub_dict['directory']}/initial_tree/treetime.nwk'"
                )]
            tree_report_nb['cells'] += initial_tree_nb['cells']
            tree_report_nb['cells'] +=  [
                nbf.v4.new_markdown_cell(
                    f"{heading_hashes}# Root-to-tip Plot & Stats.\n\n" +
                    f"![]({xml_set_sub_dict['directory']}/initial_tree/treetime_root_to_tip.png)"
                ),
                nbf.v4.new_code_cell(
                    f"cat(readLines('{xml_set_sub_dict['directory']}/initial_tree/treetime_clock_model_stats.yml'), sep = '\\n')"
                )]

    with open(output_report_path , 'w') as f:
        nbf.write(tree_report_nb, f)