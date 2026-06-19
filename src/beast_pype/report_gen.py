from PIL.Image import init
import pandas as pd
import shutil
from beast_pype.outputs import read_log_file
from beast_pype.nb_utils import execute_notebook
from beast2xml import BEAST2XML
import nbformat as nbf
import re
import os
from beast_pype.path_utils import path_to_workflow_modules, path_to_report_templates
from copy import deepcopy
from nbconvert import HTMLExporter
from tqdm.auto import tqdm


def add_unreported_outputs(notebook_template_path,
                    output_path,
                    merged_log_paths=None,
                    directory_of_merged_logs=None,
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
    output_path: str
        Path to save new report notebook file
    merged_log_paths: str or list(str), optional
        Paths to merged log file you wish to report on. If not given, all merged log files in directory_of_merged_logs will be used.
    directory_of_merged_logs: str, optional
        Path to directory containing merged log files you wish to report on.
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
    if directory_of_merged_logs is not None:
        if merged_log_paths is not None:
            raise ValueError("Only one of merged_log_paths or directory_of_merged_logs should be provided, not both.")
        merged_log_paths = [f'{directory_of_merged_logs}/{file}'
                        for file in os.listdir(directory_of_merged_logs)
                        if file.endswith(('merged.log', 'merged_logs.csv'))]
    if merged_log_paths is None and directory_of_merged_logs is None:
        raise ValueError("At least one of merged_log_paths or directory_of_merged_logs must be provided.")
    if isinstance(merged_log_paths, str):
        merged_log_paths = [merged_log_paths]
    merged_log_columns = []
    for merged_log_path in merged_log_paths:
        if merged_log_path.endswith('.log'):
            merged_log_df = read_log_file(merged_log_path)
        elif merged_log_path.endswith('.csv'):
            merged_log_df = pd.read_csv(merged_log_path)
        elif merged_log_path.endswith('.tsv'):
            merged_log_df = pd.read_csv(merged_log_path, sep='\t')
        else:
            raise ValueError(f'Only csv, tsv or log files are supported for merged_log_path ({merged_log_path}).')
        columns_already_reported += [
            col for col in merged_log_df.columns
            if col.startswith(('Unnamed',
                *report.metadata["BEAST outputs reported"]["parameters starting with"]))
        ]
        merged_log_columns += merged_log_df.columns.to_list()
    columns_to_report = set(merged_log_columns) - set(columns_already_reported)
    columns_to_report = sorted(columns_to_report)
    if 'growthRate' in columns_to_report:
        columns_to_report.remove('growthRate')
        if xml_set_comparisons:
            report['cells'].append(nbf.v4.new_markdown_cell(f"## Growth Rate\n### Per year"))
            report['cells'].append(
                    nbf.v4.new_code_cell(
                        f"growthRate_ax = plot_comparative_box_violin(df_melted_for_seaborn, 'growthRate', xml_set_label=xml_set_label)\n" +
                        f"growthRate_hdi_df = hdi_pivot(df, 'growthRate', xml_set_label=xml_set_label)\n" +
                        f"display(growthRate_hdi_df)")
                    )
            report['cells'].append(nbf.v4.new_markdown_cell(f"### Per day\n**Note:** Growth rate per day is calculated by dividing growth rate per year by 365.25."))
            report['cells'].append(
                    nbf.v4.new_code_cell(
                        f"growthRate_per_day_ax = plot_comparative_box_violin(df_melted_for_seaborn, 'Growth Rate per day', xml_set_label=xml_set_label)\n" +
                        f"growthRate_per_day_hdi_df = hdi_pivot(df, 'Growth Rate per day', xml_set_label=xml_set_label)\n" +
                        f"display(growthRate_per_day_hdi_df)")
                    )
            report['cells'].append(nbf.v4.new_markdown_cell(f"### Doublication Time in days\n**Note:** Doubling time is calculated as ln(2) divided by growth rate per day."))
            report['cells'].append(
                    nbf.v4.new_code_cell(
                        f"doubling_time_per_day_ax = plot_comparative_box_violin(df_melted_for_seaborn, 'Doubling Time in days', xml_set_label=xml_set_label)\n" +
                        f"doubling_time_per_day_hdi_df = hdi_pivot(df, 'Doubling Time in days', xml_set_label=xml_set_label)\n" +
                        f"display(doubling_time_per_day_hdi_df)")
                    )
        else:
            report['cells'].append(nbf.v4.new_markdown_cell(f"## Growth Rate\n### Per year"))
            report['cells'].append(
                    nbf.v4.new_code_cell(
                        f"growthRate_fig, growthRate_ax, growthRate_hdi_est = plot_hist_kde(trace_df=trace_df, parameter='growthRate', hdi_prob=0.95)\n" +
                        f"display(growthRate_hdi_est)"
                    )
                )
            report['cells'].append(nbf.v4.new_markdown_cell(f"### Per day\n**Note:** Growth rate per day is calculated by dividing growth rate per year by 365.25."))
            report['cells'].append(
                    nbf.v4.new_code_cell(
                        f"growthRate_per_day_fig, growthRate_per_day_ax, growthRate_per_day_hdi_est = plot_hist_kde(trace_df=trace_df, parameter='Growth Rate per day', hdi_prob=0.95)\n" +
                        f"display(growthRate_per_day_hdi_est)"
                    )
                )
            report['cells'].append(nbf.v4.new_markdown_cell(f"### Doublication Time in days\n**Note:** Doubling time is calculated as ln(2) divided by growth rate per day."))
            report['cells'].append(
                    nbf.v4.new_code_cell(
                        f"doubling_time_per_day_fig, doubling_time_per_day_ax, doubling_time_per_day_hdi_est = plot_hist_kde(trace_df=trace_df, parameter='Doubling Time in days', hdi_prob=0.95)\n" +
                        f"display(doubling_time_per_day_hdi_est)"
                    ))

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

def gen_summary_tree_notebook(directory_of_merged_trees, output_path, topology, summary_tree_low_memory=False, kernel_name='beast_pype', as_version=4):
    workflow_modules_path = path_to_workflow_modules()
    notebook = nbf.read(f"{workflow_modules_path}/Phase-5i-Gen-Summary-Trees.ipynb", as_version=as_version)
    notebook['cells'].append(nbf.v4.new_code_cell(
        f'source activate {kernel_name}'
    ))
    merged_tree_paths = [file 
        for file in os.listdir(directory_of_merged_trees) 
        if file.endswith(('merged.trees', 'merged_trees.trees'))]
    for merged_tree_path in merged_tree_paths:
        prefix = merged_tree_path.replace('merged.trees', '')
        prefix = prefix.replace('merged_trees.trees', '')
        if prefix != '':
            notebook['cells'].append(nbf.v4.new_markdown_cell(
                f"## Producing {prefix.replace('_','')} Summary Tree."
            ))
        notebook['cells'].append(nbf.v4.new_code_cell(
            f"treeannotator -burnin 0 -topology {topology} -lowMem {str(summary_tree_low_memory).lower()} " +
            f"{directory_of_merged_trees}/{merged_tree_path} {directory_of_merged_trees}/{prefix}{topology}_summary_tree.nexus"
        ))
    with open(output_path , 'w') as f:
        nbf.write(notebook, f)

def gen_tree_report(output_report_path, inputs_directory='outputs_and_reports', topology=None, collection_date_field='collection_date', plot_width=12, plot_height=6, plot_res=120, xml_set_comparisons=False, as_version=4):
    """Generate a notebook report for tree outputs of BEAST_pype.

    Parameters
    ----------
    output_report_path : str
        Path to save the report notebook to.
    inputs_directory : str, default 'outputs_and_reports'
        Path to directory containing tree outputs to report on. 
    topology : str, optional
        Topology to use for merged BEAST tree summarization and plotting, e.g. "MCC" or "CCD0". 
        See BEAST2's TreeAnnotator documentation or https://www.beast2.org/2024/06/24/what-is-new-in-v2.7.7.html for more details on tree summarization methods.
    collection_date_field : str, default 'collection_date'
        Name of field in metadata_db containing collection dates of sequences. Should be formatted YYYY-MM-DD
    plot_width : float or int, optional
        Width of tree plots in inches, by default 12. See https://search.r-project.org/CRAN/refmans/repr/html/repr-options.html.
    plot_height : float or int, optional
        Height of tree plots in inches, by default 6. See https://search.r-project.org/CRAN/refmans/repr/html/repr-options.html.
    plot_res : int, optional
        PPI for rasterization (resolution), by default 120. See https://search.r-project.org/CRAN/refmans/repr/html/repr-options.html.
    xml_set_comparisons : bool, optional
        Whether xml_set comparisons were made, by default False.
    as_version : int, optional
        Jupyter notebook version to use, by default 4
    """
    if topology is None:
        topology = ''
    if xml_set_comparisons:
        summary_tree_files = [f for f in os.listdir(inputs_directory) if f.endswith("summary_tree.nexus")]
        xml_set_dict = {f"## {file.replace(f'_{topology}_summary_tree.nexus', '')}\n" :
                        {'directory': file.replace(f'_{topology}_summary_tree.nexus', ''),
                        'summary_file':  file,
                        'summary_data_path': f"{inputs_directory}/{file.replace('_summary_tree.nexus', '')}_summary_tree_data.csv"
                        }
                        for file in summary_tree_files}
        heading_hashes = '###'
    else:
        xml_set_dict = {'' : {'directory': os.getcwd(),
         'summary_file': f"{topology}_summary_tree.nexus",
         'summary_data_path': f"{inputs_directory}/{topology}_summary_tree_data.csv"}}
        heading_hashes = '##'

    tree_report_components_path = f"{path_to_report_templates()}/tree_report_components"
    start_nb = nbf.read(f"{tree_report_components_path}/Trees_Report_Start.ipynb", as_version=as_version)
    tree_report_nb = deepcopy(start_nb)
    summary_nb = nbf.read(f"{tree_report_components_path}/Summary_Tree_Report.ipynb", as_version=as_version)
    initial_tree_nb = nbf.read(f"{tree_report_components_path}/Initial_Trees_Report.ipynb", as_version=as_version)
    downsampled_nb = nbf.read(f"{tree_report_components_path}/Downsampled_Initial_Trees_Report.ipynb", as_version=as_version)

    tree_report_nb['cells'].append(
        nbf.v4.new_code_cell(
            "# Set plot size for tree plots, see https://search.r-project.org/CRAN/refmans/repr/html/repr-options.html.\n" +
            f"options(repr.plot.width = {str(plot_width)}, repr.plot.height = {str(plot_height)}, repr.plot.res={str(plot_res)})\n" +
            f"collection_date_field <- '{collection_date_field}' # Set collection date field used in matadata."))
    for xml_set_label, xml_set_sub_dict in xml_set_dict.items():
        beast_xml = BEAST2XML(f"{xml_set_sub_dict['directory']}/beast.xml")
        youngest_tip_year_decimal = beast_xml.extract_youngest_year_decimal()
        tree_report_nb['cells'] += [
            nbf.v4.new_markdown_cell(f"{xml_set_label}{heading_hashes} {topology} Summary Tree Plot of BEAST Runs"),
            nbf.v4.new_code_cell(
                f"summary_tree_path <- '{inputs_directory}/{xml_set_sub_dict['summary_file']}'\n" +
                f"beast_youngest_tip_decimal <- {str(youngest_tip_year_decimal)}" +
                    " # Extracted using python via BEAST2XML.extract_youngest_year_decimal from the package beast2xml.\n" +
                f"summary_tree_data_path <- '{xml_set_sub_dict['summary_data_path']}'"          
            )] + summary_nb['cells']

        if os.path.exists(f"{xml_set_sub_dict['directory']}/downsampled_initial_tree"):
            tree_report_nb['cells'] += [
                nbf.v4.new_markdown_cell(f"{heading_hashes} Downsampled Initial Tree Plot"),
                nbf.v4.new_code_cell(
                    f"downsampled_tree_node_ci_path <- '{xml_set_sub_dict['directory']}/downsampled_initial_tree/treetime_node_confidence.csv'\n" +
                    f"downsampled_metadata_path <- '{xml_set_sub_dict['directory']}/downsampled_metadata.csv'\n" +
                    f"downsampled_tree_path <- '{xml_set_sub_dict['directory']}/downsampled_initial_tree/treetime.nwk'\n" +
                    f"downsampled_tree_rtp <- '{xml_set_sub_dict['directory']}/downsampled_initial_tree/treetime_root_to_tip.png'\n" +
                    f"downsampled_tree_stats <- '{xml_set_sub_dict['directory']}/downsampled_initial_tree/treetime_clock_model_stats.yml'"

                )]
            tree_report_nb['cells'] += downsampled_nb['cells']
            tree_report_nb['cells'] +=  [
                nbf.v4.new_markdown_cell(
                    f"{heading_hashes}# Root-to-tip Plot & Stats.\n\n" 
                ),
                nbf.v4.new_code_cell(
                    "img <- png::readPNG(downsampled_tree_rtp)\n" +
                    "grid::grid.raster(img)"
                ),
                nbf.v4.new_code_cell(
                    f"cat(readLines(downsampled_tree_stats), sep = '\\n')"
                )]

        if os.path.exists(f"{xml_set_sub_dict['directory']}/initial_tree"):
            tree_report_nb['cells'] += [
                nbf.v4.new_markdown_cell(f"{heading_hashes} Initial Tree Plot"),
                nbf.v4.new_code_cell(
                    f"initial_tree_node_ci_path <- '{xml_set_sub_dict['directory']}/initial_tree/treetime_node_confidence.csv'\n" +
                    f"initial_metadata_path <- '{xml_set_sub_dict['directory']}/metadata.csv'\n" +
                    f"initial_tree_path <- '{xml_set_sub_dict['directory']}/initial_tree/treetime.nwk'\n" +
                    f"initial_tree_rtp <- '{xml_set_sub_dict['directory']}/initial_tree/treetime_root_to_tip.png'\n" +
                    f"initial_tree_stats <- '{xml_set_sub_dict['directory']}/initial_tree/treetime_clock_model_stats.yml'"
                )]
            tree_report_nb['cells'] += initial_tree_nb['cells']
            tree_report_nb['cells'] +=  [
                nbf.v4.new_markdown_cell(
                    f"{heading_hashes}# Root-to-tip Plot & Stats.\n\n"
                ),
                nbf.v4.new_code_cell(
                    "img <- png::readPNG(initial_tree_rtp)\n" +
                    "grid::grid.raster(img)"
                ),
                nbf.v4.new_code_cell(
                    f"cat(readLines(initial_tree_stats), sep = '\\n')"
                )
                ]

    with open(output_report_path , 'w') as f:
        nbf.write(tree_report_nb, f)

def gen_summary_tree_report(
    output_report_path,
    summary_trees,
    beast_xml_paths,
    collection_date_field='collection_date',
    plot_width=12,
    plot_height=6,
    plot_res=120,
    as_version=4,
    kernel_name='beast_pype_R'):
    """Generate a Jupyter notebook report summarizing BEAST 2 summary trees.

    Parameters
    ----------
    output_report_path : str
        Path to save the report notebook to.
    summary_trees : str or list of str
        Path(s) to summary tree files (nexus format). A single string is
        converted to a list.
    beast_xml_paths : str or list of str
        Path(s) to corresponding BEAST 2 XML files (same order as
        summary_trees). A single string is converted to a list.
    collection_date_field : str, default 'collection_date'
        Name of field in metadata containing collection dates (YYYY-MM-DD).
    plot_width : float or int, default 12
        Width of tree plots in inches.
    plot_height : float or int, default 6
        Height of tree plots in inches.
    plot_res : int, default 120
        PPI for rasterization.
    as_version : int, default 4
        Jupyter notebook version.
    kernel_name : str, default 'beast_pype_R'
        Name of the Jupyter kernel to use when executing the notebook.

    Returns
    -------
    dict
        Paths to output files: 'notebook' and 'notebook_html'.
    """
    if isinstance(summary_trees, str):
        summary_trees = [summary_trees]
    if isinstance(beast_xml_paths, str):
        beast_xml_paths = [beast_xml_paths]

    # Convert to absolute paths so the notebook works regardless of execution directory
    summary_trees = [os.path.abspath(p) for p in summary_trees]
    beast_xml_paths = [os.path.abspath(p) for p in beast_xml_paths]

    if len(summary_trees) != len(beast_xml_paths):
        raise ValueError(
            "summary_trees and beast_xml_paths must have the same length."
        )

    tree_report_components_path = f"{path_to_report_templates()}/tree_report_components"
    start_nb = nbf.read(
        f"{tree_report_components_path}/Trees_Report_Start.ipynb", as_version=as_version
    )
    summary_nb = nbf.read(
        f"{tree_report_components_path}/Summary_Tree_Report.ipynb", as_version=as_version
    )

    tree_report_nb = deepcopy(start_nb)

    tree_report_nb['cells'].append(
        nbf.v4.new_code_cell(
            f"options(repr.plot.width = {plot_width}, repr.plot.height = {plot_height}, repr.plot.res = {plot_res})\n"
            f"collection_date_field <- '{collection_date_field}'"
        )
    )

    for summary_tree_path, beast_xml_path in zip(summary_trees, beast_xml_paths):
        beast_xml = BEAST2XML(beast_xml_path)
        youngest_tip_year_decimal = beast_xml.extract_youngest_year_decimal()

        label = os.path.basename(summary_tree_path).replace("_summary_tree.nexus", "").replace(".nexus", "")
        summary_data_path = summary_tree_path.replace(".nexus", "_data.csv")

        tree_report_nb['cells'].append(
            nbf.v4.new_markdown_cell(f"## {label} Summary Tree")
        )
        tree_report_nb['cells'].append(
            nbf.v4.new_code_cell(
                f"summary_tree_path <- '{summary_tree_path}'\n"
                f"beast_youngest_tip_decimal <- {youngest_tip_year_decimal}\n"
                f"summary_tree_data_path <- '{summary_data_path}'"
            )
        )
        tree_report_nb['cells'] += deepcopy(summary_nb['cells'])

    # --- Progress bar for save/execute/export steps ---
    steps = ['Saving notebook', 'Executing notebook', 'Exporting to HTML']
    pbar = tqdm(total=len(steps), desc='Summary tree report', unit='step', leave=False)

    # Save notebook
    pbar.set_postfix_str(steps[0])
    with open(output_report_path, 'w') as f:
        nbf.write(tree_report_nb, f)
    pbar.update(1)

    # Execute notebook
    pbar.set_postfix_str(steps[1])
    execute_notebook(
        input_path=output_report_path,
        output_path=output_report_path,
        kernel_name=kernel_name,
        progress_bar=True,
    )
    pbar.update(1)

    # Export to HTML (excluding code cells)
    pbar.set_postfix_str(steps[2])
    tree_report_nb = nbf.read(output_report_path, as_version=4)
    html_exporter = HTMLExporter(exclude_input=True)
    html_body, _ = html_exporter.from_notebook_node(tree_report_nb)
    html_path = output_report_path.replace(".ipynb", ".html")
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_body)
    pbar.update(1)
    pbar.close()

    return {
        "notebook": output_report_path,
        "notebook_html": html_path,
    }


def gen_parameters_report(
    report_template,
    output_path,
    merged_log_paths,
    beast_xml_paths,
    kernel_name='beast_pype',
    xml_set_label=None,
):
    """Generate a parameters report from a report template, then execute and export to HTML.

    This function:
    1. Validates and copies the report template to ``output_path``.
    2. Calls :func:`add_unreported_outputs` to append sections for any
       parameters found in the merged logs but not already covered by the
       template.
    3. Parameterises and executes the notebook with papermill via
       :func:`~beast_pype.nb_utils.execute_notebook`.
    4. Converts the executed notebook to HTML.

    Parameters
    ----------
    report_template : str
        Name of the report template (without the ``.ipynb`` extension).
        Must match a file in the report templates directory returned by
        :func:`~beast_pype.path_utils.path_to_report_templates`.
    output_path : str
        Path (including filename) where the output notebook will be saved.
    merged_log_paths : str or list of str
        Path(s) to merged log files (``.csv``, ``.tsv`` or ``.log``).
        A single string is converted to a one-element list.
    beast_xml_paths : str or list of str
        Path(s) to the BEAST 2 XML file(s) used to generate the outputs.
        A single string is converted to a one-element list.
    kernel_name : str, default 'beast_pype'
        Name of the Jupyter kernel to use when executing the notebook.
    xml_set_label : str, optional
        Label identifying the xml set grouping variable. Passed as a parameter
        to the executed notebook. Comparative (xml-set) visualisations are
        triggered automatically when more than one merged log path is provided.

    Returns
    -------
    dict
        Paths to output files with keys ``'notebook'`` and ``'notebook_html'``.
    """
    # --- Validate inputs ---
    if isinstance(merged_log_paths, str):
        merged_log_paths = [merged_log_paths]
    if isinstance(beast_xml_paths, str):
        beast_xml_paths = [beast_xml_paths]

    reports_path = path_to_report_templates()
    available = [f.stem for f in reports_path.glob('*.ipynb')]
    if report_template not in available:
        raise ValueError(
            f"report_template '{report_template}' not found. "
            f"Available templates: {', '.join(available)}"
        )

    template_path = f"{reports_path}/{report_template}.ipynb"

    # --- Progress bar for steps ---
    steps = [
        'Copying template',
        'Adding unreported outputs',
        'Parameterising and executing notebook',
        'Exporting to HTML',
    ]
    pbar = tqdm(total=len(steps), desc='Parameters report', unit='step', leave=False)

    # --- 1. Copy template to output_path ---
    pbar.set_postfix_str(steps[0])
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    shutil.copy2(template_path, output_path)
    pbar.update(1)

    # --- 2. Add unreported outputs ---
    pbar.set_postfix_str(steps[1])
    xml_set_comparisons = len(merged_log_paths) > 1
    add_unreported_outputs(
        notebook_template_path=output_path,
        output_path=output_path,
        merged_log_paths=merged_log_paths,
        xml_set_comparisons=xml_set_comparisons,
    )
    pbar.update(1)

    # --- 3. Parameterise and execute ---
    pbar.set_postfix_str(steps[2])
    parameters = {
        'merged_log_paths': merged_log_paths if len(merged_log_paths) > 1 else merged_log_paths[0],
        'beast_xml_paths': beast_xml_paths if len(beast_xml_paths) > 1 else beast_xml_paths[0],
    }
    if xml_set_label is not None:
        parameters['xml_set_label'] = xml_set_label

    execute_notebook(
        input_path=output_path,
        output_path=output_path,
        parameters=parameters,
        kernel_name=kernel_name,
        progress_bar=True,
    )
    pbar.update(1)

    # --- 4. Convert to HTML (excluding code cells) ---
    pbar.set_postfix_str(steps[3])
    executed_nb = nbf.read(output_path, as_version=4)
    html_exporter = HTMLExporter(exclude_input=True)
    html_body, _ = html_exporter.from_notebook_node(executed_nb)
    html_path = output_path.replace('.ipynb', '.html')
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_body)
    pbar.update(1)
    pbar.close()

    return {
        'notebook': output_path,
        'notebook_html': html_path,
    }