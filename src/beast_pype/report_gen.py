import pandas as pd
from beast_pype.outputs import read_log_file
import nbformat as nbf

def add_unreported_outputs(notebook_template_path,
                           merged_log_path,
                           output_path,
                           as_version=4):
    """
    Add sections for unreported outputs to a notebook.

    Parameters
    ----------
    notebook_template_path: str
        Path to the notebook template. Metadata of the notebook will be searched for
        the entry "BEAST outputs reported". Any parameter not listed but occurring
        in the file at merged_log_path will then be added to the notebook.
    merged_log_path: str
        Path to merged log file you wish to report on.
    output_path: str
        Path to save new report notebook file
    as_version: int
        Ipnyb notebook version.

    Returns
    -------
    None
    """
    report = nbf.read(notebook_template_path, as_version=as_version)
    if merged_log_path.endswith('.log'):
        merged_log_df = read_log_file(merged_log_path)
    elif merged_log_path.endswith('.csv'):
        merged_log_df = pd.read_csv(merged_log_path)
    elif merged_log_path.endswith('.tsv'):
        merged_log_df = pd.read_csv(merged_log_path, sep='\t')
    else:
        raise ValueError(f'Only csv, tsv or log files are supported for merged_log_path ({merged_log_path}).')
    columns_already_reported = ['xml set', 'Sample']
    columns_already_reported += report.metadata["BEAST outputs reported"]['parameters']
    columns_already_reported += [
        col for col in merged_log_df.columns
        if col.startswith(('Unnamed',
            *report.metadata["BEAST outputs reported"]["parameters starting with"]))
    ]
    columns_to_report = set(merged_log_df.columns) - set(columns_already_reported)
    columns_to_report = sorted(columns_to_report)
    for column in columns_to_report:
        report['cells'].append(nbf.v4.new_markdown_cell(f"## {column}"))
        report['cells'].append(
            nbf.v4.new_code_cell(
                f"{column.replace('.','_')}_fig, ax, {column.replace('.','_')}_hdi_est ="+
                f"plot_hist_kde(trace_df=trace_df, parameter='{column}', hdi_prob=0.95)\n"+
                f"display({column.replace('.','_')}_hdi_est)")
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