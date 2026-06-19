"""Generate a static diagnostic notebook and merged outputs from BEAST 2 runs."""

import os
import re
import glob
import nbformat as nbf
import arviz as az
from nbconvert import HTMLExporter
from tqdm.auto import tqdm
from beast_pype.nb_utils import execute_notebook

from beast_pype.diagnostics.mcmc import (
    read_log_files_as_posterior,
    burn_posterior,
    merge_logs_to_csv,
    subset_and_merge_trees,
    plot_traces,
)


def gen_static_diagnostic_nb(
    directory,
    burnin=10,
    output_prefix=None,
    parameters_per_section=1,
    kernel_name='beast_pype'):
    """Generate a static diagnostic notebook and merged outputs from a BEAST 2 run directory.

    Given a directory containing BEAST 2 `.log` and `.trees` files, this function:
    1. Loads all `.log` files as an arviz posterior.
    2. Removes the burn-in percentage from the posterior.
    3. Creates a static Jupyter notebook with trace plots and diagnostic summaries
       for each parameter.
    4. Saves the burned posterior as a merged CSV.
    5. Merges all `.trees` files at the burn-in percentage.

    Parameters
    ----------
    directory : str
        Path to directory containing BEAST 2 `.log` and `.trees` output files.
    burnin : int or float, default 10
        Burn-in percentage (0-100) to remove from the start of each chain.
    output_prefix : str, optional
        Prefix (including path) for output files. Defaults to `directory/static_diag_`.
    parameters_per_section : int, default 1
        Number of parameters to display per notebook section.
    kernel_name : str, default 'beast_pype'
        Name of the Jupyter kernel to use when executing the notebook.

    Returns
    -------
    dict
        Paths to output files: 'notebook', 'notebook_html', 'merged_log', 'merged_trees'.
    """
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")
    if not (0 <= burnin < 100):
        raise ValueError("burnin must be between 0 and 100 (exclusive).")

    if output_prefix is None:
        output_prefix = os.path.join(directory, "static_diag_")

    steps = [
        'Loading log files',
        'Applying burn-in',
        'Generating notebook',
        'Saving notebook',
        'Executing notebook',
        'Exporting to HTML',
        'Merging logs to CSV',
        'Merging .trees files',
    ]
    pbar = tqdm(total=len(steps), desc='Static diagnostic', unit='step')

    # --- 1. Load log files as posterior ---
    pbar.set_postfix_str(steps[0])
    log_files = sorted(glob.glob(os.path.join(directory, "*.log")))
    if not log_files:
        pbar.close()
        raise FileNotFoundError(f"No .log files found in: {directory}")

    log_paths = {
        re.sub(r'(-BEAST)?\.log$', '', os.path.basename(f)): os.path.abspath(f)
        for f in log_files
    }
    posterior = read_log_files_as_posterior(log_paths)
    pbar.update(1)

    # --- 2. Apply burn-in ---
    pbar.set_postfix_str(steps[1])
    burned_posterior = burn_posterior(posterior, in_percentage=burnin)
    pbar.update(1)

    # --- 3. Get parameters and draws ---
    parameters = [
        var for var in burned_posterior.posterior.data_vars if var != "draw"
    ]
    draws = posterior.posterior["draw"].values
    n_draws = len(draws)

    # --- 4. Generate static diagnostic notebook ---
    pbar.set_postfix_str(steps[2])
    nb = nbf.v4.new_notebook()
    nb["cells"] = []

    nb["cells"].append(
        nbf.v4.new_markdown_cell(
            "# Static BEAST 2 MCMC Diagnostic Report\n\n"
            f"**Directory:** `{directory}`\n\n"
            f"**Burn-in:** {burnin}%\n\n"
            f"**Log files:** {len(log_files)}\n\n"
            f"**Parameters:** {len(parameters)}"
        )
    )

    nb["cells"].append(
        nbf.v4.new_code_cell(
            "import os\n"
            "import re\n"
            "import arviz as az\n"
            "from beast_pype.diagnostics.mcmc import (\n"
            "    read_log_files_as_posterior,\n"
            "    burn_posterior,\n"
            "    plot_traces,\n"
            ")\n"
            "import warnings\n"
            "warnings.filterwarnings('ignore')\n"
        )
    )

    nb["cells"].append(
        nbf.v4.new_code_cell(
            f"log_paths = {log_paths!r}\n\n"
            f"posterior = read_log_files_as_posterior(log_paths)\n"
            f"burned_posterior = burn_posterior(posterior, in_percentage={burnin})\n"
        )
    )

    # Parameter sections
    for i in range(0, len(parameters), parameters_per_section):
        section_params = parameters[i : i + parameters_per_section]
        section_label = ", ".join(section_params)

        nb["cells"].append(
            nbf.v4.new_markdown_cell(f"## Parameters: {section_label}")
        )

        nb["cells"].append(
            nbf.v4.new_code_cell(
                f"fig, axes = plot_traces(burned_posterior, {section_params!r}, labels={list(log_paths.keys())!r})"
            )
        )

        nb["cells"].append(
            nbf.v4.new_code_cell(
                f"az.summary(burned_posterior, var_names={section_params!r}, kind='diagnostics')"
            )
        )
    pbar.update(1)

    # --- 5. Save notebook ---
    pbar.set_postfix_str(steps[3])
    notebook_path = f"{output_prefix}diagnostic.ipynb"
    with open(notebook_path, "w", encoding="utf-8") as f:
        nbf.write(nb, f)
    pbar.update(1)

    # --- 5b. Execute notebook ---
    pbar.set_postfix_str(steps[4])
    execute_notebook(
        input_path=notebook_path,
        output_path=notebook_path,
        kernel_name=kernel_name,
        progress_bar=True,
    )
    pbar.update(1)

    # --- 5c. Export to HTML (exclude code cells) ---
    pbar.set_postfix_str(steps[5])
    nb = nbf.read(notebook_path, as_version=4)
    html_exporter = HTMLExporter(exclude_input=True)
    html_body, _ = html_exporter.from_notebook_node(nb)
    html_path = f"{output_prefix}diagnostic.html"
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_body)
    pbar.update(1)

    # --- 6. Merge logs to CSV ---
    pbar.set_postfix_str(steps[6])
    merged_log_path = f"{output_prefix}merged_logs.csv"
    merge_logs_to_csv(burned_posterior, output_file=merged_log_path)
    pbar.update(1)

    # --- 7. Merge .trees files at burn-in ---
    pbar.set_postfix_str(steps[7])
    tree_files = sorted(
        os.path.abspath(f) for f in glob.glob(os.path.join(directory, "*.trees"))
    )
    merged_trees_path = f"{output_prefix}merged_trees.trees"

    if tree_files:
        in_idx = round(burnin / 100 * n_draws)
        front_idx = n_draws - 1
        in_state_number = int(draws[min(in_idx, n_draws - 1)])
        front_state_number = int(draws[min(front_idx, n_draws - 1)])

        subset_and_merge_trees(
            file_list=tree_files,
            in_number=in_state_number,
            front_number=front_state_number,
            output_file=merged_trees_path,
        )
    else:
        merged_trees_path = None
    pbar.update(1)
    pbar.close()

    return {
        "notebook": notebook_path,
        "notebook_html": html_path,
        "merged_log": merged_log_path,
        "merged_trees": merged_trees_path,
    }
