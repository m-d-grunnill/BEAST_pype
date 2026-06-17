"""Generate a static diagnostic notebook and merged outputs from BEAST 2 runs."""

import os
import glob
import nbformat as nbf
import arviz as az

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
):
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

    Returns
    -------
    dict
        Paths to output files: 'notebook', 'merged_log', 'merged_trees'.
    """
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")
    if not (0 <= burnin < 100):
        raise ValueError("burnin must be between 0 and 100 (exclusive).")

    if output_prefix is None:
        output_prefix = os.path.join(directory, "static_diag_")

    # --- 1. Load log files as posterior ---
    log_files = sorted(glob.glob(os.path.join(directory, "*.log")))
    if not log_files:
        raise FileNotFoundError(f"No .log files found in: {directory}")

    posterior = read_log_files_as_posterior(log_files)

    # --- 2. Apply burn-in ---
    burned_posterior = burn_posterior(posterior, in_percentage=burnin)

    # --- 3. Get parameters and draws ---
    parameters = [
        var for var in burned_posterior.posterior.data_vars if var != "draw"
    ]
    draws = posterior.posterior["draw"].values
    n_draws = len(draws)

    # --- 4. Generate static diagnostic notebook ---
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
            f"log_files = {log_files!r}\n\n"
            f"posterior = read_log_files_as_posterior(log_files)\n"
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
                f"fig, axes = plot_traces(burned_posterior, {section_params!r})\n"
                "fig"
            )
        )

        nb["cells"].append(
            nbf.v4.new_code_cell(
                f"az.summary(burned_posterior, var_names={section_params!r}, kind='diagnostics')"
            )
        )

    # --- 5. Save notebook ---
    notebook_path = f"{output_prefix}diagnostic.ipynb"
    with open(notebook_path, "w", encoding="utf-8") as f:
        nbf.write(nb, f)

    # --- 6. Merge logs to CSV ---
    merged_log_path = f"{output_prefix}merged_logs.csv"
    merge_logs_to_csv(burned_posterior, output_file=merged_log_path)

    # --- 7. Merge .trees files at burn-in ---
    tree_files = sorted(glob.glob(os.path.join(directory, "*.trees")))
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

    return {
        "notebook": notebook_path,
        "merged_log": merged_log_path,
        "merged_trees": merged_trees_path,
    }