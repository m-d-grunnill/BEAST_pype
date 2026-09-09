"""Generate an MCMC report notebook from BEAST 2 runs."""

import os
import re
import glob
import nbformat as nbf
from nbconvert import HTMLExporter
from tqdm.auto import tqdm
from beast_pype.nb_utils import execute_notebook, make_kernelspec

from beast_pype.diagnostics.mcmc import (
    read_log_files_as_posterior,
    burn_posterior,
    plot_traces,
)


def gen_mcmc_report_nb(
    directory,
    burnin=10,
    front_percentage=100,
    output_prefix=None,
    parameters_per_section=1,
    kernel_name='beast_pype',
    chain_names=None):
    """Generate an MCMC report notebook from a BEAST 2 run directory.

    Given a directory containing BEAST 2 `.log` files, this function:
    1. Loads the selected `.log` files as an arviz posterior.
    2. Removes the burn-in percentage from the posterior.
    3. Creates a static Jupyter notebook with trace plots and diagnostic summaries
       for each parameter.
    4. Executes the notebook and exports it to HTML.

    Parameters
    ----------
    directory : str
        Path to directory containing BEAST 2 `.log` output files.
    burnin : int or float, default 10
        Burn-in percentage (0-100) to remove from the start of each chain.
    front_percentage : int or float, default 100
        Percentage (0-100) of the posterior to keep at the front. Passed through
        to ``burn_posterior``. Must be greater than or equal to ``burnin``.
    output_prefix : str, optional
        Prefix (including path) for output files. Defaults to `directory/static_diag_`.
    parameters_per_section : int, default 1
        Number of parameters to display per notebook section.
    kernel_name : str, default 'beast_pype'
        Name of the Jupyter kernel to use when executing the notebook.
    chain_names : list of str, optional
        Names of the MCMC chains to include. For each name, the associated log file
        is assumed to be `<name>.log`, located in `directory`. If None (default),
        all `.log` files found in `directory` are used.

    Returns
    -------
    dict
        Paths to output files: 'notebook', 'notebook_html'.
    """
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")
    if not (0 <= burnin < 100):
        raise ValueError("burnin must be between 0 and 100 (exclusive).")
    if not (0 < front_percentage <= 100):
        raise ValueError("front_percentage must be between 0 and 100 (inclusive).")
    if front_percentage < burnin:
        raise ValueError("front_percentage must be greater than or equal to burnin.")

    if output_prefix is None:
        output_prefix = os.path.join(directory, "static_diag_")

    steps = [
        'Loading log files',
        'Applying burn-in',
        'Generating notebook',
        'Saving notebook',
        'Executing notebook',
        'Exporting to HTML',
    ]
    pbar = tqdm(total=len(steps), desc='MCMC report', unit='step')

    # --- 1. Load log files as posterior ---
    pbar.set_postfix_str(steps[0])
    if chain_names is not None:
        log_files = []
        for name in chain_names:
            log_file = os.path.join(directory, f"{name}.log")
            if not os.path.isfile(log_file):
                pbar.close()
                raise FileNotFoundError(f"Log file not found: {log_file}")
            log_files.append(log_file)
        log_paths = {
            name: os.path.abspath(os.path.join(directory, f"{name}.log"))
            for name in chain_names
        }
    else:
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
    burned_posterior = burn_posterior(
        posterior, in_percentage=burnin, front_percentage=front_percentage
    )
    pbar.update(1)

    # --- 3. Get parameters ---
    parameters = [
        var for var in burned_posterior.posterior.data_vars if var != "draw"
    ]

    # --- 4. Generate static diagnostic notebook ---
    pbar.set_postfix_str(steps[2])
    nb = nbf.v4.new_notebook()
    nb["metadata"]["kernelspec"] = make_kernelspec(kernel_name)
    nb["cells"] = []

    nb["cells"].append(
        nbf.v4.new_markdown_cell(
            "# Static BEAST 2 MCMC Diagnostic Report\n\n"
            f"**Directory:** `{directory}`\n\n"
            f"**Burn-in:** {burnin}%\n\n"
            f"**Keep front:** {front_percentage}%\n\n"
            f"**Log files:** {len(log_files)}\n\n"
            f"**Chains:** {', '.join(map(str, log_paths.keys()))}\n\n"
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
            f"burned_posterior = burn_posterior(posterior, in_percentage={burnin}, front_percentage={front_percentage})\n"
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
    notebook_path = f"{output_prefix}MCMC_report.ipynb"
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
    html_path = f"{output_prefix}MCMC_report.html"
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_body)
    pbar.update(1)
    pbar.close()

    return {
        "notebook": notebook_path,
        "notebook_html": html_path,
    }
