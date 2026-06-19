import os
import xarray as xr
import arviz as az
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display, clear_output
import glob
from pathlib import Path
import re
from typing import List
from tqdm.auto import tqdm
import yaml


def read_log_file_as_dataframe(file_path,
                               cut_to_first=None):
    """Read a log file into a pandas DataFrame.

    Parameters
    ----------------
    file_path: str
        Path to the log file.
    cut_to_first : int, default None
        Remove Samples/links over this number.


    Returns
    -----------
    Pandas.DataFrame
    """
    trace_df = pd.read_table(file_path, sep='\t', comment='#')
    if cut_to_first is not None:
        trace_df = trace_df[trace_df['Sample'] <= cut_to_first]
    return trace_df

def read_log_files_as_xarraydata(log_files,
                                 start=0,
                                 cut_to_first=None):
    """
    Reads log files from a list and returns a xarray DataArray.

    Parameters
    ----------
    log_files: list, tuple or dict  of str
        Contains paths of log files. If list or tuple chains are indexed from
         incrementally start argument. If dictionary keys are used as the chain index.
    start : int
        Starting index of the first log file to read.
        Remove Samples/links over this number in log files.
    cut_to_first : int, default None
        Remove Samples/links over this number.

    Returns
    -------
    xarray.DataArray
    """
    dfs = []
    if isinstance(log_files, (list,tuple)):
        for i, log_file in enumerate(log_files, start=start):
            df = read_log_file_as_dataframe(log_file, cut_to_first=cut_to_first)
            df['chain'] = i
            dfs.append(df)
    elif isinstance(log_files, dict):
        for i, log_file in log_files.items():
            df = read_log_file_as_dataframe(log_file, cut_to_first=cut_to_first)
            df['chain'] = i
            dfs.append(df)
    else:
        raise TypeError('log_files must be a list, tuple or a dictionary.')

    df = pd.concat(dfs)
    df.rename(columns={'Sample': 'draw'}, inplace=True)
    df.set_index(["chain", "draw"], inplace=True)
    return xr.Dataset.from_dataframe(df)


def read_log_files_as_posterior(log_files,
                                start=0,
                                cut_to_first=None):
    """
    Reads log files from a list and returns an arviz BEASTDiag Inference  DataArray.
    Parameters
    ----------
    log_files: list of str
        List of log files.
    start : int
        Starting index of the first log file to read.
            sample_name: str, default='Sample'
    cut_to_first : int, default None
        Remove Samples/links over this number in log files.

    Returns
    -------
    arviz.InferenceData
        An arviz BEASTDiag Inference  DataArray
    """
    xdata = read_log_files_as_xarraydata(log_files,
                                         start=start,
                                         cut_to_first=cut_to_first)
    return az.InferenceData(posterior=xdata)


def burn_posterior(posterior, in_proportion=None, in_percentage=None, in_number=None,
                 front_proportion=None, front_percentage=None, front_number=None,
                 sample_name='draw'):
    """
    Perform a burn on a posterior, from start (in) to end (front).

    Parameters
    ----------
    posterior: arviz.data.inference_data.InferenceData
        DataArray with posterior. Must have 'chain' and 'draw' dimension names for in
         posterior.
    in_proportion: float, default=None
        Proportion of posterior to burn-in.
    in_percentage: float, default=None
         Percentage of posterior to burn-in.
    in_number: int
        Number of burn-in points.
    front_proportion: float, default=None
        Proportion of posterior to keep at the front.
    front_percentage: float, default=None
        Percentage of posterior to keep at the front.
    front_number: int
        Number of points at front to keep.
    sample_name: str, default='draw'
        Name of dimension of xml_set.

    Returns
    -------
    arviz.data.inference_data.InferenceData
    """
    if in_proportion is not None:
        if not isinstance(in_proportion, float):
            raise TypeError("Proportion must be a float")
        if in_number is not None or in_percentage is not None:
            raise ValueError("Only one of in_proportion, in_percentage or in_number must be provided")
        in_number = round(in_proportion * len(posterior.posterior[sample_name]))
    elif in_percentage is not None:
        if not isinstance(in_percentage, (int,float)):
            raise TypeError("in_percentage must be a float or an integer.")
        if in_number is not None or in_proportion is not None:
            raise ValueError("Only one of in_proportion, in_percentage or in_number must be provided")
        in_number = round(in_percentage/100 * len(posterior.posterior[sample_name]))
    elif in_number is not None:
        if isinstance(in_number, int):
            raise TypeError("in_umber must be an integer")
    
    if front_proportion is not None:
        if not isinstance(front_proportion, float):
            raise TypeError("Proportion must be a float")
        if front_number is not None or front_percentage is not None:
            raise ValueError("Only one of front_proportion, front_percentage or front_number must be provided")
        front_number = round(front_proportion * len(posterior.posterior[sample_name]))
    elif front_percentage is not None:
        if not isinstance(front_percentage, (int,float)):
            raise TypeError("front_percentage must be a float or an integer.")
        if front_number is not None or front_proportion is not None:
            raise ValueError("Only one of front_proportion, front_percentage or front_number must be provided")
        front_number = round(front_percentage/100 * len(posterior.posterior[sample_name]))
    elif front_number is not None:
        if isinstance(front_number, int):
            raise TypeError("front_umber must be an integer")
    
    
    selection = {sample_name: slice(in_number, front_number)}
    return posterior.isel(**selection, groups="posterior")


def select_chains_from_posterior(posterior, selection):
    """
    Select chains from a posterior.
    Parameters
    ----------
    posterior: arviz.data.inference_data.InferenceData
        DataArray with posterior. Must have 'chain' and 'draw' dimension names for in
         posterior.
    selection: list of strings or ints
        List of chains to select from posterior.

    Returns
    -------
    arviz.data.inference_data.InferenceData
    """
    selection = {'chain': selection}
    return posterior.sel(**selection, groups="posterior")


def plot_traces(posterior, parameters, labels=None, legend=True):
    """
    Plot traces (a wrapper for arviz.plot_trace with a better positioned legend).
    
    Parameters
    ----------
    posterior: arviz.data.inference_data.InferenceData
        DataArray with posterior. Must have 'chain' and 'draw' dimension names for in
         posterior.
    parameters: list of str
        List of parameters names.
    labels: list, default=None
        List of labels to use for legend.
    legend: bool, default=True
        Add legend.

    Returns
    -------
    Numpy array of matplotlib axes.
    """
    num_params = len(parameters)
    fig, axs = plt.subplots(nrows=num_params, ncols=2, figsize=(13, 2*num_params))
    if num_params == 1:
        axs = axs.reshape((1, 2))
    plt.subplots_adjust(hspace=0.4)
    traces = az.plot_trace(posterior,
                           axes=axs,
                           var_names=parameters,
                           chain_prop="color",
                           compact=True,
                           legend=legend)
    if legend:
        if labels is not None:
             axs[0][1].legend(labels=labels)
        sns.move_legend(axs[0][1], loc="upper left", bbox_to_anchor=(1, 1), ncol=1)
    for i, parameter in enumerate(parameters):
        axs[i][0].set_title(parameter, x=1.1)
        axs[i][1].set_title('')
        axs[i][1].ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    return fig, axs

STATE_RE = re.compile(r'^\s*tree\s+STATE_(\d+)\s*=\s*(.*)$', re.IGNORECASE)
STATE_LABEL_RE = re.compile(r'(\btree\s+)STATE_\d+(\s*=)', re.IGNORECASE)

def subset_and_merge_trees(
    file_list: List[str],
    in_number: int,
    front_number: int,
    output_file: str,
    relabel_start: int = 0,
    pbar=None):
    """
    Subset and merge tree states from BEAST 2 .trees files (NEXUS format).

    Parameters
    ---------------
    file_list : List[str]
        Filenames to include, in desired order.
    in_number : int
        Minimum STATE number to keep (inclusive).
    front_number : int
        Maximum STATE number to keep (inclusive).
    output_file : str
        Output NEXUS filename.
    relabel_start : int, default=0
        New STATE numbering starts from this value.

    """
    if in_number > front_number:
        raise ValueError("in_number must be <= front_number")

    out_path = Path(output_file)

    header_lines = None
    footer_lines = []
    selected_tree_lines = []
    state_step_size = None

    owns_pbar = False
    if pbar is None:
        pbar = tqdm(
            total=len(file_list),
            desc='Selecting trees from .trees files',
            unit='file',
            leave=False
        )
        owns_pbar = True

    for fpath in file_list:
        fpath = Path(fpath)
        if not fpath.is_file():
            raise FileNotFoundError(f"File not found: {fpath}")
        if fpath.suffix.lower() not in ['.trees', '.tree']:
            raise ValueError(f"File does not have .trees or .tree extension: {fpath}")

        with fpath.open("r", encoding="utf-8") as f:
            lines = f.readlines()

        first_tree_idx = None
        second_tree_idx = None
        last_tree_idx = None
        file_selected = []

        for idx, line in enumerate(lines):
            m = STATE_RE.match(line)
            if m:
                if first_tree_idx is None:
                    first_tree_idx = idx
                    first_state_num = int(m.group(1))
                elif second_tree_idx is None:
                    second_tree_idx = idx
                    second_state_num = int(m.group(1))
                    current_state_step_size = second_state_num - first_state_num
                    if state_step_size is None:
                        state_step_size = current_state_step_size
                    elif current_state_step_size != state_step_size:
                        raise ValueError(
                            f"State step size mismatch in file {fpath} at line {idx+1}.\n" +
                            f"Expected step size: {state_step_size}, found: {current_state_step_size}." +
                            "Check that all files have consistent STATE numbering."
                        )
                last_tree_idx = idx

                state_num = int(m.group(1))

                if in_number <= state_num <= front_number:
                    file_selected.append(line)

        if first_tree_idx is None:
            raise ValueError(f"No 'tree STATE_x' lines found in: {fpath}")

        header_lines_current = lines[:first_tree_idx]

        if header_lines is None:
            header_lines = header_lines_current
            footer_lines = lines[last_tree_idx + 1 :]
        else:
            if header_lines_current != header_lines:
                raise ValueError(
                    "NEXUS header mismatch before tree STATEs.\n"
                    f"Reference file: {file_list[0]}\n"
                    f"Mismatched file: {fpath}"
                )

        selected_tree_lines.extend(file_selected)
        pbar.update(1)

    if owns_pbar:
        pbar.close()

    relabeled_lines = []
    new_state = relabel_start
    for line in selected_tree_lines:
        new_line = STATE_LABEL_RE.sub(rf"\1STATE_{new_state}\2", line, count=1)
        relabeled_lines.append(new_line)
        new_state += int(state_step_size)

    if header_lines is None:
        raise RuntimeError("No input processed; header not initialized.")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        out.writelines(header_lines)
        out.writelines(relabeled_lines)
        out.writelines(footer_lines)
    
    out.close()

def merge_logs_to_csv(posterior, output_file='merged_logs.csv', like_logcombiner=True):
    """
    Merge selected log files into one csv file.

    Parameters
    ----------
    output_file : str, default='merged.log'
        Name of output file. Saved in self.directory.
    like_logcombiner : bool, default=True
        Merged output looks similar to merged output from logcombiner.

    """
    df = posterior.to_dataframe()
    if like_logcombiner:
        step = df.iloc[1, 1] - df.iloc[0, 1]
        df[df.columns[1]] = range(0, len(df) * step, step)
        if df.columns[1] != 'Sample':
            relabel_dict = {df.columns[1]: 'Sample'}
            df.rename(columns=relabel_dict, inplace=True)
        df = df.drop(columns=['chain'])
    df.to_csv(output_file, index=False) 

class BEASTDiag:
    """
    Class to perform a burn-in and chain selection on a sampled posterior from BEAST.

    Attributes
    ---------------
    directory:  str
        Directory where log and trees files are located.
    original_posterior: arviz.data.inference_data.InferenceData
        Posterior with no chains removed and no burn-in.
    parameters: list of str
        List of parameters names.
    parameters_dict: dict {str: [strs]}
        Dictionary of parameters categories.
    original_chains: list of strs
         Chains in original posterior.
    selected_posterior: arviz.data.inference_data.InferenceData
        Posterior with chains removed and burn-in applied.
    burinin_percentage: int
        Percentage burn-in applied to self.selected_posterior.
    keep_front_percentage: int
        Percentage of posterior kept at the front in self.selected_posterior.
    selected_chains: list of strs
        Chains in selected posterior.
    diagnostics_of_selection: pandas.DataFrame
        Summary of diagnostics of self.selected_posterior. Product of
        arviz.summary(self.selected_posterior, kind='diagnostics')
    output_prefix: str
        Prefix for output files.

    """
    def __init__(self,
                 directory: str,
                 output_prefix=None):
        """

        Parameters
        ----------
        directory:  str
            Directory where log and trees files are located.
        output_prefix: str, default None
            Prefix for output files.
        """
        self.directory = directory
        if not os.path.exists(directory):
            raise FileNotFoundError('directory does not exist')
        if output_prefix is None:
            self.output_prefix = ''
        else:
            self.output_prefix = output_prefix
        log_paths = {
                re.sub(r'(-BEAST)?\.log$', '', os.path.basename(path)): path
                for path in glob.glob(directory + '/*.log')
            }
        chains = list(log_paths.keys())
        self.original_posterior = read_log_files_as_posterior(log_paths)
        self.parameters = list(self.original_posterior.posterior.data_vars)
        self.parameters_dict = None
        self.original_chains = sorted(chains)
        self.selected_posterior = None
        self.burinin_percentage = 10
        self.keep_front_percentage = 100
        self.selected_chains = chains
        self.diagnostics_of_selection = None


    def set_burnin(self, in_percentage=10, front_percentage=100, sample_name='draw'):
        """
        Perform a burn-in on a posterior.

        Parameters
        ----------------
        in_percentage: int, default=10
             Percentage of posterior to burn-in.
        front_percentage: float, default=100
            Percentage of posterior to keep at the front.
        sample_name: str, default='draw'
            Name of dimension of xml_set.

        """
        if not isinstance(in_percentage, (int,float)) or in_percentage < 0 or in_percentage > 100:
            raise TypeError('in_percentage must be an integer or a float between 0 and 100.')
        if not isinstance(front_percentage, (int,float)) or front_percentage < 0 or front_percentage > 100:
            raise TypeError('front_percentage must be an integer or a float between 0 and 100.')
        if front_percentage < in_percentage:
            raise ValueError('front_percentage must be greater than or equal to in_percentage.')
        self.burinin_percentage = in_percentage
        self.keep_front_percentage = front_percentage
        selected_posterior = burn_posterior(posterior=self.original_posterior,
                                              in_percentage=self.burinin_percentage,
                                              front_percentage=front_percentage,
                                              sample_name=sample_name)
        self.selected_posterior = select_chains_from_posterior(selected_posterior,
                                                          self.selected_chains)


    def select_chains(self, chains=None, **kwargs):
        """
        Select chains to use in posterior.

        Parameters
        ----------------
        chains: list of strings or ints, default=None
            List of chains to select from posterior.
        kwargs: dict, default=None
            Chains labels with bool value. If True chain is included in posterior.
            If False chain is excluded from posterior.

        """
        if chains is None and not kwargs:
            raise ValueError("Either chains or kwargs must be provided")
        if chains is not None and kwargs:
            raise ValueError("Either chains or kwargs must be provided")

        if kwargs:
            for chain, value in kwargs.items():
                if chain not in self.original_chains:
                    raise ValueError("Chain {0} not found in original posterior, ".format(chain) +
                                     'see self.original_chains.')
                if not isinstance(value, bool):
                    raise TypeError("Chain {0} was not given a boolean value.".format(chain))
            for chain in self.original_chains:
                if chain not in kwargs:
                    raise ValueError("Chain {0} not found in kwargs, ".format(chain))
            self.selected_chains = [chain for chain in self.original_chains if kwargs[chain]]

        if chains is not None:
            for chain in chains:
                if chain not in self.original_chains:
                    raise ValueError("Chain {0} not found in original posterior, ".format(chain) +
                                     'see self.original_chains.')
            self.selected_chains = chains

        selected_posterior = burn_posterior(posterior=self.original_posterior,
                                              in_percentage=self.burinin_percentage,
                                              front_percentage=self.keep_front_percentage)
        self.selected_posterior = select_chains_from_posterior(selected_posterior, self.selected_chains)


    def diagnose_selection(self):
        """
        Diagnose selected_posterior.

        Updates self.diagnosis_of_selection via
        arviz.summary(self.selected_posterior, kind='diagnostics')

        """
        self.diagnosis_of_selection = az.summary(self.selected_posterior,
                                                                  kind='diagnostics')

    def merge_logs_to_csv(self, output_file='merged_logs.csv', like_logcombiner=True):
        """
        Merge selected log files into one csv file.

        Parameters
        ----------
        output_file : str, default='merged.log'
            Name of output file. Saved in self.directory.
        like_logcombiner : bool, default=True
            Merged output looks similar to merged output from logcombiner.

        """
        merge_logs_to_csv(posterior=self.selected_posterior, output_file=output_file, like_logcombiner=like_logcombiner)


    def _widget_interaction(self, percentages, parameters, **kwargs):
        for chain, value in kwargs.items():
            if chain not in self.original_chains:
                raise ValueError("Chain {0} not found in original posterior, ".format(chain) +
                                 'see self.original_chains.')
            if not isinstance(value, bool):
                raise TypeError("Chain {0} was not given a boolean value.".format(chain))
        for chain in self.original_chains:
            if chain not in kwargs:
                raise ValueError("Chain {0} not found in kwargs, ".format(chain))

        selected_chains = [chain for chain in self.original_chains if kwargs[chain]]
        posterior_modified = False
        if selected_chains != self.selected_chains:
            posterior_modified = True
            self.select_chains(selected_chains)
        in_percentage, front_percentage = percentages
        if in_percentage != self.burinin_percentage or front_percentage != self.keep_front_percentage:
            posterior_modified = True
            self.set_burnin(in_percentage=in_percentage, 
                            front_percentage=front_percentage)
        if posterior_modified:
            self.diagnose_selection()
        if len(selected_chains) > 0:
            stats, trace = self._display_diagnosis(parameters)
            plt.show()
            display(stats)


    @property
    def parameters_types(self):
        return self.parameters_dict.keys()

    def _display_diagnosis(self, parameters):
        if self.selected_posterior is None:
            raise AssertionError('BEASTDiag must be modified first. ' +
                                 'Use method set_burnin_and_chains().')
        if isinstance(parameters, str) and parameters in self.parameters_dict:
            parameters = self.parameters_dict[parameters]
        else:
            parameters = [parameters]
        for parameter in parameters:
            if parameter not in self.parameters:
                raise ValueError("Parameter {0} not found in posteriors list of parameters.".format(parameter))

        traces_fig, traces_ax = plot_traces(posterior=self.selected_posterior,
                                            parameters=parameters,
                                            labels=self.selected_chains)
        return self.diagnosis_of_selection.loc[parameters], traces_fig

    def clear_all_chains(self):
        """Clear all chains."""
        for chain in self.chain_checks:
            chain.value = False

    def select_all_chains(self):
        """Select all chains."""
        for chain in self.chain_checks:
            chain.value = True
    
    def _merge_selection(self, progress_callback=None, trees_pbar=None):
        merged_log_path = f'{self.output_prefix}merged_logs.csv'
        merged_trees_path = f'{self.output_prefix}merged_trees.trees'
        selection_yaml_path = f'{self.output_prefix}merge_selection.yaml'

        if progress_callback is not None:
            progress_callback(0, 3, 'Merging selected log files...')
        self.merge_logs_to_csv(output_file=merged_log_path, like_logcombiner=True)

        selected_nexus_files = sorted(            
            os.path.join(self.directory, fname)
            for fname in os.listdir(self.directory)
            if fname.endswith(".trees")
            and any(fname.startswith(chain) for chain in self.selected_chains)
        )
        if progress_callback is not None:
            progress_callback(1, 3, 'Subsetting and merging .trees files...')

        # Convert percentage-based selection to actual STATE numbers
        draws = self.original_posterior.posterior['draw'].values
        n_draws = len(draws)
        in_idx = round(self.burinin_percentage / 100 * n_draws)
        front_idx = round(self.keep_front_percentage / 100 * n_draws) - 1
        in_state_number = int(draws[min(in_idx, n_draws - 1)])
        front_state_number = int(draws[min(front_idx, n_draws - 1)])

        subset_and_merge_trees(
            file_list=selected_nexus_files,
            in_number=in_state_number,
            front_number=front_state_number,
            output_file=merged_trees_path,
            pbar=trees_pbar
        )

        selection_info = {
            'selected_chains': self.selected_chains,
            'burinin_percentage': self.burinin_percentage,
            'keep_front_percentage': self.keep_front_percentage
        }
        if progress_callback is not None:
            progress_callback(2, 3, 'Writing selection YAML...')
        with open(selection_yaml_path, 'w', encoding='utf-8') as file:
            yaml.safe_dump(selection_info, file, sort_keys=False)

        if progress_callback is not None:
            progress_callback(3, 3, 'Merge complete.')

        return {
            'merged_logs': merged_log_path,
            'merged_trees': merged_trees_path,
            'selection_yaml': selection_yaml_path
        }



    def generate_widget(self,
                 parameters_displayed=4):
        """
        Generates widget for selecting burn-in and chains, for use in Jupyter notebook.

        Parameters
        -------------
        parameters_displayed: int
            Number of parameters to display at a time.

        Returns
        -----------
        ipywidgets.widgets.VBox

        """
        num_params = len(self.parameters)
        parameters_dict = {}
        from_index = 0
        to_index = parameters_displayed
        parameter_set = 1
        while to_index <= num_params-1:
            parameters_dict[f'parameters set {str(parameter_set)}'] = self.parameters[from_index:to_index]
            parameter_set +=1
            from_index += parameters_displayed
            to_index += parameters_displayed
        parameters_dict[f'parameters set {str(parameter_set)}'] = self.parameters[from_index:]
        self.parameters_dict = parameters_dict
        burnin_title = widgets.HTML('Burn-in & Keep Front')
        burnin_selector = widgets.IntRangeSlider(
            value=(10,100),
            min=0,
            max=100,
            step=1,
            description='%:',
            disabled=False,
            continuous_update=False,
            orientation='horizontal',
            readout=True
        )
        chain_title = widgets.HTML('Chains')
        self.chain_checks = []
        for chain in self.original_chains:
            self.chain_checks.append(widgets.Checkbox(value=True, description=chain, disabled=False))
        chain_selector = widgets.GridBox(self.chain_checks,
                                         layout=widgets.Layout(grid_template_columns="repeat(4, 250px)"))
        # select_all_buttion = widgets.Button(description='Select all Chains',
        #                                  tooltip='Select all Chains',
        #                                  button_style='')
        # select_all_buttion.on_click(self.select_all_chains)
        # clear_all_buttion = widgets.Button(description='Clear all chains',
        #                                  tooltip='Clear all Chains',
        #                                  button_style='')
        # clear_all_buttion.on_click(self.clear_all_chains)
        # select_clear_all = widgets.HBox(children=[select_all_buttion, clear_all_buttion])
        burnin_and_chain_selector = widgets.VBox(children=[
            burnin_title,
            burnin_selector,
            chain_title,
            chain_selector,
            # select_clear_all
            ],
            titles=('Burn-in', 'Chains'))
        parameter_selector = widgets.Dropdown(options=self.parameters_types,
                                              description='Parameters:')
        output_widget = widgets.interactive_output(self._widget_interaction,
                                                   controls={
                                                       'percentages': burnin_selector,
                                                       'parameters': parameter_selector,
                                                       **{chain.description: chain for chain in self.chain_checks}
                                                   })
        merge_button = widgets.Button(
            description='Merge Selection',
            tooltip='Merge selection into one log (.csv) and one trees (.trees) file.',
            button_style='primary'
        )
        merge_progress = widgets.IntProgress(
            value=0,
            min=0,
            max=3,
            step=1,
            description='Merge:',
            bar_style='',
            orientation='horizontal'
        )
        merge_progress_status = widgets.HTML('')
        merge_status = widgets.HTML('')
        merge_tqdm_output = widgets.Output()

        def _on_merge_click(_):
            merge_progress.value = 0
            merge_progress.max = 3
            merge_progress.bar_style = 'info'
            merge_progress_status.value = "<span style='color: #1565c0;'>Starting merge...</span>"
            merge_status.value = ''

            with merge_tqdm_output:
                clear_output(wait=True)
                trees_pbar = tqdm(
                    total=len(self.selected_chains),
                    desc='Subsetting/merging .trees',
                    unit='file',
                    leave=False
                )

            def _update_progress(value, total, message):
                merge_progress.max = total
                merge_progress.value = value
                merge_progress_status.value = f"<span style='color: #1565c0;'>{message}</span>"

            try:
                output_paths = self._merge_selection(progress_callback=_update_progress,
                                                     trees_pbar=trees_pbar)
                merge_progress.bar_style = 'success'
                merge_status.value = (
                    "<span style='color: #2e7d32; font-weight: 600;'>"
                    "Merge complete.<br>"
                    f"Merged logs can be found at {output_paths['merged_logs']}.<br>"
                    f"Merged trees can be found at {output_paths['merged_trees']}.<br>"
                    f"Selection info can be found at {output_paths['selection_yaml']}."
                    "</span>"
                )
            except Exception as exc:
                merge_progress.bar_style = 'danger'
                merge_progress_status.value = "<span style='color: #c62828;'>Merge failed.</span>"
                merge_status.value = (
                    "<span style='color: #c62828; font-weight: 600;'>"
                    f"Merge failed: {type(exc).__name__}: {exc}"
                    "</span>"
                )
            finally:
                trees_pbar.close()

        merge_button.on_click(_on_merge_click)

        beast_diag_widget = widgets.VBox([
            burnin_and_chain_selector,
            parameter_selector,
            output_widget,
            merge_button,
            merge_progress,
            merge_progress_status,
            merge_tqdm_output,
            merge_status
        ])

        return beast_diag_widget