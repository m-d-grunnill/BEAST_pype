import os
import xarray as xr
import arviz as az
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display
import glob


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


def burnin_posterior(posterior,
                     proportion=None,
                     percentage=None,
                     number=None,
                     sample_name='draw'):
    """
    Perform a burn-in on a posterior.

    Parameters
    ----------
    posterior: arviz.data.inference_data.InferenceData
        DataArray with posterior. Must have 'chain' and 'draw' dimension names for in
         posterior.
    proportion: float, default=None
        Proportion of posterior to burn-in.
    percentage: float, default=None
         Percentage of posterior to burn-in.
    number: int
        Number of burn-in points.
    sample_name: str, default='draw'
        Name of dimension of xml_set.

    Returns
    -------
    arviz.data.inference_data.InferenceData
    """
    if proportion is None and number is None and percentage is None:
        raise ValueError("Either proportion, percentage or number must be provided")
    if proportion is not None:
        if not isinstance(proportion, float):
            raise TypeError("Proportion must be a float")
        if number is not None or percentage is not None:
            raise ValueError("One of proportion, percentage or number must be provided")
        number = round(proportion * len(posterior.posterior[sample_name]))
    elif percentage is not None:
        if not isinstance(percentage, (int,float)):
            raise TypeError("Percentage must be a float or an integer.")
        if number is not None or proportion is not None:
            raise ValueError("One of proportion, percentage or number must be provided")
        number = round(percentage/100 * len(posterior.posterior[sample_name]))
    elif isinstance(number, int):
        raise TypeError("Number must be an integer")
    selection = {sample_name: slice(number, None)}
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
    selected_chains: list of strs
        Chains in selected posterior.
    diagnostics_of_selection: pandas.DataFrame
        Summary of diagnostics of self.selected_posterior. Product of
        arviz.summary(self.selected_posterior, kind='diagnostics')

    """
    def __init__(self,
                 directory,
                 start=0,
                 cut_to_first=None,
                 ignore_iqtree_logfiles=True):
        """

        Parameters
        ----------
        directory:  str
            Directory where log and trees files are located.
        start : int
            Starting index of the first log file to read.
                sample_name: str, default='Sample'
        cut_to_first : int, default None
            Remove Samples/links over this number in log files.
        """
        self.directory = directory
        if not os.path.exists(directory):
            raise FileNotFoundError('directory does not exist')
            
        if ignore_iqtree_logfiles:
            log_paths = {
                path.removeprefix(directory + '/').removesuffix('.log'): path
                for path in glob.glob(directory + '/*.log')
                if not path.endswith('iqtree.log')
            }
        else:
            log_paths = {
                path.removeprefix(directory + '/').removesuffix('.log'): path
                for path in glob.glob(directory + '/*.log')
            }
        chains = list(log_paths.keys())
        self.original_posterior = read_log_files_as_posterior(log_paths,
                                                              start=start,
                                                              cut_to_first=cut_to_first)

        self.parameters = list(self.original_posterior.posterior.data_vars)
        self.parameters_dict = None
        self.original_chains = sorted(chains)
        self.selected_posterior = None
        self.burinin_percentage = 0
        self.selected_chains = chains
        self.diagnostics_of_selection = None
        self.cut_to_first = cut_to_first

    def set_burnin(self, percentage=10, sample_name='draw'):
        """
        Perform a burn-in on a posterior.

        Parameters
        ----------------
        percentage: int, default=10
             Percentage of posterior to burn-in.
        sample_name: str, default='draw'
            Name of dimension of xml_set.

        """
        if not isinstance(percentage, int):
            raise TypeError('percentage must be an integer')
        self.burinin_percentage = percentage
        selected_posterior = burnin_posterior(posterior=self.original_posterior,
                                              percentage=self.burinin_percentage,
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

        selected_posterior = burnin_posterior(posterior=self.original_posterior,
                                              percentage=self.burinin_percentage)
        self.selected_posterior = select_chains_from_posterior(selected_posterior, self.selected_chains)


    def diagnose_selection(self):
        """
        Diagnose selected_posterior.

        Updates self.diagnosis_of_selection via
        arviz.summary(self.selected_posterior, kind='diagnostics')

        """
        self.diagnosis_of_selection = az.summary(self.selected_posterior,
                                                                  kind='diagnostics')

    def merge_logs_to_csv(self, output_file='merged_log.csv', like_logcombiner=True):
        """
        Merge selected log files into one csv file.

        Parameters
        ----------
        output_file : str, default='merged.log'
            Name of output file. Saved in self.directory.
        like_logcombiner : bool, default=True
            Merged output looks similar to merged output from logcombiner.

        """
        df = self.selected_posterior.to_dataframe()
        if like_logcombiner:
            step = df.iloc[1, 1] - df.iloc[0, 1]
            df[df.columns[1]] = range(0, len(df) * step, step)
            if df.columns[1] != 'Sample':
                relabel_dict = {df.columns[1]: 'Sample'}
                df.rename(columns=relabel_dict, inplace=True)
            df = df.drop(columns=['chain'])
        df.to_csv(f'{self.directory}/{output_file}', index=False)

    def logcombiner_args(self, output_pretfix, output_file='merged', suffix='.log'):
        """
        Create logcombiner args for merging selected chain files (.log or .trees) into one.

        Parameters
        ----------
        output_pretfix: str
            Prefix for output file. Saved in self.directory.
        output_file : str, default='merged'
            Name of output file. Saved in self.directory.
        suffix: str, default '.log'
            Suffix of files to merge.

        Returns
        -------
        lc_args: dict {str: str}
            Logcombiner args.
        """
        if self.cut_to_first is not None:
            raise AttributeError(
                "Use of cut_to_first and BEASTs logcombiner is not yet implemented." +
                'However, merge_logs_to_csv method of this class can perform your' +
                ' selection and produce a merged csv from the log files.')
        if suffix not in ['.log', '.trees']:
            raise ValueError("suffix must be either '.log', '.trees'.")
        output_file = f'{output_pretfix}{output_file}{suffix}'
        selected_files = [f"{self.directory}/{chain}{suffix}"
                              for chain in self.selected_chains]
        name_pre_fix = suffix.replace('.', '')
        lc_args = {
            f'{name_pre_fix}_file_burnin':  self.burinin_percentage,
            f'{name_pre_fix}_file_output': output_file,
            f'{name_pre_fix}_files_to_combine': ' '.join(selected_files)
        }
        return lc_args

    def merging_outputs_params(self, output_path, kernel_name = 'beast_pype', xml_set=None):
        if xml_set is None:
            output_prefix = f'{output_path}/'
        else:
            output_prefix = f'{output_path}/{xml_set}_'
        params = {
            **self.logcombiner_args(output_prefix, suffix='.log'),
            **self.logcombiner_args(output_prefix, suffix='.trees')
        }
        params['kernel_name'] = kernel_name
        return params


    def _widget_interaction(self, percentage, parameters, **kwargs):
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
        if percentage != self.burinin_percentage:
            posterior_modified = True
            self.set_burnin(percentage=percentage)
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
        burnin_title = widgets.HTML('Burn-in')
        burnin_selector = widgets.IntSlider(
            value=10,
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
                                                       'percentage': burnin_selector,
                                                       'parameters': parameter_selector,
                                                       **{chain.description: chain for chain in self.chain_checks}
                                                   })
        beast_diag_widget = widgets.VBox([
            burnin_and_chain_selector,
            parameter_selector,
            output_widget])

        return beast_diag_widget