import xarray as xr
import arviz as az
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display
import glob
import os
import nbformat as nbf
import importlib.resources as importlib_resources

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
    Perform a burnin on a posterior.

    Parameters
    ----------
    posterior: arviz.data.inference_data.InferenceData
        DataArray with posterior. Must have 'chain' and 'draw' dimension names for in
         posterior.
    proportion: float, default=None
        Proportion of posterior to burnin.
    percentage: float, default=None
         Percentage of posterior to burnin.
    number: int
        Number of burnin points.
    sample_name: str, default='draw'
        Name of dimension of sample.

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
        axs[i][0].set_title('')
        axs[i][1].set_title('')
        axs[i][0].set_ylabel(parameter)

    return fig, axs


class BEASTDiag:
    """
    Class to perform a burnin and chain selection on a sampled posterior from BEAST.

    Attributes
    ---------------
    directory:  str
        Directory where log and trees files are located.
    original_posterior: arviz.data.inference_data.InferenceData
        Posterior with no chains removed and no burnin.
    parameters: list of str
        List of parameters names.
    parameters_dict: dict {str: [strs]}
        Dictionary of parameters categories.
    original_chains: list of strs
         Chains in original posterior.
    selected_posterior: arviz.data.inference_data.InferenceData
        Posterior with chains removed and burnin applied.
    burinin_percentage: int
        Percentage burnin applied to self.selected_posterior.
    selected_chains: list of strs
        Chains in selected posterior.
    diagnostics_of_selection: pandas.DataFrame
        Summary of diagnostics of self.selected_posterior. Product of
        arviz.summary(self.selected_posterior, kind='diagnostics')

    """

    def __init__(self,
                 directory,
                 parameters_dict=None,
                 start=0,
                 cut_to_first=None,
                 ignore_iqtree_logfiles=True):
        """

        Parameters
        ----------
        directory:  str
            Directory where log and trees files are located.
        parameters_dict: dict {str: [strs]}
            Dictionary of parameters categories.
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
        if parameters_dict is None:
            num_params = len(self.parameters)
            parameters_dict = {}
            from_index = 0
            to_index = 5
            parameter_set = 1
            while to_index <= num_params-1:
                parameters_dict[f'parameters set {str(parameter_set)}'] = self.parameters[from_index:to_index]
                parameter_set +=1
                from_index += 5
                to_index += 5
            parameters_dict[f'parameters set {str(parameter_set)}'] = self.parameters[from_index:]
        self.parameters_dict = parameters_dict
        self.original_chains = sorted(chains)
        self.selected_posterior = None
        self.burinin_percentage = 0
        self.selected_chains = chains
        self.diagnostics_of_selection = None
        self.cut_to_first = cut_to_first

    def set_burnin(self, percentage=10, sample_name='draw'):
        """
        Perform a burnin on a posterior.

        Parameters
        ----------------
        percentage: int, default=10
             Percentage of posterior to burnin.
        sample_name: str, default='draw'
            Name of dimension of sample.

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

    def logcombiner_args(self, output_file='merged', suffix='.log'):
        """
        Create logcombiner args for merging selected chain files (.log or .trees) into one.

        Parameters
        ----------
        output_file : str, default='merged' + suffix
            Name of output file. Saved in self.directory.
        suffix: str, default '.log'
            Suffix of files to merge.

        Returns
        -------
        lc_args: str
            Logcombiner args to be passed to bash cell. 
            $1: str(self.burinin_percentage)
            $2: Merged output_file + suffix.
            ${@:3} Files to be merged.

        """
        if self.cut_to_first is not None:
            raise AttributeError(
                "Use of cut_to_first and BEASTs logcombiner is not yet implemented." +
                'However, merge_logs_to_csv method of this class can perform your' +
                ' selection and produce a merged csv from the log files.')
        if suffix not in ['.log', '.trees']:
            raise ValueError("suffix must be either '.log', '.trees'.")
        output_file = f'{self.directory}/{output_file}{suffix}'
        selected_files = [f"{self.directory}/{chain}{suffix}"
                              for chain in self.selected_chains]
        lc_args = ' '.join([str(self.burinin_percentage), output_file, *selected_files])
        return lc_args

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

    def generate_widget(self):
        """
        Generates widget for selecting burnin and chains, for use in Jupyter notebook.

        Returns
        -----------
        ipywidgets.widgets.VBox

        """
        burnin_title = widgets.HTML('Burnin')
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
        chain_checks = [
            widgets.Checkbox(value=True, description=chain, disabled=False)
            for chain in self.original_chains]
        chain_selector = widgets.HBox(children=chain_checks)
        burnin_and_chain_selector = widgets.VBox(children=[
            burnin_title,
            burnin_selector,
            chain_title,
            chain_selector],
            titles=('Burnin', 'Chains'))
        parameter_selector = widgets.Dropdown(options=self.parameters_types,
                                              description='Parameters:')
        output_widget = widgets.interactive_output(self._widget_interaction,
                                                   controls={
                                                       'percentage': burnin_selector,
                                                       'parameters': parameter_selector,
                                                       **{chain.description: chain for chain in chain_checks}
                                                   })
        beast_diag_widget = widgets.VBox([
            burnin_and_chain_selector,
            parameter_selector,
            output_widget])

        return beast_diag_widget


class BDSKYSurveillanceDiag(BEASTDiag):
    """
    A version of BEASTDiag for use with BDSKY models (preset parameters_dict).

    Attributes
    ---------------
    parameters_dict: dict {str: [strs]}
        Dictionary of parameters categories tailored for analyses of BEAST runs using
         the BDSKY model.

    """

    def __init__(self, directory, start=0, sample_name='Sample', cut_to_first=None):
        """

        Parameters
        ----------
        directory:  str
            Directory where log and trees files are located.
        start : int
            Starting index of the first log file to read.
                sample_name: str, default='Sample'
        sample_name: str, default='Sample'
            Name of the column for sample/link/iteration. Will be relabelled 'draw' to
             conform with arviz api.
        cut_to_first : int, default None
            Remove Samples/links over this number in log files.
        """
        super().__init__(directory, parameters_dict=None,
                         start=start,
                         sample_name=sample_name,
                         cut_to_first=cut_to_first)
        self.parameters_dict = {
            'general': ['posterior', 'likelihood', 'prior', 'treeLikelihood', 'BDSKY_Serial'],
            'evolutionary': ['TreeHeight', 'clockRate', 'gammaShape', 'kappa',
                             'freqParameter.1', 'freqParameter.2',
                             'freqParameter.3', 'freqParameter.4'],
            'epidemiological': ['origin_BDSKY_Serial',
                                'becomeUninfectiousRate_BDSKY_Serial',
                                'samplingProportion_BDSKY_Serial'],
            'R_e': [param for param in self.parameters
                    if param.startswith('reproductiveNumber_BDSKY_Serial.')]
        }

def gen_xml_set_diag_notebook(save_dir, directories_to_exclude = ['cache']):
    """
    Generates rest of Phase-5-Diagnosing_XML_sets_and_Generate_Report.ipynb.



    Parameters
    ----------
    save_dir: str
        Path to save copy of Phase-5-Diagnosing_XML_sets_and_Generate_Report.ipynb.
        Each subdirectory of save_dir will be assumed to contain data of an xml_set
         and used in generating thes new
         Phase-5-Diagnosing_XML_sets_and_Generate_Report.ipynb.
    directories_to_exclude:
        Any directories to exclude from save_dir when generating new
        Phase-5-Diagnosing_XML_sets_and_Generate_Report.ipynb.

    Returns
    -------
    None
    """
    directories = [directory for directory in os.listdir(save_dir)
                   if os.path.isdir(os.path.join(save_dir, directory)) and
                   directory not in directories_to_exclude]
    workflows_modules = importlib_resources.path('beast_pype', 'workflow_modules')
    diag_notebook_path = f'{workflows_modules}/Phase-5-Diagnosing_XML_sets_and_Generate_Report.ipynb'
    diag_notebook = nbf.read(diag_notebook_path, as_version=4)
    for xml_set in directories:
        diag_notebook['cells'] += [
            nbf.v4.new_markdown_cell(
                f"## Diagnosing XML Set: {xml_set}\n" +
                "### Loading data"),
            nbf.v4.new_code_cell(f"sample_diag = BEASTDiag('{xml_set}')"),
            nbf.v4.new_markdown_cell(
                "## Selecting burnin and Chains to Remove\n\n" +
                "Activating the cell below will generate an interactive widget. Widgets parts:\n" +
                "* Top interactive part: this allows you to select for a different burnin and remove chains and select the parameters used in the rest of the widget.,\n" +
                "* Middle display: KDE and trace plots, see [arviz.plot_trace documentation](https://python.arviz.org/en/stable/api/generated/arviz.plot_trace.html#arviz.plot_trace).\n" +
                "* Bottom display: A table of statistics regarding the traces, see [arviz.summary documentation](https://python.arviz.org/en/stable/api/generated/arviz.summary.html#arviz.summary). Regarding these statistics:\n" +
                "\t* Ideally the ESSs should be >= 200, see [arviz.ess documentation](https://python.arviz.org/en/stable/api/generated/arviz.ess.html#arviz.ess).\n" +
                "\t* Ideally the r_hat should be close fo 1, see [arviz.rhat documentation](https://python.arviz.org/en/stable/api/generated/arviz.rhat.html#arviz.rhat)).\n" +
                "\t* Markov Chain Standard Error MCSEs, see [arviz.mcse](https://python.arviz.org/en/stable/api/generated/arviz.mcse.html#arviz.mcse).\n\n"+
                "After making your selection click on the cell below the widget and then keep pressing shift+enter to carry on with the rest of the cells in this notebook."
            ),
            nbf.v4.new_code_cell(
                "sample_diag_widget = sample_diag.generate_widget()\n" +
                "sample_diag_widget"),
            nbf.v4.new_markdown_cell(
                "### Merge Kept chains\n" +
                "#### Log files"),
            nbf.v4.new_code_cell(
                "%%bash -l -s {sample_diag.logcombiner_args(suffix='.log')}\n" +
                "source activate beast_pype\n\n" +
                "logcombiner -b $1 -log ${@:3} -o $2"),
            nbf.v4.new_markdown_cell(
                "#### Tree files"),
            nbf.v4.new_code_cell(
                "%%bash -l -s {sample_diag.logcombiner_args(suffix='.trees')}\n" +
                "source activate beast_pype\n\n" +
                "logcombiner -b $1 -log ${@:3} -o $2"),
            nbf.v4.new_markdown_cell(
                "#### Record Chains used & burnin for this sample"),
            nbf.v4.new_code_cell(
                f"pipeline_run_info['Chains Used']['{xml_set}'] = deepcopy(sample_diag.selected_chains)\n" +
                f"pipeline_run_info['Burn-In']['{xml_set}'] = deepcopy(sample_diag.burinin_percentage)")
        ]

    diag_notebook['cells'] += [
        nbf.v4.new_markdown_cell(
            "## Update the pipeline_run_info json."),
        nbf.v4.new_code_cell(
            "with open(f'{save_dir}/pipeline_run_info.json', 'w') as fp:\n\tjson.dump(pipeline_run_info, fp, sort_keys=True, indent=4)\nfp.close()"
        ),
        nbf.v4.new_markdown_cell(
            '## Generate output Report\n' +
            'Now you can now move on to visualising outputs from BEAST using a report template.'
        ),
        nbf.v4.new_code_cell(
            "report_params = {'save_dir': save_dir, 'collection_date_field': collection_date_field}\n" +
            "output_report_path = f'{save_dir}/BEAST_pype-Report.ipynb'\n" +
            "if add_unreported_fields:\n" +
            "\tadd_unreported_outputs(report_template, f'{sample_diag.directory}/merged.log', output_report_path, xml_set_comparisons=True)\n" +
            "\tinput_path = output_report_path\n" +
            "else:\n" +
            "\tinput_path = report_template\n\n" +
            "output = pm.execute_notebook(\n" +
            "\tinput_path=input_path,\n" +
            "\toutput_path=output_report_path,\n" +
            "\tparameters=report_params,\n" +
            "\tprogress_bar=True)"),
        nbf.v4.new_markdown_cell(
            "### Convert Output Report from Jupyter Notebook to Notebook\n" +
            "This also removes code cells."
        ),
        nbf.v4.new_code_cell(
            "%%bash -l -s {output_report_path}\n" +
            "source activate beast_pype\n" +
            "jupyter nbconvert --to html --no-input $@"
        )
    ]

    with open(f'{save_dir}/Phase-5-Diagnosing_XML_sets_and_Generate_Report.ipynb', 'w') as f:
        nbf.write(diag_notebook, f)



