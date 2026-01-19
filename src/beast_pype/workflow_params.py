import time
import os
from copy import deepcopy
import yaml
from papermill.inspection import _infer_parameters
from papermill.iorw import load_notebook_node
from numpy.random import randint
import multiprocessing
import pandas as pd
from Bio import SeqIO
from datetime import datetime
from types import SimpleNamespace
from pandas.tseries.offsets import DateOffset
import importlib.resources as importlib_resources

report_templates_path = importlib_resources.path('beast_pype', 'report_templates')
workflows_path = importlib_resources.path('beast_pype', 'workflows')

def _gen_phase_4_params(
        save_dir,
        seeds,
        beast_options_without_a_value=None,
        beast_options_needing_a_value=None,
        max_threads=None,
        sbatch_options_without_a_value=None,
        sbatch_options_needing_a_value=None):
    """
    Generate phase 4 parameters.

    Parameters
    ----------
    save_dir: str
        Path to directory where workflow outputs are saved.
    seeds: list of ints or numpy.ndarray of ints
        Seeds used when running BEAST 2 in phase 4.
    beast_options_without_a_value: list of strs
        Single word arguments to pass to BEAST 2.
         For instance to use a GPU this would be ['-beagle_GPU'].
        See https://www.beast2.org/2021/03/31/command-line-options.html.
    beast_options_needing_a_value: dict
        Word followed by value arguments to pass to BEAST 2.
        For instance to use 3 threads when running BEAST 2 this would be: {'-threads': 3}.
        See https://www.beast2.org/2021/03/31/command-line-options.html.
    sbatch_options_without_a_value: list of strs
        Single word arguments to pass to sbatch.
        See https://slurm.schedmd.com/sbatch.html.
    sbatch_options_needing_a_value: dict
        Word followed by value arguments to pass to sbatch.
        See https://slurm.schedmd.com/sbatch.html.
    max_threads: int, optional
        Maximum number of threads to use.

    Returns
    -------
    Dictionary of phase 4 parameters.
    """
    if beast_options_without_a_value is None:
        beast_options_without_a_value = []
    if beast_options_needing_a_value is None:
        beast_options_needing_a_value = {}
    params_dict = {
        'save_dir': save_dir,
        'seeds': seeds}
    if sbatch_options_without_a_value is not None or sbatch_options_needing_a_value is not None:
        params_dict['sbatch_arg_string'] = _to_arg_string(
            sbatch_options_without_a_value,
            sbatch_options_needing_a_value)
    else:
        if '-threads' in beast_options_needing_a_value:
            threads_per_run = beast_options_needing_a_value['-threads']
        else:
            threads_per_run = 1
        params_dict['jobs_per_time'] = int(max_threads / threads_per_run)
    params_dict['beast_arg_string'] = _to_arg_string(
            beast_options_without_a_value,
            beast_options_needing_a_value)
    return params_dict

class WorkflowParams(SimpleNamespace):
    """Base Class responsible for setting up and checking workflow parameters."""

    workflow_name = None


    def __init__(self, **kwargs):
        """
        Initialise the workflow parameters (perform checks and set up).

        Parameters
        ----------
        kwargs:
            Parameters to use in the workflow.
        """
        kwargs = {key: value for key, value in kwargs.items() if key not in ['sys', 'NamespaceMagics', 'get_ipython', 'debugpy', 'json'] and not 'pydev' in key}
        super().__init__(**kwargs)
        # Check all the parameters are valid for this workflow.
        invalid_params = set(kwargs.keys()) - self.accepted_param_names()
        if invalid_params:
            raise ValueError(f'The following parameters that have been supplied to the {self.workflow_name} workflow are invalid:\n {", ".join(invalid_params)}')

        # Check for not being assigned
        if self.overall_save_dir is None:
            raise Exception(
                "overall_save_dir is missing from the parameters yml file, or\n" +
                "overall_save_dir has been given an 'null' value (which are converted to None in python).\n" +
                "None/null values cannot be used for overall_save_dir."
            )

        if self.number_of_beast_runs is None:
            raise Exception(
                "number_of_beast_runs is missing from the parameters yml file, or\n" +
                "number_of_beast_runs has been given an 'null' value (which are converted to None in python).\n" +
                "None/null values cannot be used for number_of_beast_runs."
            )
        if not isinstance(self.use_initial_tree, bool):
            raise TypeError('use_initial_tree must be True or False')
        if ' ' in self.overall_save_dir:
            raise ValueError('overall_save_dir cannot contain whitespace.')
        if self.fasta_path is not None and ' ' in self.fasta_path:
            raise ValueError('fasta_path cannot contain whitespace.')
        if self.log_file_basename is not None and ' ' in self.log_file_basename:
            raise ValueError('log_file_basename cannot contain whitespace.')
        if self.metadata_path is not None and not self.metadata_path.endswith(('.tsv', '.csv')):
            raise ValueError('Only .tsv and .csv are supported for metadata. Use correct file extension.')

        # Check values.
        if self.metadata_path is not None:
            if self.metadata_path.endswith('.tsv'):
                metadata_df = pd.read_csv(self.metadata_path, sep='\t')
            else:
                metadata_df = pd.read_csv(self.metadata_path)
            ids_in_metadata = metadata_df[self.sample_id_field]
            if ids_in_metadata.duplicated().any():
                duplicates_str = ', '.join(ids_in_metadata[ids_in_metadata.duplicated()].to_list())
                raise ValueError(f'There are duplicated values in the metadata field {self.sample_id_field}.\n.' +
                                 f' Check values {duplicates_str}.')
            if self.fasta_path is not None:
                sequence_lengths = [(record.id, len(record.seq)) for record in SeqIO.parse(self.fasta_path, 'fasta')]
                sequence_unique_lengths = set(record[1] for record in sequence_lengths)
                if len(sequence_unique_lengths) != 1:
                    raise ValueError(f'Fasta file {self.fasta_path} does not contain aligned fasta_path.')
                ids_in_fasta = pd.Series(record[0] for record in sequence_lengths)
                if ids_in_fasta.duplicated().any():
                    duplicates_str = ', '.join(ids_in_fasta[ids_in_fasta.duplicated()].to_list())
                    raise ValueError(f'There are duplicated id values in the fasta file.\n.' +
                                     f' Check values {duplicates_str}.')
                ids_in_metadata = set(ids_in_metadata)
                ids_in_fasta = set(ids_in_fasta)
                in_fasta_not_in_metadata = ids_in_fasta - ids_in_metadata
                in_metadata_not_in_fasta = ids_in_metadata - ids_in_fasta
                if len(in_fasta_not_in_metadata) > 0 or len(in_metadata_not_in_fasta) > 0:
                    in_fasta_not_in_metadata = ', '.join(in_fasta_not_in_metadata)
                    in_metadata_not_in_fasta = ', '.join(in_metadata_not_in_fasta)
                    raise ValueError(
                        f'There is a mismatch between the ids in the fasta at fasta_path and metadata field {self.sample_id_field}.\n' +
                        f"The following are in the metadata but not the fasta file: \n{in_metadata_not_in_fasta}.\n" +
                        f"The following are in the fasta file but not the metadata: \n{in_fasta_not_in_metadata}."
                    )
                try:
                    pd.to_datetime(metadata_df[self.collection_date_field], errors='raise')
                except:
                    raise ValueError(
                        f'The field {self.collection_date_field} in {self.metadata_path} is cannot be read as a date.\n' +
                        f' Check the value in {self.collection_date_field}')
                self.number_of_sequences = len(ids_in_fasta)

        if self.root_strain_names is None and self.remove_root:
            raise ValueError('Root strain names must be specified if remove_root is True.')

        if self.max_threads is None:
            if 'SLURM_CPUS_PER_TASK' in os.environ:
                self.max_threads = int(os.environ['SLURM_CPUS_PER_TASK'])
            else:
                self.max_threads = multiprocessing.cpu_count() - 1

        # Checking MCMC options.
        if self.chain_length is not None:
            if not (isinstance(self.chain_length, int) and self.chain_length > 1):
                raise ValueError('If specified chain_length must be an integer greater than 1.')
        if self.trace_log_every is not None:
            if isinstance(self.trace_log_every, float) and self.trace_log_every < 1 and self.trace_log_every > 0:
                if not isinstance(self.chain_length, int):
                    raise AssertionError('If trace_log_every is a float, chain_length must be specified.')
                self.trace_log_every = int(round(self.trace_log_every*self.chain_length))
            elif not (isinstance(self.tree_log_every, int) and self.tree_log_every > 1):
                raise AssertionError('tree_log_every must be an integer greater than 1 or a float between 0 and 1.')
        if self.tree_log_every is not None:
            if isinstance(self.tree_log_every, float) and self.tree_log_every < 1 and self.tree_log_every > 0:
                if not isinstance(self.chain_length, int):
                    raise AssertionError('If tree_log_every is a float, chain_length must be specified.')
                self.tree_log_every = int(round(self.tree_log_every*self.chain_length))
            elif not (isinstance(self.tree_log_every, int) and self.tree_log_every > 1):
                raise AssertionError('tree_log_every must be an integer greater than 1 or a float between 0 and 1.')
        if self.screen_log_every is not None:
            if isinstance(self.screen_log_every, float) and self.screen_log_every < 1 and self.screen_log_every > 0:
                if not isinstance(self.chain_length, int):
                    raise AssertionError('If screen_log_every is a float, chain_length must be specified.')
                self.screen_log_every = int(round(self.screen_log_every*self.chain_length))
            elif not (isinstance(self.screen_log_every, int) and self.screen_log_every > 1):
                raise AssertionError('screen_log_every must be an integer greater than 1 or a float between 0 and 1.')
        if self.store_state_every is not None:
            if isinstance(self.store_state_every, float) and self.store_state_every < 1 and self.store_state_every > 0:
                if not isinstance(self.chain_length, int):
                    raise AssertionError('If store_state_every is a float, chain_length must be specified.')
                self.store_state_every = int(round(self.store_state_every*self.chain_length))
            elif not (isinstance(self.store_state_every, int) and self.store_state_every > 1):
                raise AssertionError('store_state_every must be an integer greater than 1 or a float between 0 and 1.')

        self.save_dir = self.overall_save_dir + '/' + self.specific_run_save_dir
        if self.metadata_path is not None:
            metadata_df.to_csv(f'{self.save_dir}/metadata.csv', index=False)

    def record_parameters(self, extras_dict={}):
        """
        Record the parameters used in workflow in `pipeline_run_info.yml`.

        Parameters
        ----------
        extras_dict: dict
            Dictionary of extra parameters to be recorded.
        """
        pipeline_run_info_path = f'{self.save_dir}/pipeline_run_info.yml'
        if os.path.exists(pipeline_run_info_path):
            with open(pipeline_run_info_path, 'r') as fp:
                to_record  = yaml.safe_load(fp)
            fp.close()
            if 'parameters' in to_record:
                to_record['parameters'].update(self.__dict__)
            to_record.update(extras_dict)
        else:
            to_record = {'parameters': self.__dict__, **extras_dict}
        with open(pipeline_run_info_path, 'w') as fp:
            yaml.dump(to_record, fp, sort_keys=True)
        fp.close()

    def retrieve_params(self, parameter_names):
        """
        Retrieve parameters held in an instance of this class.

        Parameters
        ----------
        parameter_names: list of strings
            List of parameter names to retrieve.

        Returns
        -------
        Dictionary of parameter names and values.
        """
        return {key: value for key, value in self.__dict__.items() if key in parameter_names}

    def retrieve_phase_4_params(self):
        """
        Retrieve parameters used in phase 4 of the workflow.

        Returns
        -------
        Dictionary of parameter names and values.
        """
        phase_4_params_dict = _gen_phase_4_params(**self.retrieve_params(['save_dir',
                                                                          'seeds',
                                                                          'beast_options_without_a_value',
                                                                          'beast_options_needing_a_value',
                                                                          'max_threads',
                                                                          'sbatch_options_without_a_value',
                                                                          'sbatch_options_needing_a_value']))
        phase_4_params_dict.update(self.retrieve_params(['kernel_name']))
        return phase_4_params_dict

    def accepted_param_names(self):
        """
        Returns a set of accepted parameter names.

        Returns
        -------
        Set of strings.
        """
        nb = load_notebook_node(self.associated_workflow)
        parameters = _infer_parameters(nb)
        parameter_names = [parameter.name
                           for parameter in parameters]  + ['parent_workflow_parameters']
        return set(parameter_names)


class SimpleWorkflowParams(WorkflowParams):
    """Setup and check parameters for a simple workflows. """


    def __init__(self, **kwargs):
        """
        Initialise the workflow parameters (perform checks and set up).

        Parameters
        ----------
        kwargs:
            Parameters to use in the workflow.
        """
        super().__init__(**kwargs)

        # Check for not being assigned
        if self.ready_to_go_xml is None:
            if self.fasta_path is None:
                raise Exception(
                    "A ready_to_go_xml has not being given and" +
                    "fasta_path is missing from the parameters yml file, or\n" +
                    "fasta_path has been given an 'null' value (which are converted to None in python)."
                )
            if self.template_xml_path is None:
                raise Exception(
                    "A ready_to_go_xml has not being given and" +
                    "template_xml_path is missing from the parameters yml file, or\n" +
                    "template_xml_path has been given an 'null' value (which are converted to None in python)."
                )
            if self.log_file_basename is None:
                raise Exception(
                    "A ready_to_go_xml has not being given and" +
                    "log_file_basename is missing from the parameters yml file, or\n" +
                    "log_file_basename has been given an 'null' value (which are converted to None in python)."
                )
            if self.metadata_path is None:
                raise Exception(
                    "A ready_to_go_xml has not being given and" +
                    "metadata_path is missing from the parameters yml file, or\n" +
                    "metadata_path has been given an 'null' value (which are converted to None in python).\n" +
                    "None/null values cannot be used for metadata_path."
                )
        else:
            if self.fasta_path is not None:
                raise Exception(
                    "A ready_to_go_xml has being given and:"
                    "fasta_path is in the parameters yml file, or\n" +
                    "fasta_path has not been given an 'null' value (which are converted to None in python)."
                )
            if self.template_xml_path is not None:
                raise Exception(
                    "A ready_to_go_xml has being given and:"
                    "template_xml_path is in the parameters yml file, or\n" +
                    "template_xml_path has not been given an 'null' value (which are converted to None in python)."
                )
            if self.log_file_basename is not None:
                raise Exception(
                    "A ready_to_go_xml has being given and:"
                    "log_file_basename is in the parameters yml file, or\n" +
                    "log_file_basename has not been given an 'null' value (which are converted to None in python)."
                )

        if self.ready_to_go_xml is not None or not self.use_initial_tree:  # None of the parameters below are needed if a ready_to_go_xml is provided.
            if self.initial_tree_type not in ['Temporal', 'Distance', None]:
                raise ValueError('initial_tree_type must be either "Temporal" or "Distance" or None.')

            if self.initial_tree_path is not None and self.initial_tree_type is None:
                raise ValueError('initial_tree_type must be specified if initial_tree_path is given.')

            if self.down_sample_to is not None and (self.initial_tree_path is not None or self.initial_tree_type == 'Distance'):
                raise ValueError("Currently beast_pype's down_sampling method is tied to its Tree Time tree building." +
                                 "Therefore, to use this down sampling method an initial_tree_path should not be given and an initial_tree_type should be set to 'Temporal'.")

        # If seeds not given generate them.
        # BEASTs random number seed can select the same seed for multiple runs if they are launched close together in time (such as programmatically). Therefore, numpy is used to generate seeds for running BEAST.
        if self.seeds is None:
            number_of_seeds = self.number_of_beast_runs
            self.seeds = randint(low=1, high=int(1e6), size=number_of_seeds).tolist()

    def retrieve_phase_2i_params(self):
        """
        Retrieve parameters used in phase 2i of the workflow.

        Returns
        -------
        Dictionary of parameter names and values.
        """
        return self.retrieve_params(['save_dir', 'fasta_path' , 'max_threads'])

    def retrieve_phase_2ii_params(self):
        """
        Retrieve parameters used in phase 2ii of the workflow.

        Returns
        -------
        Dictionary of parameter names and values.
        """
        return self.retrieve_params(['save_dir',
                                     'fasta_path',
                                     'metadata_path',
                                     'sample_id_field',
                                     'collection_date_field',
                                     'down_sample_to',
                                     'root_strain_names',
                                     'remove_root'])

    def retrieve_phase_5_params(self):
        """
        Retrieve parameters used in phase 5 of the workflow.

        Returns
        -------
        Dictionary of parameter names and values.
        """
        phase_5_params_dict = self.retrieve_params(['kernel_name'])
        if self.report_template is None:
            phase_5_params_dict['report_template'] = self.default_report_template
        else:
            phase_5_params_dict['report_template'] = self.report_template
        return phase_5_params_dict

    def tip_date_range(self):
        """
        Generate tip date range.
        Returns
        -------
        oldest_tip_date: datetime
            Oldest tip date.
        youngest_tip_date: datetime
            Youngest tip date.

        """
        if self.down_sample_to is None:
            metadata_path = self.metadata_path
        else:
            metadata_path = f'{self.save_dir}/down_sampled_metadata.csv'
        if metadata_path.endswith('.tsv'):
            sep = '\t'
        elif metadata_path.endswith('.csv'):
            sep = ','
        else:
            raise ValueError('Only .tsv and .csv are supported for metadata. Use correct file extension.')
        metadata = pd.read_csv(metadata_path, parse_dates=[self.collection_date_field], sep=sep)
        youngest_tip_date = metadata[self.collection_date_field].max()
        oldest_tip_date = metadata[self.collection_date_field].min()
        return oldest_tip_date, youngest_tip_date


class GenericWorkflowParams(SimpleWorkflowParams):
    """Set up and check Generic workflow parameters."""


    workflow_name = 'Generic'
    associated_workflow = f"{workflows_path}/Generic.ipynb"
    default_report_template = "Generic"

    def __init__(self, **kwargs):
        """
        Initialise the workflow parameters (perform checks and set up).

        Parameters
        ----------
        kwargs:
            Parameters to use in the workflow.
        """
        super().__init__(**kwargs)
        if self.ready_to_go_xml is not None and self.metadata_path is not None:
                raise Exception(
                    "A ready_to_go_xml has being given and"
                    "metadata_path is in the parameters yml file, or\n" +
                    "metadata_path has not been given an 'null' value (which are converted to None in python)."
                )
        self.record_parameters()
        self.sampling_prop_partition_dates = None
        self.rt_partition_dates = None

    def retrieve_phase_3_params(self):
        """
        Retrieve parameters used in phase 3 of the workflow.

        Returns
        -------
        Dictionary of parameter names and values.
        """
        return self.retrieve_params(['save_dir',
                                     'template_xml_path',
                                     'use_initial_tree',
                                     'initial_tree_path',
                                     'collection_date_field',
                                     'sample_id_field',
                                     'log_file_basename',
                                     'chain_length',
                                     'trace_log_every',
                                     'tree_log_every',
                                     'screen_log_every',
                                     'store_state_every'])

def _bdksy_serial_errors(rt_dims,
                         rt_partitions,
                         sampling_prop_dims,
                         sampling_prop_partitions,
                         zero_sampling_before_first_sample,
                         ready_to_go_xml,
                         use_initial_tree,
                         origin_start_addition,
                         initial_tree_type,
                         origin_upper_addition):
    if rt_dims is not None and rt_partitions is not None:
        raise TypeError('rt_dims and rt_partitions cannot be used together')
    if sampling_prop_dims is not None and sampling_prop_partitions is not None:
        raise TypeError('sampling_prop_dims and sampling_prop_partitions cannot be used together')
    if sampling_prop_dims is not None and zero_sampling_before_first_sample:
        raise TypeError('Currently zero_sampling_before_first_sample can only be used '+
                        'on its own or with sampling_prop_partitions, but NOT with sampling_prop_dims.')
    if ready_to_go_xml is not None or not use_initial_tree:
        if origin_start_addition is not None:
            if initial_tree_type != 'Temporal':
                raise ValueError('origin_start_addition is reliant on the initial_tree_type being "Temporal".')
            if origin_upper_addition is None:
                raise ValueError('origin_start_addition is reliant on origin_upper_addition being given as well.')

        if origin_upper_addition is not None:
            if initial_tree_type != 'Temporal':
                raise ValueError('origin_upper_addition is reliant on the initial_tree_type being "Temporal".')
            if origin_upper_addition is None:
                raise ValueError('origin_upper_addition is reliant on origin_start_addition being given as well.')

def _partition_dates_dict_to_list(change_dates_dict, youngest_tip_date, oldest_tip_date, parameter_name, inclusive_of_end=False):
    if 'end' in change_dates_dict:
        end = change_dates_dict['end']
        if isinstance(end, str):
            end = datetime.strptime(end, '%Y-%m-%d')
        if not isinstance(end, datetime):
            raise TypeError(f'{parameter_name}["end"] should be a datetime object or string of format YYYY-MM-DD')
    else:
        end = oldest_tip_date
    change_dates_list = []
    offsets = {change_dates_dict['unit']: change_dates_dict['every']}
    date_to_append = youngest_tip_date - DateOffset(**offsets)
    if not date_to_append > end:
        raise ValueError(
            f'Change dates dictionary for {parameter_name} makes no sense.\n' +
            'First partition date generated is not after the "end" date.\n' +
            'Suggest using smaller "units" or "every" values or pushing "end" date back.'
                         )
    while date_to_append > end:
        change_dates_list.append(date_to_append.strftime('%Y-%m-%d'))
        date_to_append = date_to_append - DateOffset(**offsets)
    if inclusive_of_end:
        change_dates_list.append(end.strftime('%Y-%m-%d'))
    return change_dates_list

def _phase_3_partition_params(rt_partitions,
                              sampling_prop_partitions,
                              sampling_prop_include_oldest_tip_date,
                              youngest_tip_date,
                              oldest_tip_date):
    if isinstance(rt_partitions, dict):
        rt_partition_dates = _partition_dates_dict_to_list(rt_partitions,
                                                        youngest_tip_date,
                                                        oldest_tip_date,
                                                        parameter_name='rt_partitions')
    else:
        rt_partition_dates = rt_partitions
    if isinstance(sampling_prop_partitions, dict):
        sampling_prop_partition_dates = _partition_dates_dict_to_list(sampling_prop_partitions,
                                                                   youngest_tip_date,
                                                                   oldest_tip_date,
                                                                   parameter_name='sampling_prop_partitions',
                                                                   inclusive_of_end=sampling_prop_include_oldest_tip_date
                                                                   )
    else:
        sampling_prop_partition_dates = sampling_prop_partitions
    return rt_partition_dates, sampling_prop_partition_dates

_accepted_bdsky_param = {'origin_start_addition',
                         'origin_upper_addition',
                         'origin_prior',
                         'rt_dims',
                         'rt_partitions',
                         'sampling_prop_dims',
                         'sampling_prop_partitions',
                         'zero_sampling_before_first_sample'
                         }



class BDSKYSerialWorkflowParams(SimpleWorkflowParams):
    """Set up and check BDSKY-Serial workflow parameters."""

    workflow_name = 'BDSKY-Serial'
    associated_workflow = f"{workflows_path}/BDSKY-Serial.ipynb"
    default_report_template = "BDSKY-Serial"

    def __init__(self, **kwargs):
        """
        Initialise the workflow parameters (perform checks and set up).

        Parameters
        ----------
        kwargs:
            Parameters to use in the workflow.
        """
        super().__init__(**kwargs)
        _bdksy_serial_errors(
            rt_dims=self.rt_dims,
            rt_partitions=self.rt_partitions,
            sampling_prop_dims=self.sampling_prop_dims,
            sampling_prop_partitions=self.sampling_prop_partitions,
            zero_sampling_before_first_sample=self.zero_sampling_before_first_sample,
            ready_to_go_xml=self.ready_to_go_xml,
            use_initial_tree=self.use_initial_tree,
            origin_start_addition=self.origin_start_addition,
            initial_tree_type=self.initial_tree_type,
            origin_upper_addition=self.origin_upper_addition
        )
        self.record_parameters()
        self.sampling_prop_partition_dates = None
        self.rt_partition_dates = None

    def retrieve_phase_3_params(self):
        """
        Retrieve parameters used in phase 3 of the workflow.

        Returns
        -------
        Dictionary of parameter names and values.
        """
        phase_3_params = self.retrieve_params(['save_dir',
                                               'template_xml_path',
                                               'use_initial_tree',
                                               'rt_dims',
                                               'sampling_prop_dims',
                                               'zero_sampling_before_first_sample',
                                               'collection_date_field',
                                               'sample_id_field',
                                               'origin_start_addition',
                                               'origin_upper_addition',
                                               'log_file_basename',
                                               'origin_prior',
                                               'chain_length',
                                               'trace_log_every',
                                               'tree_log_every',
                                               'screen_log_every',
                                               'store_state_every'])
        if isinstance(self.rt_partitions, dict) or isinstance(self.sampling_prop_partitions, dict):
            oldest_tip_date, youngest_tip_date = self.tip_date_range()
            rt_partition_dates, sampling_prop_partition_dates = _phase_3_partition_params(
                self.rt_partitions,
                self.sampling_prop_partitions,
                self.zero_sampling_before_first_sample,
                youngest_tip_date,
                oldest_tip_date)
            phase_3_params['rt_partitions'] = rt_partition_dates
            phase_3_params['sampling_prop_partitions'] = sampling_prop_partition_dates
        else:
            phase_3_params['rt_partitions'] = self.rt_partitions
            phase_3_params['sampling_prop_partitions'] = self.sampling_prop_partitions
            if self.zero_sampling_before_first_sample:
                oldest_tip_date, youngest_tip_date = self.tip_date_range()
                if self.sampling_prop_partitions is not None:
                    partition_dates = pd.to_datetime(self.sampling_prop_partitions)
                    if oldest_tip_date > min(partition_dates):
                        raise ValueError('If using zero_sampling_before_first_sample oldest partition date should be before oldest date in the list from of sampling_prop_partitions.')
                    phase_3_params['sampling_prop_partitions'] = [oldest_tip_date.strftime('%Y-%m-%d')] + phase_3_params['sampling_prop_partitions']
                else:
                    phase_3_params['sampling_prop_partitions'] = [oldest_tip_date.strftime('%Y-%m-%d')]

        with open(f'{self.save_dir}/pipeline_run_info.yml', 'r') as file:
            pipeline_run_info = yaml.safe_load(file)
        self.rt_partition_dates = phase_3_params['rt_partitions']
        pipeline_run_info['rt_partitions'] = phase_3_params['rt_partitions']
        self.sampling_prop_partition_dates = phase_3_params['sampling_prop_partitions']
        pipeline_run_info['sampling_prop_partitions'] = phase_3_params['sampling_prop_partitions']
        with open(f'{self.save_dir}/pipeline_run_info.yml', 'w') as fp:
            yaml.dump(pipeline_run_info, fp, sort_keys=True)
        fp.close()
        return phase_3_params


class ComparativeWorkflowParams(WorkflowParams):
    """Setup and check parameters for comparative workflows. """

    def __init__(self, **kwargs):
        """
        Initialise the workflow parameters (perform checks and set up).

        Parameters
        ----------
        kwargs:
            Parameters to use in the workflow.
        """
        super().__init__(**kwargs)
        #Check for not being assigned
        if self.fasta_path is None:
            raise Exception(
                "A ready_to_go_xml has not being given and" +
                "fasta_path is missing from the parameters yml file, or\n" +
                "fasta_path has been given an 'null' value (which are converted to None in python).\n" +
                "None/null values cannot be used for fasta_path."
            )
        if self.template_xml_path is None:
            raise Exception(
                "template_xml_path is missing from the parameters yml file, or\n" +
                "template_xml_path has been given an 'null' value (which are converted to None in python).\n" +
                "None/null values cannot be used for template_xml_path."
            )
        if self.log_file_basename is None:
            raise Exception(
                "log_file_basename is missing from the parameters yml file, or\n" +
                "log_file_basename has been given an 'null' value (which are converted to None in python).\n" +
                "None/null values cannot be used for log_file_basename."
            )
        if self.metadata_path is None:
            raise Exception(
                "metadata_path is missing from the parameters yml file, or\n" +
                "metadata_path has been given an 'null' value (which are converted to None in python).\n" +
                "None/null values cannot be used for metadata_path."
            )
        if self.xml_set_definitions is None:
            raise Exception(
                "xml_set_definitions is missing from the parameters yml file, or\n" +
                "xml_set_definitions has been given an 'null' value (which are converted to None in python).\n" +
                "None/null values cannot be used for xml_set_definitions."
            )
        if not self.use_initial_tree:  # None of the parameters below are needed if you are not using an inial tree.
            if self.initial_tree_type not in ['Temporal', 'Distance']:
                raise ValueError('initial_tree_type must be either "Temporal" or "Distance".')

            if self.down_sample_to is not None and (self.initial_tree_type != 'Temporal'):
                raise ValueError("Currently beast_pype's down_sampling method is tied to its Tree Time tree building." +
                                 "Therefore, to use this down sampling method an use_initial_tree = True and initial_tree_type = 'Temporal'.")

        # If seeds not given generate them.
        # BEASTs random number seed can select the same seed for multiple runs if they are launched close together in time (such as programmatically). Therefore, numpy is used to generate seeds for running BEAST.
        if self.seeds is None:
            number_of_seeds = self.number_of_beast_runs * len(self.xml_set_definitions)
            self.seeds = randint(low=1, high=int(1e6), size=number_of_seeds).tolist()

    def gen_xml_set_directories(self):
        """
        Generate the XML set directories for use in the workflow.
        """
        xml_set_directories = {}
        for xml_set in self.xml_set_definitions.keys():
            xml_set_directory = f'{os.getcwd()}/{self.save_dir}/{xml_set}'
            os.makedirs(xml_set_directory)
            xml_set_directories[xml_set] = xml_set_directory

        return xml_set_directories


    def retrieve_phase_1_params(self):
        """
        Retrieve parameters used in phase 1 of the workflow.

        Returns
        -------
        Dictionary of parameter names and values.
        """
        return self.retrieve_params(['save_dir',
                                     'xml_set_definitions',
                                     'xml_set_directories',
                                     'metadata_path',
                                     'data_filter',
                                     'sample_id_field',
                                     'collection_date_field',
                                     'fasta_path'])

    def retrieve_phase_2i_params(self):
        """
        Retrieve parameters used in phase 2i of the workflow.

        Returns
        -------
        Dictionary of parameter names and values.
        """
        parameters = self.retrieve_params(['save_dir', 'max_threads'])
        parameters['fasta_path'] = 'sequences.fasta'
        return parameters

    def retrieve_phase_2ii_params(self, xml_set_directory):
        """
        Retrieve parameters used in phase 2ii of the workflow.

        Parameters
        -------------
        xml_set_directory: str
            Path to the XML set directory to include under 'save_dir' of the returned
            dictionary.

        Returns
        -------
        Dictionary of parameter names and values.
        """
        parameters = self.retrieve_params([
            'sample_id_field',
            'collection_date_field',
            'down_sample_to',
            'root_strain_names',
            'remove_root'])
        parameters['save_dir'] = xml_set_directory
        return parameters

    def retrieve_phase_5_params(self):
        """
        Retrieve parameters used in phase 5 of the workflow.

        Returns
        -------
        Dictionary of parameter names and values.
        """
        if self.xml_set_label is not None:
            phase_5_params_dict = self.retrieve_params([
                'xml_set_label', 'kernel_name'])
        else:
            phase_5_params_dict = {'xml_set_label': 'xml set',
                                   **self.retrieve_params(['kernel_name'])}
        if self.report_template is None:
            phase_5_params_dict['report_template'] = self.default_report_template
        else:
            phase_5_params_dict['report_template'] = self.report_template
        return phase_5_params_dict

    def tip_date_range(self, xml_set_directory):
        """
        Generate tip date range.

        Parameters
        -------------
        xml_set_directory: str
            Path to XML set directory where metadata is held.

        Returns
        -------
        oldest_tip_date: datetime
            Oldest tip date.
        youngest_tip_date: datetime
            Youngest tip date.
        """
        if self.down_sample_to is None:
            metadata_path = f'{xml_set_directory}/metadata.csv'
        else:
            metadata_path = f'{xml_set_directory}/down_sampled_metadata.csv'
        if metadata_path.endswith('.tsv'):
            sep = '\t'
        elif metadata_path.endswith('.csv'):
            sep = ','
        else:
            raise ValueError('Only .tsv and .csv are supported for metadata. Use correct file extension.')
        metadata = pd.read_csv(metadata_path, parse_dates=[self.collection_date_field], sep=sep)
        youngest_tip_date = metadata[self.collection_date_field].max()
        oldest_tip_date = metadata[self.collection_date_field].min()
        return oldest_tip_date, youngest_tip_date


class GenericComparativeWorkflowParams(ComparativeWorkflowParams):
    """Set up and check Generic Comparative workflow parameters."""

    workflow_name = 'Generic Comparative'
    associated_workflow = f"{workflows_path}/Generic-Comparative.ipynb"
    default_report_template = "Generic-Comparative"

    def __init__(self, **kwargs):
        """
        Initialise the workflow parameters (perform checks and set up).

        Parameters
        ----------
        kwargs:
            Parameters to use in the workflow.
        """
        super().__init__(**kwargs)
        xml_set_directories = self.gen_xml_set_directories()
        self.record_parameters(extras_dict={'xml set directories': xml_set_directories})
        self.sampling_prop_partition_dates = None
        self.rt_partition_dates = None
        self.xml_set_directories = xml_set_directories

class BDSKYSerialComparativeWorkflowParams(ComparativeWorkflowParams):
    """Set up and check BDSKY-Serial Comparative workflow parameters."""

    workflow_name = 'BDSKY-Serial Comparative'
    associated_workflow = f"{workflows_path}/BDSKY-Serial-Comparative.ipynb"
    default_report_template = "BDSKY-Serial-Comparative"

    def __init__(self, **kwargs):
        """
        Initialise the workflow parameters (perform checks and set up).

        Parameters
        ----------
        kwargs:
            Parameters to use in the workflow.
        """
        super().__init__(**kwargs)
        _bdksy_serial_errors(
            rt_dims=self.rt_dims,
            rt_partitions=self.rt_partitions,
            sampling_prop_dims=self.sampling_prop_dims,
            sampling_prop_partitions=self.sampling_prop_partitions,
            zero_sampling_before_first_sample=self.zero_sampling_before_first_sample,
            ready_to_go_xml=None,
            use_initial_tree=self.use_initial_tree,
            origin_start_addition=self.origin_start_addition,
            initial_tree_type=self.initial_tree_type,
            origin_upper_addition=self.origin_upper_addition
        )
        xml_set_directories = self.gen_xml_set_directories()
        self.record_parameters(extras_dict={'xml set directories': xml_set_directories})
        self.xml_set_directories = xml_set_directories
        self.sampling_prop_partition_dates = None
        self.rt_partition_dates = None

    def retrieve_phase_3_params(self, xml_set, xml_set_directory):
        """
        Retrieve parameters used in phase 3 of the workflow.

        Parameters
        ----------
        xml_set: str
            Which xml_set are phase 3 parameters being retrieved for.
        xml_set_directory: str
            Path to XML set directory.

        Returns
        -------
        Dictionary of parameter names and values.
        """
        if self.rt_partitions is not None and xml_set in self.rt_partitions:
            rt_partitions_to_use = self.rt_partitions[xml_set]
        else:
            rt_partitions_to_use = self.rt_partitions

        if self.sampling_prop_partitions is not None and xml_set in self.sampling_prop_partitions:
            sampling_prop_partitions_to_use = self.sampling_prop_partitions[xml_set]
        else:
            sampling_prop_partitions_to_use = self.sampling_prop_partitions
        phase_3_xml_set_params = self.retrieve_params([
            'template_xml_path',
            'use_initial_tree',
            'rt_dims',
            'sampling_prop_dims',
            'collection_date_field',
            'sample_id_field',
            'origin_start_addition',
            'origin_upper_addition',
            'log_file_basename',
            'origin_prior',
            'chain_length',
            'trace_log_every',
            'tree_log_every',
            'screen_log_every',
            'store_state_every'])
        phase_3_xml_set_params['save_dir'] = xml_set_directory
        if isinstance(rt_partitions_to_use, dict) or isinstance(sampling_prop_partitions_to_use, dict):
            oldest_tip_date, youngest_tip_date = self.tip_date_range(xml_set_directory)
            rt_partition_dates, sampling_prop_partition_dates = _phase_3_partition_params(
                rt_partitions_to_use,
                sampling_prop_partitions_to_use,
                self.zero_sampling_before_first_sample,
                youngest_tip_date,
                oldest_tip_date)
            phase_3_xml_set_params['rt_partitions'] = rt_partition_dates
            phase_3_xml_set_params['sampling_prop_partitions'] = sampling_prop_partition_dates
        else:
            phase_3_xml_set_params['rt_partitions'] = rt_partitions_to_use
            phase_3_xml_set_params['sampling_prop_partitions'] = sampling_prop_partitions_to_use
            if self.zero_sampling_before_first_sample:
                oldest_tip_date, youngest_tip_date = self.tip_date_range(xml_set_directory)
                if sampling_prop_partitions_to_use is not None:
                    partition_dates = pd.to_datetime(sampling_prop_partitions_to_use)
                    if oldest_tip_date > min(partition_dates):
                        raise ValueError('If using zero_sampling_before_first_sample oldest partition date should be before oldest date in the list from of sampling_prop_partitions.')
                    phase_3_xml_set_params['sampling_prop_partitions'] = [oldest_tip_date.strftime('%Y-%m-%d')] + phase_3_xml_set_params['sampling_prop_partitions']
                else:
                    phase_3_xml_set_params['sampling_prop_partitions'] = [oldest_tip_date.strftime('%Y-%m-%d')]

        with open(f'{self.save_dir}/pipeline_run_info.yml', 'r') as file:
            pipeline_run_info = yaml.safe_load(file)

        if 'rt_partitions' not in pipeline_run_info:
            self.rt_partition_dates = {xml_set: phase_3_xml_set_params['rt_partitions']}
            pipeline_run_info['rt_partitions'] = {xml_set: phase_3_xml_set_params['rt_partitions']}
        else:
            self.rt_partition_dates[xml_set] = phase_3_xml_set_params['rt_partitions']
            pipeline_run_info['rt_partitions'][xml_set] = phase_3_xml_set_params['rt_partitions']

        if 'sampling_prop_partitions' not in pipeline_run_info:
            self.sampling_prop_partition_dates = {
                xml_set: phase_3_xml_set_params['sampling_prop_partitions']}
            pipeline_run_info['sampling_prop_partitions'] = {
                xml_set: phase_3_xml_set_params['sampling_prop_partitions']}
        else:
            self.sampling_prop_partition_dates[xml_set] = phase_3_xml_set_params[
                'sampling_prop_partitions']
            pipeline_run_info['sampling_prop_partitions'][xml_set] = phase_3_xml_set_params[
                'sampling_prop_partitions']

        with open(f'{self.save_dir}/pipeline_run_info.yml', 'w') as fp:
            yaml.dump(pipeline_run_info, fp, sort_keys=True)
        fp.close()
        return phase_3_xml_set_params


def check_file_for_phrase(file_path,
                          phrase='slurm_job_complete',
                          wait_time=60):
    """
    If file exists check file for presence of phrase. Not found wait_time in seconds.

    Used to check a slurm .out file if the slurm job is complete,
    provided slurm sbatch ends with `; echo slurm_job_complete`.

    Parameters
    ----------
    file_path: str
       Path to text file (slurm .out file).
    phrase: str (default='slurm_job_complete')
       Phrase to check for file present at file_path.
    wait_time : int (default=60)
       Time in seconds to wait between each check.

    Returns
    -------
    None

    """
    # Section is silenced as the slurm job is not started immediately.
    #if not os.path.exists(file_path):
    #    raise FileNotFoundError(f'file_path "{file_path}" does not exist.')
    #if not os.path.isfile(file_path):
    #    raise FileNotFoundError(f'file_path "{file_path}" is not a file.')
        
    complete_phrase_found = False
    while not complete_phrase_found:
        if os.path.isfile(file_path):
            with open(file_path) as file:
                s = file.read()
                complete_phrase_found = phrase in s
            file.close()
        if not complete_phrase_found:
            time.sleep(wait_time)

def _to_arg_string(options_without_a_value, options_needing_a_value):
    if options_without_a_value is None:
        list_all_args = []
    else:
        list_all_args = deepcopy(options_without_a_value)
    if options_needing_a_value is None:
        options_needing_a_value = {}
    list_all_args += [key + ' ' + str(value) for key, value in options_needing_a_value.items()]
    return ' '.join(list_all_args)

def setup_optimising_config(name,
                            configuration,
                            save_dir,
                            ready_to_go_xml,
                            threads_arg=None,
                            instances_arg=None,
                            cpu_arg=None
                            ):
    """
    Set up a configuration to test in Optimising-BEAST-Runs workflow.

    Parameters
    ----------
    name: str
    configuration: dict
    save_dir: str
    ready_to_go_xml: str
        Path of ready to go xml file.
    threads_arg: int
    instances_arg: int
    cpu_arg: int

    Returns
    -------
    Dictionary of parameter names and values.
    """
    configuration = deepcopy(configuration)
    save_name = deepcopy(name)
    if threads_arg is not None:
        if 'beast_options_needing_a_value' not in configuration:
            configuration['beast_options_needing_a_value'] = {}
        if '-threads' in configuration['beast_options_needing_a_value']:
            raise ValueError('thread_options should not be being if -threads inside is in beast_options_needing_a_value.')
        configuration['beast_options_needing_a_value']['-threads']=threads_arg
        save_name = f'{save_name}_threads_{str(threads_arg)}'
    if instances_arg is not None:
        if 'beast_options_needing_a_value' not in configuration:
            configuration['beast_options_needing_a_value'] = {}
        if '-instances' in configuration['beast_options_needing_a_value']:
            raise ValueError('thread_options should not be being if -instances inside is in beast_options_needing_a_value.')
        configuration['beast_options_needing_a_value']['-instances']=instances_arg
        save_name = f'{save_name}_instances_{str(instances_arg)}'
    if cpu_arg is not None:
        if 'sbatch_options_needing_a_value' not in configuration:
            configuration['sbatch_options_needing_a_value'] = {}
        if '--cpus-per-task' in configuration['sbatch_options_needing_a_value']:
            raise ValueError('cpu_options should not be being if --cpus-per-task inside is in sbatch_options_needing_a_value.')
        configuration['sbatch_options_needing_a_value']['--cpus-per-task']=cpu_arg
        save_name = f'{save_name}_cpu_{str(cpu_arg)}'
    if 'sbatch_options_needing_a_value' in configuration and '--job-name' not in configuration['sbatch_options_needing_a_value']:
        configuration['sbatch_options_needing_a_value']['--job-name'] = save_name
    if not ('sbatch_options_needing_a_value' in configuration or 'sbatch_single_wor_args'  in configuration):
        if 'max_threads' not in configuration:
            configuration['max_threads'] = multiprocessing.cpu_count() - 1


    if 'seeds' not in configuration:
        configuration['seeds'] = randint(low=1, high=int(1e6), size=configuration['number_of_beast_runs']).tolist()
    del configuration['number_of_beast_runs']
    save_dir = f'{save_dir}/{save_name}'
    os.makedirs(save_dir)
    parameters = _gen_phase_4_params(
        save_dir=save_dir,
        **configuration
    )

    parameters['beast_xml'] = ready_to_go_xml
    with open(f'{save_dir}/configuration_run_info.yml', 'w') as fp:
        yaml.dump(parameters, fp, sort_keys=True)
    fp.close()
    return parameters