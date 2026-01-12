"""
Functions for the generation of beast xmls for specific BEAST models
"""
from dark.fasta import FastaReads
from beast2xml.beast2 import BEAST2XML
import ete3
import pandas as pd
from beast_pype.date_utilities import date_to_decimal
import warnings

def two_df_cols_to_dict(df, key, value):
    return df[[key, value]].set_index(key).to_dict()[value]

def gen_xml_from_any_template(template_path,
                              sequences_path,
                              metadata_path,
                              output_path,
                              sample_id_field='strain',
                              collection_date_field='date',
                              initial_tree_path=None,
                              log_file_basename=None,
                              chain_length=None,
                              trace_log_every=None,
                              tree_log_every=None,
                              screen_log_every=None,
                              store_state_every=None):
    """
    Generate a BEAST 2 xml from a BEAST 2 xml template for any model.

    Parameters
    ----------
    template_path: str
        Path to template_xml_path.
    sequences_path: str
        Path to sequences must be fasta_path.
    metadata_path: str
        Path to metadata_update must be csv.
    sample_id_field: str, default 'strain'
        Field used to identify samples.
    collection_date_field: str, default 'date'
        Field to use as sequence collection date.
    output_path: str
        Path to save output xml to.
    initial_tree_path: str, optional
        Path to initial_tree    must be Newick file (nwk).
    log_file_basename: str, optional
            The base filename to write logs to. A .log or .trees suffix will be appended
            to this to make the actual log file names.  If None, the log file names in
            the template will be retained.
    chain_length : int, optional
        The length of the MCMC chain. If C{None}, the value in the template will
         be retained.
    trace_log_every: int, optional
        Specifying how often to write to the trace log file. If None, the value in the
        template will be retained.
    tree_log_every: int, optional
        Specifying how often to write to the file_path log file. If None, the value in the
        template will be retained.
    screen_log_every: int, optional
        Specifying how often to write to the terminal (screen) log. If None, the
        value in the template will be retained.
    store_state_every : int, optional
        Specifying how often to write MCMC state file. If None, the
        value in the template will be retained.

    """
    if metadata_path.endswith('.tsv'):
        delimiter = '\t'
    elif metadata_path.endswith('.csv'):
        delimiter = ','
    else:
        raise TypeError(
            f"metadata_path must be a csv or tsv file, ending with the appropriate file extension. Value given is {metadata_path}")
    beast2xml = BEAST2XML(template=template_path)
    seqs = FastaReads([sequences_path])
    beast2xml.add_sequences(seqs)
    beast2xml.add_dates(date_data=metadata_path,
                        seperator=delimiter,
                        sample_id_field=sample_id_field,
                        collection_date_field=collection_date_field)
    if initial_tree_path is not None:
        beast2xml.add_initial_tree(initial_tree_path)
    beast2xml.to_xml(
        output_path,
        chain_length=chain_length,
        log_file_basename=log_file_basename,
        trace_log_every=trace_log_every,
        tree_log_every=tree_log_every,
        screen_log_every=screen_log_every,
        store_state_every=store_state_every,
    )


def gen_bdsky_serial_xml(template_path,
                         sequences_path,
                         metadata_path,
                         output_path,
                         collection_date_field='date',
                         sample_id_field='strain',
                         initial_tree_path=None,
                         origin_upper_height_addition=None,
                         origin_start_addition=None,
                         origin_prior=None,
                         rt_dims=None,
                         rt_change_dates=None,
                         sampling_prop_dims=None,
                         sampling_prop_change_dates=None,
                         zero_sampling_before_first_sample=True,
                         log_file_basename=None,
                         chain_length=None,
                         trace_log_every=None,
                         tree_log_every=None,
                         screen_log_every=None,
                         store_state_every=None):
    """
    Generate a BDSKY xml with initial tree inserted.

    Parameters
    ----------
    template_path: str
        Path to template_xml_path.
    sequences_path: str
        Path to sequences must be fasta_path.
    metadata_path: str
        Path to metadata_update must be csv.
    output_path: str
        Path to save output xml to.
    sample_id_field: str, default 'strain'
        Field used to identify samples.
    collection_date_field: str, default 'date'
        Field to use as sequence collection date.
    initial_tree_path: str, optional
        Path to initial_tree    must be Newick file (nwk).
    origin_upper_height_addition: int or float, optional
        Value to add to tree height for upper limit of origin prior. Origin prior is
         uniformly distributed.
    origin_start_addition: int or float, optional
        Value to add to tree height for starting value. . Origin prior is uniformly
         distributed.
    origin_prior: dict {'lower': float, 'upper': float, 'start': float}, optional
        Details of the origin prior assumed to be uniformly distributed.
    rt_dims: int, optional
        Number of Rt dimensions (time periods).
    rt_change_dates: : list, tuple, pd.Series or pd.DatetimeIndex of datetimes
        Internal partitions of Rt estimation periods.
    sampling_prop_dims: int, optional
        Number of sampling proportion dimensions (time periods).
    sampling_prop_change_dates: : list, tuple, pd.Series or pd.DatetimeIndex of dates
        Internal partitions of sampling proportion estimation periods.
    zero_sampling_before_first_sample: bool, default True
        Whether to have zero sampling before date the first sample was collected.
    log_file_basename: str, optional
            The base filename to write logs to. A .log or .trees suffix will be appended
            to this to make the actual log file names.  If None, the log file names in
            the template will be retained.
    chain_length : int, optional
        The length of the MCMC chain. If C{None}, the value in the template will
         be retained.
    trace_log_every: int, optional
        Specifying how often to write to the trace log file. If None, the value in the
        template will be retained.
    tree_log_every: int, optional
        Specifying how often to write to the file_path log file. If None, the value in the
        template will be retained.
    screen_log_every: int, optional
        Specifying how often to write to the terminal (screen) log. If None, the
        value in the template will be retained.
    store_state_every : int, optional
        Specifying how often to write MCMC state file. If None, the
        value in the template will be retained.

    """
    if metadata_path.endswith('.tsv'):
        delimiter = '\t'
    elif metadata_path.endswith('.csv'):
        delimiter = ','
    else:
        raise TypeError(
            f"metadata_path must be a csv or tsv file, ending with the appropriate file extension. Value given is {metadata_path}")
    metadata_df = pd.read_csv(metadata_path, parse_dates=[collection_date_field], sep=delimiter)
    metadata_df['year_decimal'] = metadata_df[collection_date_field].map(date_to_decimal)
    if origin_prior is None:
        if origin_upper_height_addition is not None and  origin_start_addition is not None:
            if initial_tree_path is None:
                raise ValueError('If parameterising the origin prior via ' +
                                 'origin_upper_height_addition and origin_start_addition an ' +
                                 'initial tree must be provided.')
            tree = ete3.Tree(initial_tree_path, format=1)
            furthest_leaf, tree_height = tree.get_farthest_leaf()
            youngest_tip = metadata_df.year_decimal.max()
            oldest_tip = metadata_df.year_decimal.min()
            tip_distance = youngest_tip - oldest_tip
            if tip_distance > tree_height:
                raise ValueError('tree_height must be greater than distance between youngest_tip_date and oldest_tip.')
            origin_prior = {
                'lower': tip_distance,
                'upper': tree_height + origin_upper_height_addition,
                'start': tree_height + origin_start_addition}
    else:
        warnings.warn("If using your own origin prior there is a chance" +
                      " that an origin value will be less than the tree " +
                      " height when BEAST 2 is running." +
                      " If this happens BEAST 2 will crash.\n"+
                      "We recommend using one generated by supplying:\n" +
                      "* An initial temporal tree. \n" +
                      "* origin_upper_height_addition \n" +
                      "* origin_start_addition")
    beast2xml = BEAST2XML(template=template_path)
    seqs = FastaReads([sequences_path])
    beast2xml.add_sequences(seqs)
    beast2xml.add_dates(date_data=metadata_path,
                        seperator=delimiter,
                        sample_id_field=sample_id_field,
                        collection_date_field=collection_date_field)
    if origin_prior is not None:
        # Change Origin starting value, lower and upper limit on state node
        beast2xml.change_parameter_state_node("origin", value=origin_prior["start"])
        del origin_prior["start"]
        beast2xml.change_prior("origin", "uniform", **origin_prior)
    else:
        warnings.warn("If using the Origin prior in the template xml there is a chance" +
                      " that an origin value will be less than the tree " +
                      " height when BEAST 2 is running.." +
                      " If this happens BEAST 2 will crash.\n"+
                      "We recommend using one generated by supplying:\n" +
                      "* An initial temporal tree. \n" +
                      "* origin_upper_height_addition \n" +
                      "* origin_start_addition")
    if rt_change_dates is not None:
        if rt_dims is not None:
            raise AssertionError("Either rt_dims or rt_partitions can be given but not both.")
        beast2xml.add_rate_change_dates(
            parameter="birthRateChangeTimes",
            dates=rt_change_dates)
    if rt_dims is not None:
        beast2xml.change_parameter_state_node(parameter='reproductiveNumber',
                                              dimension=rt_dims)
    if sampling_prop_change_dates is not None:
        if sampling_prop_dims is not None:
            raise AssertionError("Either sampling_prop_dims or sampling_prop_partition_freq can be given but not both.")
        if zero_sampling_before_first_sample:
            beast2xml.add_rate_change_dates(
                parameter="samplingRateChangeTimes",
                dates=sampling_prop_change_dates,
                offset_earliest=1e-6)
            beast2xml.set_dimension_values_to_0(parameter="samplingProportion")
        else:
            beast2xml.add_rate_change_dates(
                parameter="samplingRateChangeTimes",
                dates=sampling_prop_change_dates)
    if sampling_prop_dims is not None:
        beast2xml.change_parameter_state_node(parameter='samplingProportion',
                                              dimension=sampling_prop_dims)
    if initial_tree_path is not None:
        beast2xml.add_initial_tree(initial_tree_path)
    beast2xml.to_xml(
        output_path,
        chain_length=chain_length,
        log_file_basename=log_file_basename,
        trace_log_every=trace_log_every,
        tree_log_every=tree_log_every,
        screen_log_every=screen_log_every,
        store_state_every=store_state_every,
    )