"""
Functions for the generation of beast xmls for specific BEAST models
"""
from dark.fasta import FastaReads
from beast2xml.beast2 import BEAST2XML
import ete3
import pandas as pd
from beast_pype.date_utilities import date_to_decimal
import warnings

def gen_xml_from_any_template(template_path,
                 sequences_path,
                 metadata_path,
                 output_path,
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
    metadata_df = pd.read_csv(metadata_path, parse_dates=['date'], sep=delimiter)
    metadata_df['year_decimal'] = metadata_df['date'].map(date_to_decimal)
    beast2xml.add_ages(metadata_df)
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
                         initial_tree_path=None,
                         origin_upper_height_addition=None,
                         origin_start_addition=None,
                         origin_prior=None,
                         rt_dimensions=None,
                         rt_change_dates=None,
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
    rt_dimensions: int, optional
        Number of Rt dimensions (time periods).
    rt_change_dates: : list, tuple, pd.Series or pd.DatetimeIndex of dates
        Rt estimation periods.
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
    metadata_df = pd.read_csv(metadata_path, parse_dates=['date'], sep=delimiter)
    metadata_df['year_decimal'] = metadata_df['date'].map(date_to_decimal)
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
            origin_prior = {
                'lower': youngest_tip - oldest_tip,
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
    beast2xml.add_ages(metadata_df)
    if origin_prior is not None:
        # Change Origin starting value, lower and upper limit on statenode
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
        if rt_dimensions is not None:
            raise AssertionError("Either rt_dimensions or rt_change_dates can be given but not both.")
        beast2xml.add_rate_change_dates(
            parameter="birthRateChangeTimes",
            dates=rt_change_dates)
    if rt_dimensions is not None:
        beast2xml.change_parameter_state_node(parameter='reproductiveNumber',
                                              dimension=rt_dimensions)
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