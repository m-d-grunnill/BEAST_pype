"""
This code has for functions timescale and write_nodes_CI been adapted from
 https://gitlab.in2p3.fr/ete/CoV-flow/-/blob/master/scripts/tree_time_scale.py?ref_type=heads.
On 2024-01-13

"""
import copy
from treetime import TreeTime
import numpy as np
from treetime.utils import parse_dates
import pandas as pd
import matplotlib.pyplot as plt
from beast_pype.fig_utils import year_decimal_to_date_tick_labels



def timescale(ftree, falignment, fdates, reroot='least-squares', clock_rate=None, clock_std=None,
              clock_filter=3.0, clock_filter_method='local', remove_root=True, coalescent_tc="opt",
              sample_id_field = None,
              collection_date_field = 'date',
              rng_seed=None,
              negative_tolerance=0.001,
               **kwargs):
    r"""
    Timescale a phylogenetic tree using tree time.

    Parameters
    ----------
    ftree: str
        Path to newick tree file.
    falignment: str
        Path to fasta_path alignment file.
    fdates: str
        Path to dates file.
    reroot: str or list of str, default `least-squares` 
        Method to reroot the tree see `treetime.TreeTime.reroot`:
        Which method should be used to find the best root. Available methods are:

            :code:`best`, `least-squares` - minimize squared residual or likelihood of
             root-to-tip regression.

            :code:`min_dev` - minimize variation of root-to-tip distance.

            :code:`oldest` - reroot on the oldest node.

            :code:`<node_name>` - reroot to the node with name :code:`<node_name>`.

            :code:`[<node_name1>, <node_name2>, ...]` - reroot to the MRCA of these
             nodes.
    clock_rate: float
        Mutation rate (substitutions per position per year).
    clock_std: float
        Standard deviation of the mutation (rates substitutions per position per year).
    clock_filter: float or None, default 3.0
        Threshold for clock filtering. The interpretation depends on `clock_filter_method`:
        - If `clock_filter_method='local'`: this is the z-score threshold. Sequences
          whose temporal signal deviates by more than this many z-scores are marked
          as outliers and pruned.
        - If `clock_filter_method='residual'`: this is n_iqd, the number of
          interquartile distance intervals for outlier detection.
        - If `clock_filter` is None, clock filtering is disabled.
        Set to None to disable clock filtering.
    clock_filter_method: str, default 'local'
        Method to use for clock filtering. Options are:
        - 'local': uses local z-score based outlier detection
          (see treetime.clock_filter_methods.local_filter).
        - 'residual': uses residual-based IQD filtering
          (see treetime.clock_filter_methods.residual_filter).
        - None: disables clock filtering.
    remove_root: bool, default True
        If True, remove the root after rerooting. This is useful if the root is not a real sample and is only used for rooting purposes.
    coalescent_tc: : float, str
           Value used in
           If not None, use coalescent model to correct the branch lengths by
           introducing merger costs.
           If Tc is float, it is interpreted as the coalescence timescale.
           If Tc is str, it should be one of (:code:`opt`, :code:`const`, :code:`skyline`)
    sample_id_field : str, optional
        Name of column containing taxon names in fdates. If None, will use
        first column that contains 'name', 'strain', 'accession'
    collection_date_field : str, default='date'
        Name of column containing taxon names in fdates. If None, will use
        a column that contains the substring 'date'.
    rng_seed: int, optional
        Seed for random number generator.
    negative_tolerance: float, default 0.001
        Tolerance for negative branch lengths. Branch lengths that are negative but have an absolute value less than this
        value will be set to 0.0.

    kwargs: dict, default None
        Key word arguments to pass to TreeTime.run.


    Returns
    -------
    time_tree: treetime.TreeTime
    bad_tips : list of str
    """
    dates = parse_dates(fdates, name_col=sample_id_field, date_col=collection_date_field)

    time_tree = TreeTime(tree=ftree, aln=falignment, dates=dates,
                         verbose=1, use_fft=True, precision='auto', rng_seed=rng_seed)

    # Let run() handle clock filtering internally via n_iqd and clock_filter_method.
    # After run(), outlier tips are marked bad_branch=True and time_tree.outliers is set.
    run_kwargs = dict(
        infer_gtr=True,
        root=reroot,
        Tc=coalescent_tc,
        time_marginal='always',
        branch_length_mode='joint',
        resolve_polytomies=True,
        max_iter=2,
        fixed_pi=None,
        fixed_clock_rate=clock_rate,
        stochastic_resolve=True,
        vary_rate=clock_std,
        use_covariation=False,
        raise_uncaught_exceptions=True,
    )
    if clock_filter is not None:
        run_kwargs['n_iqd'] = clock_filter
        run_kwargs['clock_filter_method'] = clock_filter_method
    run_kwargs.update(kwargs)

    time_tree.run(**run_kwargs)

    # Prune outlier tips marked by clock_filter
    bad_tips = [n.name for n in time_tree.tree.get_terminals() if n.bad_branch]
    if bad_tips:
        for tip_name in bad_tips:
            time_tree.tree.prune(tip_name)

    # Prune root strains if requested
    if remove_root and (reroot != 'least-squares' and reroot != 'best'):
        for root_name in reroot:
            if root_name not in bad_tips:
                time_tree.tree.prune(root_name.strip())

    time_tree.convert_dates()
    time_tree.branch_length_to_years()
    for node in time_tree.tree.find_clades():
        if hasattr(node, "branch_length") and node.branch_length is not None:
            if node.branch_length < 0 and abs(node.branch_length) < negative_tolerance:
                node.branch_length = 0.0

    return time_tree, bad_tips


def tree_nodes_ci(time_tree, fraction=0.95):
    """
    Get node confidence intervals from tree time tree.

    Parameters
    ----------
    time_tree: treetime.TreeTime
    fraction: float
        Confidence interval fraction

    Returns
    -------
    pd.Dataframe
    """
    records = []
    for n in time_tree.tree.find_clades():
        conf = time_tree.get_max_posterior_region(n, fraction=fraction)
        record = {
            'node': n.name,
            'date': n. date,
            'year_decimal': n.numdate,
            'interval_' + str(1-fraction): conf[0],
            'upper_' + str(fraction): conf[1]
        }
        records.append(record)
    return pd.DataFrame.from_records(records)

def plot_root_to_tip(time_tree, label=True, x_tick_freq='automatic'):
    """
    Plot root-to-tip regression, with outliers shown as orange dots.

    Parameters
    ----------
    time_tree: treetime.TreeTime
    label: bool, default True
        If true, label the plot.
    x_tick_freq: str, default='automatic'
        Suggested tick frequency. Options are 'automatic', 'yearly', 'quarterly',
        'monthly', 'half month' or 'weekly'.

    Returns
    -------
    fig, ax: matplotlib.figure.Figure
    """
    fig, ax = plt.subplots(1, 1)
    time_tree.plot_root_to_tip(ax=ax, label=label)

    # Plot outliers as orange dots if available
    if hasattr(time_tree, 'outliers') and time_tree.outliers is not None:
        outlier_df = time_tree.outliers
        if 'given_date' in outlier_df.columns and 'apparent_date' in outlier_df.columns:
            # Use given_date (x) vs dist2root proxy (apparent_date mapped to branch length)
            # The root-to-tip plot uses date vs dist2root, so we need the original residual info
            # outlier_df has 'given_date' and 'apparent_date' columns
            clock_rate = time_tree.clock_model['slope']
            intercept = time_tree.clock_model['intercept']
            outlier_dates = outlier_df['given_date'].values
            outlier_dist2root = clock_rate * outlier_df['apparent_date'].values + intercept
            ax.scatter(outlier_dates, outlier_dist2root, color='orange', zorder=5,
                       label='Clock Filter\nRemoved Outliers', edgecolors='black', linewidths=0.5, s=50)
            ax.legend()

    x_year_decimal = np.array([tip.raw_date_constraint
                               for tip in time_tree.tree.get_terminals()])
    # Include outlier dates in tick range calculation
    if hasattr(time_tree, 'outliers') and time_tree.outliers is not None:
        if 'given_date' in time_tree.outliers.columns:
            x_year_decimal = np.concatenate([
                x_year_decimal,
                time_tree.outliers['given_date'].values
            ])

    tick_year_decimals, tick_labels = year_decimal_to_date_tick_labels(x_year_decimal, tick_freq=x_tick_freq)
    ax.xaxis.set_ticks(tick_year_decimals)
    ax.set_xticklabels(tick_labels)
    ax.tick_params(axis='x', labelrotation=45)
    plt.setp(ax.get_xticklabels(), ha='right')
    return fig, ax



def temporal_pruning_sampler(time_tree: TreeTime, sample_size: int, draws=1, seed=None):
    """
    Sample a time tree object.

    Removal of tips is based on normalised residual of temporal signal.
    Parameters
    ----------
    time_tree : treetime.TreeTime
        Time tree to remove tips from.
    sample_size: int
        Sample size of desired tree.
    draws : int
        Number of draws.
    seed : int
        Random seed for sampling.

    Returns
    -------
   If draws>1:
        A list of lists of strain ids/names.
   If draws==1:
        A list of strain ids/names.
    """
    if not isinstance(draws, int):
        raise ValueError("draws must be an integer > 0.")
    tips = time_tree.tree.get_terminals()
    tip_names = [tip.name for tip in tips]
    n_tips = len(tips)
    to_prune = n_tips - sample_size
    if to_prune == 0:
        raise ValueError('time_tree provided has tips equal to sample_size')
    if to_prune < 0:
        raise ValueError('time_tree provided has less tips than sample_size')
    # get values of terminals
    x_year_decimal = np.array([tip.raw_date_constraint for tip in tips])
    root_to_tip_actual = np.array([tip.dist2root for tip in tips])
    root_to_tip_expected = time_tree.clock_model['slope'] * x_year_decimal + time_tree.clock_model['intercept']
    abs_residuals = np.absolute(root_to_tip_actual - root_to_tip_expected)
    prune_prob = abs_residuals / abs_residuals.sum()
    selections = []
    rng = np.random.default_rng(seed=seed)
    for draw in range(draws):
        selection = copy.deepcopy(tip_names)
        to_prune_index = rng.choice(n_tips, size=to_prune, p=prune_prob, replace=False)
        selection = [ele for index, ele in enumerate(selection) if index not in to_prune_index]
        selections.append(selection)

    if draws==1:
        return selections[0]
    else:
        return selections





