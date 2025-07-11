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



def timescale(ftree, falignment, fdates, reroot='least-squares', clock_rate=None, clock_std=None,
              clock_filter=None, remove_root=True, coalescent_tc="opt",
              node_confidence_dir=None, sample_id_field = None,
              collection_date_field = None,
               **kwargs):
    """
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
    clock_filter: float or None
        If given clock filter is applies (treetime.TreeTime.clock_filter).
        This value is then used in
         n_iqd:  float
            The number of iqd intervals. The outlier nodes are those which do not fall
             into :math:`IQD\cdot n_iqd` interval (:math:`IQD` is the interval between
            75\ :sup:`th` and 25\ :sup:`th` percentiles)
    remove_root: bool, default True
    coalescent_tc: : float, str
           Value used in
           If not None, use coalescent model to correct the branch lengths by
           introducing merger costs.

           If Tc is float, it is interpreted as the coalescence time scale.

           If Tc is str, it should be one of (:code:`opt`, :code:`const`, :code:`skyline`)
    node_confidence_dir: str or None (default None)
        Path to directory to write csv containing node date confidence intervals
         (node_confidence.csv).
    sample_id_field : str, optional
        Name of column containing taxon names in fdates. If None, will use
        first column that contains 'name', 'strain', 'accession'
    collection_date_field : str, optional
        Name of column containing taxon names in fdates. If None, will use
        a column that contains the substring 'date'

    kwargs: dict, default None
        Key word arguments to pass to TreeTime.run.


    Returns
    -------
    time_tree: treetime.TreeTime
    bad_tips : list of str
    """
    bad_tips = list()
    dates = parse_dates(fdates, name_col=sample_id_field, date_col=collection_date_field)

    time_tree = TreeTime(gtr='JC69', tree=ftree, aln=falignment, dates=dates,
                         verbose=1, use_fft=True, precision='auto', rng_seed=None)

    if clock_filter is not None:
        time_tree.clock_filter(reroot=reroot, n_iqd=clock_filter, plot=False)
        tips = [x for x in time_tree.tree.get_terminals()]
        for n in tips:
            if n.bad_branch:
                time_tree.tree.prune(n)
                bad_tips.append(n.name)
        time_tree.prepare_tree()
        # remove bad tips
        if len(bad_tips):
            print("Prunning leaves :\n", "\n".join(bad_tips))
            with open('results/treetime_ignored_tips.txt', 'w') as f:
                for tip in bad_tips:
                    count = f.write("%s\n" % tip)

        time_tree.stl = {n.name: n for n in time_tree.tree.get_terminals()}
        time_tree.stl.update({n.name.upper(): n for n in time_tree.tree.get_terminals()})
        time_tree.stl.update({n.name.lower(): n for n in time_tree.tree.get_terminals()})
        time_tree.tree.root.up = None
        for node in time_tree.tree.get_nonterminals():
            for n in node.clades:
                n.up = node
    marginal = 'always'  # 'assign' # estimate confidence intervals via marginal ML and assign
    branch_length_inference = 'joint'  # auto
    resolve_polytomies = True
    max_iter = 2
    covariance = False

    time_tree.run(infer_gtr=True, root=reroot, Tc=coalescent_tc, time_marginal=marginal,
                  branch_length_mode=branch_length_inference, resolve_polytomies=resolve_polytomies,
                  max_iter=max_iter, fixed_pi=None, fixed_clock_rate=clock_rate,
                  stochastic_resolve=resolve_polytomies, vary_rate=clock_std, use_covariation=covariance,
                  raise_uncaught_exceptions=True, **kwargs)

    if remove_root and (reroot != 'least-squares' or reroot != 'best'):
        for root in reroot:
            if root not in bad_tips:
                root = root.strip()
                time_tree.tree.prune(root)

    time_tree.convert_dates()
    _write_nodes_CI(time_tree, node_confidence_dir)
    time_tree.branch_length_to_years()

    return time_tree, bad_tips


def _write_nodes_CI(time_tree, output_dir=None, fraction=0.95, to_node=True):
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
        if to_node:
            n.num_date_confidence = list(conf)
    if output_dir is not None:
        df = pd.DataFrame.from_records(records)
        df.to_csv(output_dir+'/node_confidence.csv', index=False)


def temporal_pruning_sampler(time_tree: TreeTime, sample_size: int, draws=1):
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
    for draw in range(draws):
        selection = copy.deepcopy(tip_names)
        to_prune_index = np.random.choice(n_tips, size=to_prune, p=prune_prob, replace=False)
        selection = [ele for index, ele in enumerate(selection) if index not in to_prune_index]
        selections.append(selection)

    if draws==1:
        return selections[0]
    else:
        return selections





