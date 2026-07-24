r"""
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
from Bio import Phylo, SeqIO
from beast_pype.fig_utils import year_decimal_to_date_tick_labels
import io
import os



def timescale(ftree, falignment, fdates, reroot='least-squares', clock_rate=None, clock_std=None,
              clock_filter=3.0, clock_filter_method='local', remove_root=True, coalescent_tc="opt",
              sample_id_field=None,
              collection_date_field='date',
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
    # Handle None defaults (when passed from workflow notebooks via papermill)
    if clock_filter_method is None:
        clock_filter_method = 'local'

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


def iterative_timescale(ftree, falignment, fdates,
                        reroot='least-squares',
                        clock_rate=None, clock_std=None,
                        clock_filter=3.0, clock_filter_method='local',
                        remove_future_tips=True,
                        remove_root=True, coalescent_tc="opt",
                        sample_id_field=None,
                        collection_date_field='date',
                        rng_seed=None,
                        negative_tolerance=0.001,
                        max_iterations=50,
                        **kwargs):
    r"""
    Iteratively run TreeTime clock filtering until no outliers remain.

    Mirrors the ITERATIVE_CLOCK_FILTER pattern: runs TreeTime, identifies
    outliers, prunes them from the tree/alignment/metadata, and repeats
    until convergence (no outliers in an iteration) or max_iterations reached.

    Parameters
    ----------
    ftree : str
        Path to newick tree file.
    falignment : str
        Path to fasta alignment file.
    fdates : str
        Path to dates file (CSV or TSV with sample_id_field and collection_date_field).
    reroot : str or list of str, default 'least-squares'
        Method or node(s) to reroot the tree.
    clock_rate : float, optional
        Fixed clock rate (substitutions/site/year).
    clock_std : float, optional
        Standard deviation of the clock rate. When provided, TreeTime uses a
        relaxed clock model (vary_rate parameter).
    clock_filter : float or None, default 3.0
        Threshold for clock filtering (z-score for 'local', n_iqd for 'residual').
        Set to None to disable clock filtering (runs once with no filtering).
    clock_filter_method : str, default 'local'
        Method for clock filtering: 'local' (z-score) or 'residual' (IQD).
    remove_future_tips : bool, default True
        If True, tips that TreeTime places in the future (numdate > youngest
        real sample date) are removed at each iteration. If False, only clock
        filter outliers are removed.
    remove_root : bool, default True
        If True and reroot is a list of tip names, prune those tips after the
        final iteration.
    coalescent_tc : float or str, default 'opt'
        Coalescent time constant for branch length correction.
    sample_id_field : str, optional
        Column name for taxon IDs in fdates.
    collection_date_field : str, default 'date'
        Column name for collection dates in fdates.
    rng_seed : int, optional
        Random seed for TreeTime.
    negative_tolerance : float, default 0.001
        Branch lengths negative but with abs < this value are set to 0.
    max_iterations : int, default 50
        Maximum number of clock filter iterations.
    **kwargs
        Additional keyword arguments passed to TreeTime.run().

    Returns
    -------
    time_tree : treetime.TreeTime
        The final time-scaled tree (after all outlier removal).
    all_outliers_df : pd.DataFrame
        DataFrame of all outliers removed across iterations, with columns:
        'iteration', 'name', and any additional info from time_tree.outliers.
    fasta_path : str
        Path to the final (filtered) alignment FASTA.
    metadata_path : str
        Path to the final (filtered) metadata file.
    """
    # Handle None defaults (when passed from workflow notebooks via papermill)
    if max_iterations is None:
        max_iterations = 50
    if clock_filter_method is None:
        clock_filter_method = 'local'
    if remove_future_tips is None:
        remove_future_tips = True

    # Work with copies so we don't modify original files
    current_tree = ftree
    current_fasta = falignment
    current_metadata = fdates

    all_outliers = []
    iteration = 0

    while iteration < max_iterations:
        iteration += 1
        print(f"=== Clock filter iteration {iteration} ===")

        # Run timescale (without remove_root — we handle that at the end)
        time_tree, bad_tips = timescale(
            ftree=current_tree,
            falignment=current_fasta,
            fdates=current_metadata,
            reroot=reroot,
            clock_rate=clock_rate,
            clock_std=clock_std,
            clock_filter=clock_filter,
            clock_filter_method=clock_filter_method,
            remove_root=False,  # Don't remove root until final iteration
            coalescent_tc=coalescent_tc,
            sample_id_field=sample_id_field,
            collection_date_field=collection_date_field,
            rng_seed=rng_seed,
            negative_tolerance=negative_tolerance,
            **kwargs
        )

        # Detect tips placed in the future (only if remove_future_tips is enabled)
        # Only consider tips that are actual sequences (exist in the alignment),
        # not internal nodes that became terminals through pruning.
        future_tips = []
        if remove_future_tips:
            seq_ids = {rec.id for rec in SeqIO.parse(current_fasta, 'fasta')}
            tips = time_tree.tree.get_terminals()
            root_names = reroot if isinstance(reroot, list) else []
            dated_tips = [t for t in tips if hasattr(t, 'numdate') and t.numdate is not None
                          and t.name not in root_names
                          and t.name in seq_ids]
            if dated_tips:
                raw_dates = [t.raw_date_constraint for t in dated_tips
                             if hasattr(t, 'raw_date_constraint')
                             and t.raw_date_constraint is not None]
                if raw_dates:
                    youngest_date = max(raw_dates)
                    future_tips = [t.name for t in dated_tips
                                   if t.numdate > youngest_date
                                   and t.name not in bad_tips]

        # Exclude root strains from outlier lists (they're kept for rooting)
        if isinstance(reroot, list):
            bad_tips = [t for t in bad_tips if t not in reroot]
            future_tips = [t for t in future_tips if t not in reroot]

        # Combine clock filter outliers and future-placed tips
        all_bad_this_iter = list(set(bad_tips + future_tips))

        # Check convergence: no clock filter outliers AND no future-placed tips
        if clock_filter is None and not future_tips:
            print(f"Converged after {iteration} iteration(s). No outliers or future-placed tips.")
            break

        if not all_bad_this_iter:
            print(f"Converged after {iteration} iteration(s). No outliers or future-placed tips.")
            break

        print(f"  Clock filter outliers: {len(bad_tips)}")
        if remove_future_tips:
            print(f"  Tips placed in future: {len(future_tips)}")

        # Collect outlier info for this iteration
        outlier_info = getattr(time_tree, 'outliers', None)
        if outlier_info is not None and isinstance(outlier_info, pd.DataFrame):
            iter_df = outlier_info.copy()
            if iter_df.index.name != 'name' and 'name' not in iter_df.columns:
                iter_df.index.name = 'name'
                iter_df = iter_df.reset_index()
            # Only keep rows for tips actually in bad_tips
            if 'name' in iter_df.columns:
                iter_df = iter_df[iter_df['name'].isin(bad_tips)]
        else:
            iter_df = pd.DataFrame({'name': bad_tips, 'iteration': iteration,
                                    'diagnosis': 'clock_filter'})

        # Add future-placed tips to the outlier DataFrame
        if future_tips:
            future_df = pd.DataFrame({
                'name': future_tips,
                'iteration': iteration,
                'diagnosis': 'placed_in_future'
            })
            iter_df = pd.concat([iter_df, future_df], ignore_index=True)

        iter_df['iteration'] = iteration
        all_outliers.append(iter_df)

        # Use all_bad_this_iter for filtering (both clock outliers and future-placed)
        bad_tips = all_bad_this_iter

        # Filter FASTA — remove outliers
        filtered_seqs = [rec for rec in SeqIO.parse(current_fasta, 'fasta')
                         if rec.id not in bad_tips]
        filtered_fasta = current_fasta.replace('.fasta', f'_iter{iteration}.fasta')
        if filtered_fasta == current_fasta:
            filtered_fasta = current_fasta + f'.iter{iteration}'
        with open(filtered_fasta, 'w') as handle:
            SeqIO.write(filtered_seqs, handle, 'fasta')
        current_fasta = filtered_fasta

        # Filter metadata — remove outliers
        if current_metadata.endswith('.tsv'):
            meta_df = pd.read_csv(current_metadata, sep='\t')
        else:
            meta_df = pd.read_csv(current_metadata)

        id_col = sample_id_field if sample_id_field and sample_id_field in meta_df.columns else meta_df.columns[0]
        meta_df = meta_df[~meta_df[id_col].isin(bad_tips)]
        filtered_meta = current_metadata.replace('.csv', f'_iter{iteration}.csv').replace('.tsv', f'_iter{iteration}.tsv')
        if filtered_meta == current_metadata:
            filtered_meta = current_metadata + f'.iter{iteration}'
        meta_df.to_csv(filtered_meta, index=False)
        current_metadata = filtered_meta

        # Use the pruned timetree as input for next iteration
        pruned_tree_path = current_tree.replace('.nwk', f'_iter{iteration}.nwk')
        if pruned_tree_path == current_tree:
            pruned_tree_path = current_tree + f'.iter{iteration}'
        Phylo.write(time_tree.tree, pruned_tree_path, format='newick',
                    format_branch_length='%1.8f')
        current_tree = pruned_tree_path

    else:
        print(f"WARNING: Reached max_iterations ({max_iterations}) without convergence.")

    # Now prune root strains if requested (only on the final tree)
    if remove_root and isinstance(reroot, list):
        for root_name in reroot:
            targets = [t for t in time_tree.tree.get_terminals()
                       if t.name == root_name.strip()]
            for t in targets:
                time_tree.tree.prune(t)

    # Combine all outlier DataFrames
    if all_outliers:
        all_outliers_df = pd.concat(all_outliers, ignore_index=True)
        # Reorder columns so 'iteration' and 'name' are first
        cols = ['iteration', 'name'] + [c for c in all_outliers_df.columns
                                         if c not in ('iteration', 'name')]
        all_outliers_df = all_outliers_df[cols]
    else:
        all_outliers_df = pd.DataFrame(columns=['iteration', 'name'])

    return time_tree, all_outliers_df, current_fasta, current_metadata


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

def plot_root_to_tip(time_tree, outliers_df=None, x_tick_freq='automatic'):
    """
    Plot root-to-tip regression with clock model line and outliers as orange dots.

    Parameters
    ----------
    time_tree : treetime.TreeTime
        A fitted TreeTime object with clock_model attribute.
    outliers_df : pd.DataFrame, optional
        DataFrame of outliers with 'numdate' and 'dist2root' columns.
        If None, checks time_tree.outliers for single-iteration outlier data.
    x_tick_freq : str, default='automatic'
        Suggested tick frequency. Options are 'automatic', 'yearly', 'quarterly',
        'monthly', 'half month' or 'weekly'.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    # Extract tip data from time_tree
    tips = time_tree.tree.get_terminals()
    xi = np.array([tip.raw_date_constraint for tip in tips])
    yi = np.array([tip.dist2root for tip in tips])

    # Clock model parameters
    rate = time_tree.clock_model.get('slope', None)
    intercept = time_tree.clock_model.get('intercept', None)
    r_val = time_tree.clock_model.get('r_val', None)
    r_sq = r_val**2 if r_val is not None else None
    t_mrca = -intercept / rate if (rate and intercept is not None) else None

    # Fall back to time_tree.outliers if no outliers_df provided
    if outliers_df is None and hasattr(time_tree, 'outliers') and time_tree.outliers is not None:
        outliers_df = time_tree.outliers.copy()
        # Standardise column names if needed
        if 'given_date' in outliers_df.columns and 'numdate' not in outliers_df.columns:
            outliers_df['numdate'] = outliers_df['given_date']
        if 'dist2root' not in outliers_df.columns and 'apparent_date' in outliers_df.columns:
            # Estimate dist2root from apparent_date using clock model
            outliers_df['dist2root'] = rate * outliers_df['apparent_date'] + intercept

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    if len(xi) > 0:
        # Retained tips (blue)
        ax.scatter(xi, yi, s=20, alpha=0.7, zorder=2, label='tips')

        # Outliers as orange dots
        n_outliers_total = 0
        if outliers_df is not None and 'numdate' in outliers_df.columns and 'dist2root' in outliers_df.columns:
            out_x = pd.to_numeric(outliers_df['numdate'], errors='coerce')
            out_y = pd.to_numeric(outliers_df['dist2root'], errors='coerce')
            mask = out_x.notna() & out_y.notna()
            n_outliers_total = int(mask.sum())
            if mask.any():
                ax.scatter(out_x[mask], out_y[mask], s=30, c='orange',
                           edgecolors='black', linewidths=0.5,
                           zorder=3, label='outliers removed')

        # Regression line from TMRCA to youngest tip
        if rate and intercept is not None:
            youngest = xi.max()
            if outliers_df is not None and 'numdate' in outliers_df.columns:
                out_dates = pd.to_numeric(outliers_df['numdate'], errors='coerce').dropna()
                if len(out_dates) > 0:
                    youngest = max(youngest, out_dates.max())
            time_span = youngest - t_mrca
            x_lo = t_mrca - 0.02 * time_span
            x_hi = youngest + 0.02 * time_span
            x_line = np.array([x_lo, x_hi])
            y_line = rate * x_line + intercept

            label_str = f'rate={rate:.2e} subs/site/yr\nroot date: {t_mrca:.1f}'
            if r_sq:
                label_str += f'\nR\u00b2={r_sq:.4f}'
            ax.plot(x_line, y_line, c='k', lw=2, label=label_str)
            ax.set_xlim([x_lo, x_hi])

        ax.set_ylabel('root-to-tip distance', fontsize=14)
        ax.set_xlabel('date', fontsize=14)
        ax.ticklabel_format(useOffset=False)
        ax.tick_params(labelsize=11)
        ax.set_ylim([0, 1.1 * np.max(yi)])

        # Count outliers outside the plot window and update legend
        n_outliers_outside = 0
        if n_outliers_total > 0:
            y_upper = ax.get_ylim()[1]
            x_lims = ax.get_xlim()
            out_x_vals = pd.to_numeric(outliers_df['numdate'], errors='coerce')
            out_y_vals = pd.to_numeric(outliers_df['dist2root'], errors='coerce')
            valid = out_x_vals.notna() & out_y_vals.notna()
            outside = valid & ((out_y_vals > y_upper) | (out_x_vals < x_lims[0]) | (out_x_vals > x_lims[1]))
            n_outliers_outside = int(outside.sum())
        handles, labels = ax.get_legend_handles_labels()
        if n_outliers_outside > 0:
            for i, lbl in enumerate(labels):
                if 'outliers removed' in lbl:
                    labels[i] = f'outliers removed ({n_outliers_outside}/{n_outliers_total} outside window)'
        ax.legend(handles, labels, fontsize=11)

        # Apply date tick labels
        all_dates = xi.copy()
        if outliers_df is not None and 'numdate' in outliers_df.columns:
            out_dates = pd.to_numeric(outliers_df['numdate'], errors='coerce').dropna().values
            if len(out_dates) > 0:
                all_dates = np.concatenate([all_dates, out_dates])
        tick_year_decimals, tick_labels = year_decimal_to_date_tick_labels(all_dates, tick_freq=x_tick_freq)
        ax.xaxis.set_ticks(tick_year_decimals)
        ax.set_xticklabels(tick_labels)
        ax.tick_params(axis='x', labelrotation=45)
        plt.setp(ax.get_xticklabels(), ha='right')
    else:
        ax.text(0.5, 0.5, 'No tip data available for root-to-tip plot',
                ha='center', va='center', transform=ax.transAxes)

    plt.tight_layout()
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





