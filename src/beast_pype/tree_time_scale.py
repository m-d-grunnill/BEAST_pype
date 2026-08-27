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
from matplotlib.patches import Rectangle
import io
import os


def _parse_dates_robust(fdates, name_col=None, date_col='date'):
    """Parse dates from CSV/TSV, bypassing TreeTime's parse_dates quirks.

    Reads the metadata file with the correct separator based on extension,
    resolves column names flexibly, and converts date strings to numeric
    dates (year decimals) as TreeTime expects.

    Parameters
    ----------
    fdates : str
        Path to dates file (CSV or TSV).
    name_col : str, optional
        Column name for taxon IDs. If None, auto-detected.
    date_col : str, default 'date'
        Column name for collection dates.

    Returns
    -------
    dict
        Mapping of taxon name to numeric date (year decimal).
    """
    from beast_pype.date_utilities import date_to_decimal

    # Read with correct separator
    if fdates.endswith('.tsv'):
        df = pd.read_csv(fdates, sep='	')
    else:
        df = pd.read_csv(fdates)

    # Strip whitespace from column names
    df.columns = df.columns.str.strip()

    # Resolve name column
    if name_col is None:
        for col in df.columns:
            if any(x in col.lower() for x in ['name', 'strain', 'accession']):
                name_col = col
                break
        if name_col is None:
            name_col = df.columns[0]
    elif name_col not in df.columns:
        for col in df.columns:
            if col.lower() == name_col.lower():
                name_col = col
                break

    # Resolve date column (case-insensitive + partial match)
    if date_col not in df.columns:
        for col in df.columns:
            if col.lower() == date_col.lower():
                date_col = col
                break
        if date_col not in df.columns:
            for col in df.columns:
                if date_col.lower() in col.lower():
                    date_col = col
                    break

    if date_col not in df.columns:
        raise ValueError(
            f"Date column '{date_col}' not found. "
            f"Available columns: {', '.join(df.columns)}"
        )

    # Convert dates to year decimals
    dates = {}
    for _, row in df.iterrows():
        name = str(row[name_col]).strip()
        date_val = row[date_col]
        if pd.isna(date_val):
            continue
        try:
            dt = pd.to_datetime(date_val)
            dates[name] = date_to_decimal(dt)
        except Exception:
            continue

    return dates


def timescale(ftree, falignment, fdates, reroot='least-squares', clock_rate=None, clock_std=None,
              clock_filter=3.0, clock_filter_method='local', remove_root=True, coalescent_tc="opt",
              sample_id_field=None,
              collection_date_field='date',
              rng_seed=None,
              negative_tolerance=0.001,
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
    run_kwargs = {
        'infer_gtr': True,
        'root': reroot,
        'Tc': coalescent_tc,
        'time_marginal': 'always',
        'branch_length_mode': 'joint',
        'resolve_polytomies': True,
        'max_iter': 2,
        'fixed_pi': None,
        'fixed_clock_rate': clock_rate,
        'stochastic_resolve': True,
        'vary_rate': clock_std,
        'use_covariation': False,
        'raise_uncaught_exceptions': True,
    }
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
                        remove_future_tips=None,
                        remove_root=True, coalescent_tc="opt",
                        sample_id_field=None,
                        collection_date_field='date',
                        rng_seed=None,
                        negative_tolerance=0.001,
                        max_iterations=50,
                        cache_dir=None,
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
    remove_future_tips : float, int, or None, default None
        Buffer in days beyond the youngest sample date for detecting tips
        placed in the future. Tips whose projected date exceeds
        (youngest_date + buffer_in_year_decimal) are removed. Set to 0 for
        no buffer (remove any tip projected beyond the youngest date).
        Set to None to disable future-tip removal entirely.
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
    cache_dir : str, optional
        Path to a directory for storing temporary intermediate files (pruned
        trees, filtered FASTAs and metadata) created during iterations. If
        None, intermediate files are written alongside the original input
        files. The directory is created if it does not exist.
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
    # remove_future_tips: None disables, numeric value is buffer in days
    # (converted to year decimal buffer below)

    # Set up cache directory for intermediate files
    if cache_dir is not None:
        cache_dir = os.path.abspath(cache_dir)
        os.makedirs(cache_dir, exist_ok=True)

    # Convert to absolute paths to avoid issues when papermill
    # changes the working directory during notebook execution
    ftree = os.path.abspath(ftree)
    falignment = os.path.abspath(falignment)
    fdates = os.path.abspath(fdates)

    # Work with copies so we don't modify original files
    current_fasta = falignment
    current_metadata = fdates

    # Store base names for clean iteration suffixes (avoid stacking like _iter1_iter2_iter3)
    # When cache_dir is set, write intermediate files there instead of alongside originals
    if cache_dir is not None:
        fasta_base = os.path.join(cache_dir, os.path.splitext(os.path.basename(falignment))[0])
        fasta_ext = os.path.splitext(falignment)[1]
        meta_base = os.path.join(cache_dir, os.path.splitext(os.path.basename(fdates))[0])
        meta_ext = os.path.splitext(fdates)[1]
        tree_base = os.path.join(cache_dir, os.path.splitext(os.path.basename(ftree))[0])
        tree_ext = os.path.splitext(ftree)[1]
    else:
        fasta_base, fasta_ext = os.path.splitext(falignment)
        meta_base, meta_ext = os.path.splitext(fdates)
        tree_base, tree_ext = os.path.splitext(ftree)

    all_outliers = []
    intermediate_files = []  # Track intermediate files for cleanup
    iteration = 0

    # Read original tree ONCE before the loop
    original_tree = Phylo.read(ftree, 'newick')

    while iteration < max_iterations:
        iteration += 1
        print(f"=== Clock filter iteration {iteration} ===")

        # Prune a copy of the original tree to match the current alignment.
        # This avoids MissingDataError when tips in the tree have no
        # corresponding sequence in the filtered FASTA.
        seq_ids_for_tree = {rec.id for rec in SeqIO.parse(current_fasta, 'fasta')}
        tree_copy = copy.deepcopy(original_tree)
        tips_to_remove = [t.name for t in tree_copy.get_terminals()
                          if t.name not in seq_ids_for_tree]
        for tip_name in tips_to_remove:
            tree_copy.prune(tip_name)
        pruned_tree_path = f"{tree_base}_pruned{tree_ext}"
        Phylo.write(tree_copy, pruned_tree_path, format='newick',
                    format_branch_length='%1.8f')

        # Run timescale (without remove_root — we handle that at the end)
        # If clock filter rerooting fails (e.g. least-squares cannot find a
        # positive-rate root), fall back to running without clock filtering.
        try:
            time_tree, bad_tips = timescale(
                ftree=pruned_tree_path,
                falignment=current_fasta,
                fdates=current_metadata,
                reroot=reroot,
                clock_rate=clock_rate,
                clock_std=clock_std,
                clock_filter=clock_filter,
                clock_filter_method=clock_filter_method,
                remove_root=False,
                coalescent_tc=coalescent_tc,
                sample_id_field=sample_id_field,
                collection_date_field=collection_date_field,
                rng_seed=rng_seed,
                negative_tolerance=negative_tolerance,
                **kwargs
            )
        except Exception as e:
            err_str = str(e).lower()
            recoverable = clock_filter is not None and (
                "rerooting failed" in err_str
                or "infinity" in err_str
                or "overflow" in err_str
            )
            if recoverable:
                print(f"  TreeTime failed with clock filter: {e}")
                print("  Falling back to running without clock filter.")
                time_tree, bad_tips = timescale(
                    ftree=pruned_tree_path,
                    falignment=current_fasta,
                    fdates=current_metadata,
                    reroot=reroot,
                    clock_rate=clock_rate,
                    clock_std=clock_std,
                    clock_filter=None,
                    remove_root=False,
                    coalescent_tc=coalescent_tc,
                    sample_id_field=sample_id_field,
                    collection_date_field=collection_date_field,
                    rng_seed=rng_seed,
                    negative_tolerance=negative_tolerance,
                    **kwargs
                )
            else:
                raise

        # Detect tips placed in the future using clock model projection.
        # Uses (dist2root - intercept) / rate to compute projected date,
        # which is more reliable than TreeTime's inferred numdate.
        # Only considers actual sequences (in alignment), not internal nodes
        # that became terminals through pruning (NODE_ prefix).
        future_tips = []
        future_tip_data = {}
        if remove_future_tips is not None:
            # Convert buffer from days to year decimal
            future_buffer = float(remove_future_tips) / 365.25
            rate = time_tree.clock_model.get('slope', None)
            intercept_val = time_tree.clock_model.get('intercept', None)
            if rate and intercept_val is not None and rate > 0:
                # Reuse seq_ids_for_tree from tree pruning step (already parsed current_fasta)
                tips = time_tree.tree.get_terminals()
                root_names = reroot if isinstance(reroot, list) else []
                dated_tips = [t for t in tips
                              if hasattr(t, 'raw_date_constraint') and t.raw_date_constraint is not None
                              and t.name not in root_names
                              and t.name in seq_ids_for_tree
                              and not str(t.name).startswith('NODE_')]
                if dated_tips:
                    youngest_date = max(t.raw_date_constraint for t in dated_tips)
                    threshold_date = youngest_date + future_buffer
                    for t in dated_tips:
                        dist = getattr(t, 'dist2root', None)
                        if dist is not None:
                            projected_date = (dist - intercept_val) / rate
                            if projected_date > threshold_date and t.name not in bad_tips:
                                future_tips.append(t.name)
                                future_tip_data[t.name] = {
                                    'numdate': t.raw_date_constraint,
                                    'dist2root': dist,
                                    'projected_date': projected_date,
                                }

        # Exclude root strains from outlier lists (they're kept for rooting)
        if isinstance(reroot, list):
            bad_tips = [t for t in bad_tips if t not in reroot]
            future_tips = [t for t in future_tips if t not in reroot]

        # Combine clock filter outliers and future-placed tips
        all_bad_this_iter = list(set(bad_tips + future_tips))

        # Check convergence: no clock filter outliers AND no future-placed tips
        if not all_bad_this_iter:
            print(f"Converged after {iteration} iteration(s). No outliers or future-placed tips.")
            break

        print(f"  Clock filter outliers: {len(bad_tips)}")
        if remove_future_tips is not None:
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
            future_records = []
            for tip_name in future_tips:
                record = {'name': tip_name, 'iteration': iteration, 'diagnosis': 'placed_in_future'}
                if tip_name in future_tip_data:
                    record.update(future_tip_data[tip_name])
                future_records.append(record)
            future_df = pd.DataFrame(future_records)
            iter_df = pd.concat([iter_df, future_df], ignore_index=True)

        iter_df['iteration'] = iteration
        all_outliers.append(iter_df)

        # Use all_bad_this_iter for filtering (both clock outliers and future-placed)
        bad_tips = all_bad_this_iter

        # Filter FASTA — remove outliers
        filtered_seqs = [rec for rec in SeqIO.parse(current_fasta, 'fasta')
                         if rec.id not in bad_tips]
        filtered_fasta = f"{fasta_base}_iter{iteration}{fasta_ext}"
        with open(filtered_fasta, 'w') as handle:
            SeqIO.write(filtered_seqs, handle, 'fasta')
        intermediate_files.append(filtered_fasta)
        current_fasta = filtered_fasta

        # Filter metadata — remove outliers
        if current_metadata.endswith('.tsv'):
            meta_df = pd.read_csv(current_metadata, sep='\t')
        else:
            meta_df = pd.read_csv(current_metadata)

        id_col = sample_id_field if sample_id_field and sample_id_field in meta_df.columns else meta_df.columns[0]
        meta_df = meta_df[~meta_df[id_col].isin(bad_tips)]
        filtered_meta = f"{meta_base}_iter{iteration}{meta_ext}"
        # Preserve original separator: TSV files must be written with tab separator
        meta_sep = '\t' if filtered_meta.endswith('.tsv') else ','
        meta_df.to_csv(filtered_meta, index=False, sep=meta_sep)
        intermediate_files.append(filtered_meta)
        current_metadata = filtered_meta



    else:
        print(f"WARNING: Reached max_iterations ({max_iterations}) without convergence.")

    # Clean up pruned tree file
    pruned_tree_path = f"{tree_base}_pruned{tree_ext}"
    if os.path.exists(pruned_tree_path):
        os.remove(pruned_tree_path)

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

    # Clean up intermediate files (keep only the final versions)
    for f in intermediate_files:
        if f != current_fasta and f != current_metadata:
            if os.path.exists(f):
                os.remove(f)

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

def plot_root_to_tip(time_tree, outliers_df=None):
    """
    Plot root-to-tip regression with clock model line, outliers (orange), and
    future-placed tips (red).

    Parameters
    ----------
    time_tree : treetime.TreeTime
        A fitted TreeTime object with clock_model attribute.
    outliers_df : pd.DataFrame, optional
        DataFrame of outliers with columns for date (numdate/given_date) and
        optionally dist2root. Must have a 'diagnosis' column to distinguish
        clock outliers from future-placed tips.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    tips = time_tree.tree.get_terminals()
    xi = np.array([tip.raw_date_constraint for tip in tips
                   if hasattr(tip, 'raw_date_constraint') and tip.raw_date_constraint is not None])
    yi = np.array([tip.dist2root for tip in tips
                   if hasattr(tip, 'raw_date_constraint') and tip.raw_date_constraint is not None])

    rate = time_tree.clock_model.get('slope', None)
    intercept = time_tree.clock_model.get('intercept', None)
    r_val = time_tree.clock_model.get('r_val', None)
    r_sq = r_val ** 2 if r_val is not None else None
    t_mrca = -intercept / rate if (rate and intercept is not None and rate != 0) else None

    # Separate clock outliers from future tips
    clock_outliers_df = None
    future_tips_df = None
    if outliers_df is not None and 'diagnosis' in outliers_df.columns:
        future_mask = outliers_df['diagnosis'] == 'placed_in_future'
        future_tips_df = outliers_df[future_mask] if future_mask.any() else None
        clock_outliers_df = outliers_df[~future_mask] if (~future_mask).any() else None
    elif outliers_df is not None:
        clock_outliers_df = outliers_df

    # Resolve numdate column — could be 'numdate' or 'given_date'
    def get_numdate_col(df):
        if df is None:
            return None
        if 'numdate' in df.columns:
            return 'numdate'
        if 'given_date' in df.columns:
            return 'given_date'
        if 'raw_date_constraint' in df.columns:
            return 'raw_date_constraint'
        return None

    def get_dist2root_col(df):
        if df is None:
            return None
        if 'dist2root' in df.columns:
            return 'dist2root'
        return None

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    if len(xi) > 0:
        ax.scatter(xi, yi, s=20, alpha=0.7, zorder=2, label='tips')

        # Clock outliers (orange)
        n_clock_total = 0
        if clock_outliers_df is not None:
            numdate_col = get_numdate_col(clock_outliers_df)
            dist_col = get_dist2root_col(clock_outliers_df)
            if numdate_col:
                out_x = pd.to_numeric(clock_outliers_df[numdate_col], errors='coerce')
                if dist_col:
                    out_y = pd.to_numeric(clock_outliers_df[dist_col], errors='coerce')
                elif rate and intercept is not None:
                    out_y = rate * out_x + intercept
                else:
                    out_y = pd.Series([np.nan] * len(out_x))
                mask = out_x.notna() & out_y.notna()
                n_clock_total = int(mask.sum())
                if mask.any():
                    ax.scatter(out_x[mask], out_y[mask], s=30, c='orange',
                               edgecolors='black', linewidths=0.5,
                               zorder=3, label='clock outliers removed')

        # Future tips (red)
        n_future_total = 0
        if future_tips_df is not None:
            numdate_col = get_numdate_col(future_tips_df)
            dist_col = get_dist2root_col(future_tips_df)
            if numdate_col:
                fut_x = pd.to_numeric(future_tips_df[numdate_col], errors='coerce')
                if dist_col:
                    fut_y = pd.to_numeric(future_tips_df[dist_col], errors='coerce')
                elif rate and intercept is not None:
                    fut_y = rate * fut_x + intercept
                else:
                    fut_y = pd.Series([np.nan] * len(fut_x))
                mask_f = fut_x.notna() & fut_y.notna()
                n_future_total = int(mask_f.sum())
                if mask_f.any():
                    ax.scatter(fut_x[mask_f], fut_y[mask_f], s=30, c='red',
                               edgecolors='black', linewidths=0.5,
                               zorder=3, label='future tips removed')

        # Regression line
        if rate and intercept is not None and t_mrca is not None:
            youngest = xi.max()
            if n_clock_total > 0:
                out_dates = pd.to_numeric(clock_outliers_df[get_numdate_col(clock_outliers_df)], errors='coerce').dropna()
                if len(out_dates) > 0:
                    youngest = max(youngest, out_dates.max())
            if n_future_total > 0:
                fut_dates = pd.to_numeric(future_tips_df[get_numdate_col(future_tips_df)], errors='coerce').dropna()
                if len(fut_dates) > 0:
                    youngest = max(youngest, fut_dates.max())
            time_span = youngest - t_mrca
            x_lo = t_mrca - 0.2 * time_span
            x_hi = youngest + 0.2 * time_span
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

        # Update legend with outside-window counts
        handles, labels = ax.get_legend_handles_labels()
        y_upper = ax.get_ylim()[1]
        x_lims = ax.get_xlim()

        if n_clock_total > 0:
            numdate_col = get_numdate_col(clock_outliers_df)
            dist_col = get_dist2root_col(clock_outliers_df)
            out_x_vals = pd.to_numeric(clock_outliers_df[numdate_col], errors='coerce')
            if dist_col:
                out_y_vals = pd.to_numeric(clock_outliers_df[dist_col], errors='coerce')
            elif rate and intercept is not None:
                out_y_vals = rate * out_x_vals + intercept
            else:
                out_y_vals = pd.Series([np.nan] * len(out_x_vals))
            valid = out_x_vals.notna() & out_y_vals.notna()
            outside = valid & ((out_y_vals > y_upper) | (out_x_vals < x_lims[0]) | (out_x_vals > x_lims[1]))
            n_outside = int(outside.sum())
            if n_outside > 0:
                for i, lbl in enumerate(labels):
                    if 'clock outliers removed' in lbl:
                        labels[i] = f'clock outliers removed ({n_outside}/{n_clock_total} outside window)'

        if n_future_total > 0:
            numdate_col = get_numdate_col(future_tips_df)
            dist_col = get_dist2root_col(future_tips_df)
            fut_x_vals = pd.to_numeric(future_tips_df[numdate_col], errors='coerce')
            if dist_col:
                fut_y_vals = pd.to_numeric(future_tips_df[dist_col], errors='coerce')
            elif rate and intercept is not None:
                fut_y_vals = rate * fut_x_vals + intercept
            else:
                fut_y_vals = pd.Series([np.nan] * len(fut_x_vals))
            valid_f = fut_x_vals.notna() & fut_y_vals.notna()
            outside_f = valid_f & ((fut_y_vals > y_upper) | (fut_x_vals < x_lims[0]) | (fut_x_vals > x_lims[1]))
            n_f_outside = int(outside_f.sum())
            if n_f_outside > 0:
                for i, lbl in enumerate(labels):
                    if 'future tips removed' in lbl:
                        labels[i] = f'future tips removed ({n_f_outside}/{n_future_total} outside window)'

        ax.legend(handles, labels, fontsize=11)

        # Date tick labels — prefer axis range, fall back to data points if overflow
        x_lims = ax.get_xlim()
        try:
            tick_year_decimals, tick_labels = year_decimal_to_date_tick_labels(
                np.array([x_lims[0], x_lims[1]])
            )
        except (OverflowError, ValueError):
            # Axis limits produce invalid dates — use actual data range instead
            all_dates = xi.copy()
            if clock_outliers_df is not None:
                numdate_col = get_numdate_col(clock_outliers_df)
                if numdate_col:
                    out_dates = pd.to_numeric(clock_outliers_df[numdate_col], errors='coerce').dropna().values
                    if len(out_dates) > 0:
                        all_dates = np.concatenate([all_dates, out_dates])
            if future_tips_df is not None:
                numdate_col = get_numdate_col(future_tips_df)
                if numdate_col:
                    fut_dates = pd.to_numeric(future_tips_df[numdate_col], errors='coerce').dropna().values
                    if len(fut_dates) > 0:
                        all_dates = np.concatenate([all_dates, fut_dates])
            tick_year_decimals, tick_labels = year_decimal_to_date_tick_labels(all_dates)

        ax.xaxis.set_ticks(tick_year_decimals)
        ax.set_xticklabels(tick_labels)
        ax.tick_params(axis='x', labelrotation=45)
        plt.setp(ax.get_xticklabels(), ha='right')

    plt.tight_layout()
    return fig, ax


def plot_temporal_tree(time_tree):
    """
    Plot temporal tree (branch lengths already in years from branch_length_to_years).

    Parameters
    ----------
    time_tree : treetime.TreeTime
        A fitted TreeTime object after branch_length_to_years() has been called.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    tips = time_tree.tree.get_terminals()
    nleafs = len(tips)

    # Root date
    root = time_tree.tree.root
    root_date = root.numdate if hasattr(root, 'numdate') and root.numdate is not None else None

    # Set root branch length to 0 for plotting
    root.branch_length = 0.0

    fig_height = max(8, nleafs * 0.2)
    fig, ax = plt.subplots(1, 1, figsize=(12, fig_height))

    label_func = lambda x: x.name if (x.is_terminal() and nleafs < 30) else ''
    Phylo.draw(time_tree.tree, axes=ax, do_show=False, label_func=label_func)

    if root_date is not None:
        # x-axis: Phylo.draw places root at x=0, tips at x=sum(branch_lengths)
        # After branch_length_to_years, x positions are years from root
        # So date = root_date + x_position
        xlim = ax.get_xlim()
        # Generate ticks covering the full x-axis range (converted to dates)
        date_lo = root_date + xlim[0]
        date_hi = root_date + xlim[1]
        tick_year_decimals, tick_labels = year_decimal_to_date_tick_labels(
            np.array([date_lo, date_hi])
        )
        # Convert date ticks to x-position (x = date - root_date)
        xticks = tick_year_decimals - root_date
        ax.set_xticks(xticks)
        ax.set_xticklabels(tick_labels, rotation=45, ha='right')
        ax.set_xlim([0, xlim[1]])

        # Add alternating shaded bands
        tip_dates = [t.numdate for t in tips if hasattr(t, 'numdate') and t.numdate is not None]
        ylim = ax.get_ylim()
        xlim_new = ax.get_xlim()
        if tip_dates:
            date_range = max(tip_dates) - root_date
            if date_range > 5:
                step = 1.0
            elif date_range > 2:
                step = 0.25
            elif date_range > 1:
                step = 1.0 / 12
            else:
                step = 1.0 / 52

            for idx, year_offset in enumerate(np.arange(0, xlim_new[1] + step, step)):
                shade = 0.92 if idx % 2 == 0 else 0.96
                r = Rectangle((year_offset, ylim[0]), step, ylim[1] - ylim[0],
                               facecolor=[shade] * 3, edgecolor='none', zorder=0)
                ax.add_patch(r)

    ax.set_xlabel('date', fontsize=12)
    ax.set_ylabel('')
    ax.set_title('Time-scaled phylogeny', fontsize=14)
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





