from warnings import warn
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from beast_pype.date_utilities import decimal_to_date, date_to_decimal
from beast_pype.fig_utils import year_decimal_to_date_tick_labels
from scipy.interpolate import interp1d
from copy import deepcopy
import seaborn as sns
import re
from arviz import hdi


def read_log_file(file_path,
                  remove_0_first_sampling_prop=True,
                  youngest_tip_date=None,
                  convert_become_uninfectious_rate=False,
                  convert_growth_rate=True
                  ):
    """Read a log file into a pandas DataFrame.

    Parameters
    ----------------
    file_path: str
        Path to the log file.
    remove_0_first_sampling_prop : bool, default True
        If all the values in the first column starting with 'samplingProportion' are equal
        to 0, remove that column.
    youngest_tip_date : datetime, default None
        Youngest tip date. If supplied an attempt are made to:
         * determine the treeheight column and convert to TMRCA (through subtraction from youngest_tip_date).
         * determine the BDSKY Origin column and convert to Origin (through subtraction from youngest_tip_date).
    convert_become_uninfectious_rate:   bool, default False
        If true column 'becomeUninfectiousRate_BDSKY_Serial' will be used to
         calculate columns 'Rate of Becoming Uninfectious (per day)' and
          'Infection period (per day)'.
    convert_growth_rate: bool, default True
        If true column 'growthRate' will be used to calculate columns 'Growth Rate per day' and 'Doubling Time in days'.

    Returns
    -----------
    Pandas.DataFrame
    """
    if file_path.endswith('.log'):
        trace_df = pd.read_table(file_path, sep='\t', comment='#')
    elif file_path.endswith('.csv'):
        trace_df = pd.read_csv(file_path)
    elif file_path.endswith('.tsv'):
        trace_df = pd.read_csv(file_path, sep='\t')
    else:
        raise ValueError('file_path should end with .log, .csv or .tsv')
    if remove_0_first_sampling_prop:
        sampling_prop_columns = [column for column in trace_df.columns if column.startswith('samplingProportion')]
        if len(sampling_prop_columns) > 0 and (trace_df[sampling_prop_columns[0]] == 0).all():
            trace_df.drop(sampling_prop_columns[0], axis=1, inplace=True)
    if youngest_tip_date is not None:
        if not isinstance(youngest_tip_date, (int, float)):
            youngest_tip_date = date_to_decimal(youngest_tip_date)
        possible_origins = [column
                            for column in trace_df.columns
                            if 'origin' in column.lower()]
        if len(possible_origins) >1:
            warn(f'More than one column contains the substring "origin" (after conversion to lowercase) in trace at {file_path}.' +
                 'Therefore conversion to "Origin" column by subtraction from youngest tip date skipped.')
        if len(possible_origins) ==1:
            column_starting_origin = possible_origins[0]
            trace_df['Origin'] = youngest_tip_date - trace_df[column_starting_origin]
            # trace_df.drop(columns=column_starting_origin, inplace=True)
            # If this column is deleted it screws up the _gridded_skyline function
            # needed in plotting smoothed skylines. I do not have time to debug the
            # resulting break.

        possible_tree_heights = [column
                                    for column in trace_df.columns
                                    if 'treeheight' in re.sub(r'[^a-z]', '', column.lower())]
        if len(possible_tree_heights) >1:
            warn(f'More than one column contains the substring "treeheight" (after conversion to lowercase and removal of non-letters) in trace at {file_path}.' +
                 'Therefore conversion to "TMRCA" column by subtraction from youngest tip date skipped.')
        elif len(possible_tree_heights) ==1:
            tree_height = possible_tree_heights[0]
            trace_df['TMRCA'] = youngest_tip_date - trace_df[tree_height]
            # trace_df.drop(columns=tree_height, inplace=True)
            # Deleting this column caused a few bugs. Don't have time to sortout.
        else:
            warn(
                f'No column contains the substring "treeheight" (after conversion to lowercase and removal of non-letters) in trace at {file_path}.' +
                'Therefore conversion to "TMRCA" column by subtraction from youngest tip date skipped.')
        if convert_growth_rate and 'growthRate' in trace_df.columns:
            trace_df['Growth Rate per day'] = trace_df['growthRate'] / 365.25
            trace_df['Doubling Time in days'] = np.log(2) / trace_df['Growth Rate per day']
        
        if convert_become_uninfectious_rate:
            trace_df['Rate of Becoming Uninfectious (per day)'] = trace_df['becomeUninfectiousRate_BDSKY_Serial'] / 365
            trace_df['Mean Infection period (days)'] = 1 / trace_df['Rate of Becoming Uninfectious (per day)']
    return trace_df


def read_xml_set_logs_for_plotting(file_path_dict,
                                   xml_set_label='xml_set',
                                   remove_0_first_sampling_prop=True,
                                   convert_become_uninfectious_rate=False,
                                   youngest_tip_dates_dict=None):
    """
    Read xml set log files into a pandas DataFrame and a melted dataframe.

    Parameters
    ----------
    file_path_dict: dict {str: str}
        Dictionary of paths to log files. Keys are strain names.
    xml_set_label: str default 'xml_set'
        Label for xml_set divisions e.g. 'country' or 'strain'.
    remove_0_first_sampling_prop : bool, default True
        If all the values in the first column starting with 'samplingProportion' are equal
        to 0, remove that column.
    convert_become_uninfectious_rate:   bool, default False
        If true column 'becomeUninfectiousRate_BDSKY_Serial' will be used to
         calculate columns 'Rate of Becoming Uninfectious (per day)' and
          'Infection period (per day)'.
    youngest_tip_dates_dict: dict {str: float}
        Dictionary of youngest tip dates. Keys are xml_set names.


    Returns
    -------

    """
    dfs_list = []
    for xml_set, file_path in file_path_dict.items():
        if youngest_tip_dates_dict is not None:
            youngest_tip = youngest_tip_dates_dict[xml_set]
        else:
            youngest_tip = youngest_tip_dates_dict
        temp_df = read_log_file(file_path,
                                remove_0_first_sampling_prop=remove_0_first_sampling_prop,
                                youngest_tip_date=youngest_tip,
                                convert_become_uninfectious_rate=convert_become_uninfectious_rate)
        temp_df.insert(loc=0, column=xml_set_label, value=xml_set)
        dfs_list.append(temp_df)

    df = pd.concat(dfs_list)

    id_vars = [xml_set_label]
    value_vars = [col for col in df.columns if col not in id_vars]
    df_melt = df.melt(id_vars=id_vars,
                      value_vars=value_vars,
                      value_name='Estimate')
    return df, df_melt


def percentile_5th(g):
    return np.percentile(g, 5)

def percentile_95th(g):
    return np.percentile(g, 95)

def percentile_pivot(df, column, xml_set_label='xml_set'):
    """
    Create a percentile based pivot tale.
    
    """
    df_to_return = pd.pivot_table(df, values=[column], index=[xml_set_label],
                                  aggfunc= [percentile_5th, np.median, percentile_95th])
    if column == 'Origin':
        df_to_return = df_to_return.map(decimal_to_date).map(lambda x: x.strftime('%Y-%m-%d'))
    df_to_return.index = df_to_return.index.str.split('_', expand=True)
    df_to_return.reset_index(inplace=True)
    df_to_return.columns = ['Type of Strain', 'Strain', 'Sample Size','5th Percentile', 'Median (50th Percentile)', '95th Percentile']
    return df_to_return


def hdi_columns_starting_with(df, starting_with, prob=0.95):
    """
    Calculate HDI for columns starting with a given string.

    Parameters
    ----------
    df: pd.DataFrame
    starting_with: str
        C
    hdi_prob: float, default=0.95
        HDI probability to use.

    Returns
    -------
    pd.DataFrame

    """
    records = []
    cols_starting_with = [col for col in df.columns if col.startswith(starting_with)]
    for column in cols_starting_with:
        selected_value = df[column].to_numpy()
        lower_interval, upper_interval = hdi(selected_value, prob=prob)
        median_val = np.median(selected_value)
        records.append({'Parameter': column,
                        f'Lower {str(prob)} HDI': lower_interval,
                        'Median': median_val,
                        f'Upper {str(prob)} HDI': upper_interval})
    return pd.DataFrame.from_records(records)

def hdi_pivot(df, column, xml_set_label='xml_set', hdi_prob=0.95):
    """
    Create a hdi based pivot tale.

    Parameters
    ----------
    df: pd.DataFrame
        Data Frame to be pivoted.
    column: str
        Column to calculate HDI for.
    xml_set_label: str, optional
        Label for XML set, if not provided 'xml_set' is used.
    hdi_prob: float, default=0.95
        HDI probability to use.

    Returns
    -------
    pd.DataFrame
    """
    records = []
    for selection in df[xml_set_label].unique():
        selected_value = df[column][df[xml_set_label]==selection].to_numpy()
        lower_interval, upper_interval = hdi(selected_value, prob=hdi_prob)
        median_val = np.median(selected_value)
        records.append(
            {xml_set_label: selection,
             f'Lower {str(hdi_prob)} HDI': lower_interval,
             'Median':median_val,
             f'Upper {str(hdi_prob)} HDI': upper_interval
             })
    
    return pd.DataFrame.from_records(records)

def summary_stats_and_plot(df,
                           x,
                           y,
                           convert_plot_to_seconds=False,
                           include_grid=True,
                           violinplot_kwargs={},
                           boxplot_kwargs={}):
    """
    Calculate summary statistics and plot as box violin plots.

    Parameters
    ----------
    df: pandas.DataFrame
        Melted DataFrame for use in Seaborne.
    x: str
        Column to plot on x-axis.
    y: str
        Column to plot on y-axis.
    convert_plot_to_seconds: bool, default=False
        Whether to convert plot to seconds.
    include_grid: bool, default=True
        Whether to include grid.
    violinplot_kwargs: dict
        Keyword arguments to pass to seaborn.violinplot.
    boxplot_kwargs: dict
        Keyword arguments to pass to seaborn.boxplot.

    Returns
    -------
    fig: Seaborn plot
        Seaborn box violin plots.
    y_df: pandas.DataFrame
        DataFrame with columns specified in y and x.
        If convert_plot_to_second ==True an additional column where y is converted
        to seconds is included.
    summary_stats: pandas.DataFrame
        Summary statistics.
    """
    df = deepcopy(df)
    if convert_plot_to_seconds:
        df[f'{y} seconds'] = df[y].dt.total_seconds()
        columns = [x, f'{y} seconds', y]
        plot_y = f'{y} seconds'
    else:
        columns = [x, y]
        plot_y = y

    y_df = df[columns]
    summary_stats = y_df.groupby(x).describe().transpose()
    fig = plot_box_violin(df=y_df,
                          x=x,
                          y=plot_y,
                          include_grid=include_grid,
                          violinplot_kwargs=violinplot_kwargs,
                          boxplot_kwargs=boxplot_kwargs
                          )
    return fig, y_df, summary_stats



def plot_box_violin(df, x, y, include_grid=True, violinplot_kwargs={}, boxplot_kwargs={}):
    """
    Plot box violin plots.

    Parameters
    ----------
    df: pandas.DataFrame
        Melted DataFrame for use in Seaborne.
    x: str
        Column to plot on x-axis.
    y: str
        Column to plot on y-axis.
    include_grid: bool, default=True
        Whether to include grid.
    violinplot_kwargs: dict
        Keyword arguments to pass to seaborn.violinplot.
    boxplot_kwargs: dict
        Keyword arguments to pass to seaborn.boxplot.

    Returns
    -------
    Seaborn box violin plots.
    """
    fig = sns.violinplot(data=df, y=y, x=x, inner=None, hue=x, **violinplot_kwargs)
    sns.boxplot(data=df, y=y, x=x, ax=fig, color='grey', width=0.3, **boxplot_kwargs)
    fig.set(ylabel=y)
    fig.set_xticklabels(fig.get_xticklabels(), rotation=90)
    plt.tight_layout()
    if include_grid:
        plt.grid(linestyle='-.')
    return fig


def plot_comparative_box_violin(df_melted, parameter, xml_set_label='xml_set', prior_draws=None,
                                include_grid=True, violinplot_kwargs={}, boxplot_kwargs={}):
    """
    Plot box violin plots comparing xml_sets.

    Parameters
    ----------
    df_melted : pandas.DataFrame
        Melted DataFrame for use in Seaborne.  Contains columns 'xml_set', 'variable' and
        'Estimate'.
    parameter   :   str
        Parameter name.
    xml_set_label: str, optional
        Label for x-axis, if not provided this will be 'xml_set'.
    prior_draws :   numpy.ndarray,default None
        If given included in plot as 'Draws from Prior'
    include_grid:   bool, default=True
        Whether to include grid in violin plot.
    violinplot_kwargs: dict
        Keyword arguments to pass to seaborn.violinplot.
    boxplot_kwargs: dict
        Keyword arguments to pass to seaborn.boxplot.

    Returns
    -------
    Seaborn box violin plots.
    """
    df_melted = df_melted[df_melted['variable'] == parameter]
    if prior_draws is not None:
        temp_data = {xml_set_label: 'Draws from Prior',
                     'Sample': None,
                     'variable': parameter,
                     'Estimate': prior_draws}
        temp_df = pd.DataFrame(temp_data)
        df_melted = pd.concat([df_melted, temp_df])

    fig = plot_box_violin(df=df_melted,
                          x=xml_set_label,
                          y='Estimate',
                          include_grid=include_grid,
                          violinplot_kwargs=violinplot_kwargs,
                          boxplot_kwargs=boxplot_kwargs)
    return fig


def _set_dates_skyline_plotting_df(log_df, parameter_start, style, partition_year_decimals, hdi_prob=0.95):
    """Generate a DataFrame for plotting HDI values from a skyline model log.

    Parameters
    ---------------
    log_df: pd.DataFrame
        BDSKY log dataframe.
    hdi_prob: float
        Highest density interval to be used.


    Returns
    -----------
    parameter_plotting_df: pd.DataFrame
        HDI and median R_t values with deci dates (year fractions) and dates.
    """
    parameter_df = log_df.loc[:, log_df.columns.str.startswith(parameter_start)]
    parameter_df.dropna(axis=1, how='all', inplace=True)
    parameter_hdi = pd.DataFrame({label: hdi(series.to_numpy(), prob=hdi_prob) for label, series in parameter_df.items()})
    parameter_median = parameter_df.median()
    parameter_values = pd.concat([parameter_hdi.iloc[0, :], parameter_median, parameter_hdi.iloc[-1, :]], axis=1)
    parameter_values.columns = ['lower', 'median', 'upper']
    if parameter_df.columns[0].startswith('sampling') and not parameter_df.columns[0].endswith('.1'):
        # This only happens if the firs sampling proportion has been dropped in read_log_file
        # for being fixed at 0.
        year_decimals = pd.Series(np.sort(partition_year_decimals))
    else:
        start = log_df['Origin'].median()
        year_decimals = pd.Series(np.sort(np.append(partition_year_decimals, start)))
    mid_points = [np.mean([i, j])
                  for i, j in zip(year_decimals[:-1], year_decimals[1:])]
    if style.startswith('skyline'):
        parameter_plotting_df = pd.DataFrame({
            'year_decimal': year_decimals.repeat(2)[1:-1].to_list(),
            'lower': parameter_values['lower'].repeat(2).to_list(),
            'median': parameter_values['median'].repeat(2).to_list(),
            'upper': parameter_values['upper'].repeat(2).to_list()
        })
    elif style.startswith('smooth spline'):
        kind = 'quadratic'
        lower_spline = interp1d(mid_points, parameter_values['lower'], kind=kind)
        median_spline = interp1d(mid_points, parameter_values['median'], kind=kind)
        upper_spline = interp1d(mid_points, parameter_values['upper'], kind=kind)
        new_x_vals = np.linspace(mid_points[0], mid_points[-1], 1000)
        parameter_plotting_df = pd.DataFrame({
            'year_decimal': new_x_vals,
            'lower': lower_spline(new_x_vals),
            'median': median_spline(new_x_vals),
            'upper': upper_spline(new_x_vals)
        })
    else:
        raise ValueError('style start with "skyline" or "smooth spline".')
    if style.endswith('with mid-points'):
        mid_point_df = pd.DataFrame({
            'year_decimal': mid_points,
            'lower': parameter_values['lower'],
            'median': parameter_values['median'],
            'upper': parameter_values['upper']
        })
        return parameter_plotting_df, mid_point_df
    else:
        return parameter_plotting_df


def _gridded_skyline(trace_df,
                     youngest_tip,
                     parameter_start='reproductiveNumber',
                     grid_size=100):
    """
    Create gridded skyline parameters from Birth-Death Skyline trace over time.

    Parameters
    ----------
    trace_df : DataFrame
        Data frame of BEAST2 MCMC runs.
    youngest_tip : float
        Year decimal of youngest tip.
    parameter_start : str
        Starting string of parameters.
    grid_size: int, default 100
        Grid size for smoothing skyline interpolation.


    Returns
    -------
    Pandas.Dataframe

    Notes
    --------
    Skyline griding interpolation creation adapted from the R function gridSkyline
    from https://github.com/laduplessis/bdskytools/tree/master
    """
    parameter_subset = trace_df.loc[:, trace_df.columns.str.startswith(parameter_start)]
    origin_column = [column for column in trace_df.columns if column.startswith('origin')][0]
    origins = trace_df[origin_column]
    origin_med = origins.median()
    grid_times = np.linspace(start=0, stop=origin_med, num=grid_size)
    links, sky_grids = parameter_subset.shape
    sky_matrix = parameter_subset.to_numpy()

    # skyline_gridded creation adapted from the R function gridSkyline from
    # https://github.com/laduplessis/bdskytools/tree/master
    skyline_gridded = []
    for i in range(links):
        skyline_indices = np.maximum(1,
                                     sky_grids - np.floor(grid_times / origins[i] * sky_grids)
                                     ).astype(int) - 1
        skyline_gridded.append(sky_matrix[i, skyline_indices])

    skyline_gridded = pd.DataFrame(skyline_gridded, columns=grid_times)
    skyline_gridded_hpd = skyline_gridded.apply(lambda x: hdi(x.to_numpy(),
                                                              prob=0.95), axis=0, result_type='expand')
    skyline_gridded_hpd = skyline_gridded_hpd.transpose()
    skyline_gridded_hpd.columns = ['lower', 'upper']
    skyline_gridded_hpd.insert(1, 'median', skyline_gridded.median())
    skyline_gridded_hpd['year_decimal'] = youngest_tip - grid_times
    return skyline_gridded_hpd

def _plot_df_to_hdi_df(plot_df):
    hdi_est_df = plot_df[plot_df.index % 2 == 0]
    end_df = plot_df[plot_df.index % 2 != 0]
    hdi_est_df['Start of Period'] = hdi_est_df.year_decimal.map(decimal_to_date).map(
        lambda x: x.strftime('%Y-%m-%d'))
    hdi_est_df['End of Period'] = end_df.year_decimal.map(decimal_to_date).map(
        lambda x: x.strftime('%Y-%m-%d')).to_list()
    return hdi_est_df[['Start of Period', 'End of Period', 'lower', 'median', 'upper']]
    
    

def plot_skyline(traces_df,
                 youngest_tip_year_decimal,
                 parameter_start='reproductiveNumber',
                 x_tick_freq='automatic',
                 style='skyline',
                 y_label='$R_t$',
                 palette=sns.color_palette(),
                 include_grid=True,
                 grid_size=100,
                 partition_year_decimals=None,
                 xml_set_label='xml_set'
                 ):
    """
    Plot parameters skylines with set dates from Birth-Death Skyline traces over time.

    Parameters
    ----------
    traces_df : DataFrame
        Data frame of BEAST2 MCMC runs.
    partition_year_decimals : Series or list of dates
        Dates of partitioning skyline.
    youngest_tip_year_decimal : float
        Year decimal of youngest tip.
    parameter_start : str
        Starting string of parameters.
    x_tick_freq : str, default 'automatic'
        Frequency of x ticks. Options are 'automatic', 'yearly', 'quarterly', 'monthly', 'half month' or 'weekly'.
    style : str, default 'skyline'
        Style of plotting.
    y_label : str, default '$R_t$'
        Label of y-axis.
    palette : str, default sns.color_palette()
        Palette of plotting.
    include_grid : bool, default True
        Whether to include grid lines.
    grid_size : int, default 100
        Grid size for smoothing skyline interpolation.
    xml_set_label : str, default None
        Title for legend, if not provided 'xml_set' is used instead.

    Returns
    -------
    fig : Matplotlib figure
        Matplotlib figure.
    ax : Matplotlib axes
        Matplotlib axes.
    hdi_df : pd.DataFrame
        Table of HDIs with dates of boundaries.
    """
    if xml_set_label in traces_df.columns:
        xml_sets = traces_df[xml_set_label].unique()
    else:
        xml_sets = []
    if len(xml_sets) > 1:
        if partition_year_decimals is None:
            partition_year_decimals = {xml_set: None for xml_set in xml_sets}
        for name, value in {'partition_year_decimals': partition_year_decimals,
                            'youngest_tip_date': youngest_tip_year_decimal}.items():
            if not isinstance(value, dict):
                raise TypeError(name + ' must be a dict if there are multiple xml_sets.')
            if any(stain not in xml_sets for stain in value.keys()):
                raise ValueError('Each key of ' + name + ' should be the name of a xml_set.')
            if any(strain not in value.keys() for strain in xml_sets):
                raise ValueError('All of the xml_sets should be represented in ' + name + '.')
        plot_dfs = {}
        hdi_df = []
        for xml_set in xml_sets:
            if partition_year_decimals[xml_set] is not None:
                plot_df = _set_dates_skyline_plotting_df(traces_df[traces_df[xml_set_label] == xml_set],
                                                         parameter_start=parameter_start,
                                                         style=style,
                                                         partition_year_decimals=partition_year_decimals[xml_set])
                plot_dfs[xml_set] = deepcopy(plot_df)
                hdi_est = _plot_df_to_hdi_df(plot_df=plot_df)
                hdi_est.insert(0, xml_set_label, xml_set)
                hdi_df.append(hdi_est)
            else:
                if style != 'skyline':
                    raise ValueError('Currently only "skyline" style is supported when plotting without partition dates.')
                plot_df = _gridded_skyline(
                        traces_df[traces_df[xml_set_label] == xml_set],
                        parameter_start=parameter_start,
                        youngest_tip=youngest_tip_year_decimal[xml_set],
                        grid_size=grid_size)
                plot_dfs[xml_set] = deepcopy(plot_df)
                plot_df.insert(loc=0, column=xml_set_label, value=xml_set)
                hdi_df.append(plot_df)
        hdi_df = pd.concat(hdi_df)
    else:
        if partition_year_decimals is not None:
            plot_dfs = _set_dates_skyline_plotting_df(traces_df,
                                                      parameter_start=parameter_start,
                                                      style=style,
                                                      partition_year_decimals=partition_year_decimals)
            hdi_df = _plot_df_to_hdi_df(plot_dfs)
        else:
            if style != 'skyline':
                raise ValueError('Currently only "skyline" style is supported when plotting without partition dates.')
            plot_dfs = _gridded_skyline(
                traces_df,
                parameter_start=parameter_start,
                youngest_tip=youngest_tip_year_decimal,
                grid_size=grid_size)
            hdi_df = deepcopy(plot_dfs)

    if style.endswith('with mid-points'):
        x_vals = pd.concat([item[0] for item in plot_dfs.values()])['year_decimal']
    else:
        if isinstance(plot_dfs, dict):
            x_vals = pd.concat(plot_dfs.values())['year_decimal']
        else:
            x_vals = plot_dfs['year_decimal']
    tick_year_decimals, tick_labels = year_decimal_to_date_tick_labels(x_vals, tick_freq=x_tick_freq)
    fig, ax = plt.subplots(1, 1, layout='constrained', figsize=(7.5, 5))
    if not isinstance(plot_dfs, dict):
        plot_dfs = {None: plot_dfs}
    for i, (xml_set, plot_df) in enumerate(plot_dfs.items()):
        if style.endswith('with mid-points'):
            mid_point_df = plot_df[1]
            plot_df = plot_df[0]
            ax.scatter(mid_point_df['year_decimal'], mid_point_df['upper'], marker='v', c=palette[i])
            ax.scatter(mid_point_df['year_decimal'], mid_point_df['median'],marker='D', c=palette[i])
            ax.scatter(mid_point_df['year_decimal'], mid_point_df['lower'], marker="^", c=palette[i])

        ax.plot(plot_df['year_decimal'], plot_df['median'], color=palette[i], label=xml_set)
        ax.fill_between(plot_df['year_decimal'], plot_df['lower'], plot_df['upper'], color=palette[i],
                        alpha=0.2)

    ax.xaxis.set_ticks(tick_year_decimals)
    ax.set_xticklabels(tick_labels)
    ax.tick_params(axis='x', labelrotation=45)
    plt.setp(ax.get_xticklabels(), ha='right')
    ax.set_xlabel('Date')
    ax.set_ylabel(y_label)
    if len(xml_sets) > 1:
        plt.legend()
    plt.tight_layout()
    if include_grid:
        plt.grid(linestyle='-.')

    return fig, ax, hdi_df

def outputs_for_possible_skyline(trace_df, youngest_tip_year_decimal, param_start, axis_label, grid_size=100, include_grid=True):
    """
    Determine if output figure should be a skyline if so plot with

    Parameters
    ----------
    trace_df : DataFrame
        Data frame of BEAST2 MCMC runs.
    youngest_tip_year_decimal : float
        Year decimal of youngest tip.
    param_start : str
        Starting string of parameters.
    axis_label:
    include_grid : bool, default True
        Whether to include grid lines.
    grid_size : int, default 100
        Grid size for smoothing skyline interpolation.

    Returns
    -------
    fig: matploplib figure
     ax: matploplib.axes
     hdi: pd.DataFrame
        Table of HD intervals
    """
    cols_starting_with = [col for col in trace_df.columns if col.startswith(param_start)]
    if len(cols_starting_with) > 1:
        fig, ax = plot_skyline(trace_df,
                               youngest_tip_year_decimal,
                               parameter_start=param_start,
                               y_label=axis_label,
                               grid_size=grid_size,
                               include_grid=include_grid)
        hdi = hdi_columns_starting_with(trace_df, param_start)
    elif len(cols_starting_with) == 1:
        fig, ax, hdi = plot_hist_kde(trace_df=trace_df, parameter=param_start, hdi_prob=0.95, x_label=axis_label)
    else:
        raise ValueError(f'Trace file does not contain any columns starting with "{param_start}".')
    return fig, ax, hdi


def plot_hist_kde(trace_df,
                  parameter,
                  x_label=None,
                  color=None,
                  hdi_prob=0.95,
                  tight_layout=True):
    """
    Plot histogram with kde.

    Parameters
    ----------
    trace_df: pandas.DataFrame
        DataFrame with column headed with parameters.
    parameter: str
        Parameter name.
    x_label: str, default None
        Label of x-axis. If None, use parameters name.
    color: str, default None
        Color of histogram and kde line.
    hdi_prob: float, default 0.95
        Prob for which the highest density interval will be computed. If not none
        lower and HDI lines plotted.
    tight_layout : bool, default True
        Whether to use tight layout or not.

    Returns
    -------
    if hdi_prob is None:
        Matplotlib figure and axis.
    else:
        Matplotlib figure, axis and dictionary of 'lower' and 'upper' hdi with median.
    """
    if x_label is None:
        x_label = parameter
    fig, ax = plt.subplots()  # initializes figure and plots
    sns.histplot(trace_df[parameter], ax=ax, kde=True, color=color)
    if hdi_prob is not None:
        hdi_est = hdi(trace_df[parameter].to_numpy(), prob=hdi_prob)
        upper_key = f'Upper {str(hdi_prob)} HDI'
        lower_key = f'Lower {str(hdi_prob)} HDI'
        hdi_est = {lower_key: float(hdi_est[0]),
                   'Median': float(trace_df[parameter].median()),
                   upper_key: float(hdi_est[1])}
        ax.axvline(hdi_est['Median'], color='k', lw=2)
        ax.axvline(hdi_est[lower_key], color='k', ls='--', lw=1)
        ax.axvline(hdi_est[upper_key], color='k', ls='--', lw=1)
    ax.set_xlabel(x_label)
    if tight_layout:
        plt.tight_layout()

    if hdi_prob is None:
        return fig, ax
    else:
        return fig, ax, hdi_est


def plot_origin_or_tmrca(trace_df, parameter, x_tick_freq='automatic', hdi_prob=None, color=None):
    """Plot histograms of origin or TMRCA, with kde line.

    Parameters
    ----------
    trace_df: pandas.DataFrame
        DataFrame with column headed with parameters.
    parameter: str
        Parameter name.
    x_tick_freq : str, default 'automatic'
        Frequency of x ticks. Options are 'automatic', 'yearly', 'quarterly', 'monthly', 'half month' or 'weekly'.
    hdi_prob: float, default None
        Prob for which the highest density interval will be computed. If not none
        lower and HDI lines plotted.
    color: str or RGB tuple, default None
        Color of histogram and kde line.

    Returns
    -------
    if hdi_prob is None:
        Matplotlib figure and axis.
    else:
        Matplotlib figure, axis and dictionary of 'lower' and 'upper' hdi with median.
    """
    if parameter not in ['Origin', 'TMRCA']:
        raise ValueError('parameter must be "Origin" or "TMRCA".')
    tick_year_decimals, tick_labels = year_decimal_to_date_tick_labels(trace_df[parameter], tick_freq=x_tick_freq)
    output = plot_hist_kde(trace_df,
                           parameter,
                           x_label=f'{parameter} (date)',
                           color=color,
                           hdi_prob=hdi_prob,
                           tight_layout=False)
    if hdi_prob is None:
        fig, ax = output
    else:
        fig, ax, hdi_est = output
    ax.xaxis.set_ticks(tick_year_decimals)
    ax.set_xticklabels(tick_labels)
    ax.tick_params(axis='x', labelrotation=45)
    plt.setp(ax.get_xticklabels(), ha='right')
    plt.tight_layout()
    if hdi_prob is None:
        return fig, ax
    else:
        hdi_est = {key: decimal_to_date(value) for key, value in hdi_est.items()}
        return fig, ax, hdi_est


def plot_comparative_origin_or_tmrca(df_melted, parameter, xml_set_label='xml_set',
                                     x_tick_freq='automatic', one_figure=False, palette=None):
    """Plot histograms of origins for each strain, with kde line.

    Parameters
    ----------
    df_melted : pandas.DataFrame
        Melted DataFrame for use in Seaborne.  Contains columns 'xml_set', 'variable' and
        'Estimate'.
    parameter : str
        Parameter to plot must be either 'Origin' or 'TMRCA'.
    xml_set_label: str, default None
        Legend title to use, if none provided this will be 'xml_set'.
    x_tick_freq: str, default 'automatic'
        Suggested x tick frequency. Options are 'automatic', 'yearly', 'quarterly', 
        'monthly', 'half month' or 'weekly'.
    one_figure: bool, default False
        If True, then one figure will be created with the histograms for each strain
        overlaid on top of each other.
        If False, then sub-figures will be created with a histogram for each strain per
        row.
    palette : list or dict of RGB colors (tuples), default None
        List of colors to use for each strain.
        Or dict of RGB colors to use for each strain.

    Returns
    -------
    If one_figure is True, then a single figure will be created.
    If one_figure is False, then a Seaborn.FactGrid will be created.
    """
    if parameter not in ['Origin', 'TMRCA']:
        raise ValueError('parameter must be "Origin" or "TMRCA".')
    df = df_melted[df_melted.variable == parameter]
    df.rename(columns={'Estimate': parameter}, inplace=True)
    tick_year_decimals, tick_labels = year_decimal_to_date_tick_labels(df[parameter], tick_freq=x_tick_freq)
    if one_figure:
        ax = sns.histplot(data=df, x=parameter, hue=xml_set_label, kde=True, palette=palette)
        if palette is None:
            palette = sns.color_palette()
        for colour, xml_set in zip(palette, df[xml_set_label].unique()):
            strain_df = df[df[xml_set_label] == xml_set]
            ax.axvline(strain_df[parameter].median(), color=colour)
            ax.axvline(strain_df[parameter].quantile(0.05), color=colour, ls='--', lw=1)
            ax.axvline(strain_df[parameter].quantile(0.95), color=colour, ls='--', lw=1)
        fig = ax
    else:
        fig = sns.FacetGrid(df, row=xml_set_label, hue=xml_set_label, margin_titles=True, aspect=4)
        fig.map_dataframe(sns.histplot, x=parameter, kde=True)
        fig.set_titles(row_template='')
        fig.add_legend()
        for ax, xml_set in zip(fig.axes.flat, df[xml_set_label].unique()):
            strain_df = df[df[xml_set_label] == xml_set]
            ax.axvline(strain_df[parameter].median(), color='k', lw=2)
            ax.axvline(strain_df[parameter].quantile(0.05), color='k', ls='--', lw=1)
            ax.axvline(strain_df[parameter].quantile(0.95), color='k', ls='--', lw=1)

    ax.xaxis.set_ticks(tick_year_decimals)
    ax.set_xticklabels(tick_labels)
    ax.tick_params(axis='x', labelrotation=45)
    plt.setp(ax.get_xticklabels(), ha='right')
    plt.tight_layout()
    return fig


def read_metadata(metadata_path, collection_date_field):
    """Read a metadata CSV/TSV file and parse the collection date field.

    Parameters
    ----------
    metadata_path : str
        Path to the metadata file. Must end with ``.csv`` or ``.tsv``.
    collection_date_field : str
        Name of the column containing sample collection dates. This column is
        parsed to ``datetime`` (unparseable values become ``NaT``).

    Returns
    -------
    pandas.DataFrame
        The metadata with ``collection_date_field`` cast to datetime.
    """
    if metadata_path.endswith('.tsv'):
        sep = '\t'
    elif metadata_path.endswith('.csv'):
        sep = ','
    else:
        raise ValueError(
            'Only .csv and .tsv files are supported for metadata. '
            f'Got: {metadata_path}'
        )
    metadata = pd.read_csv(metadata_path, sep=sep)
    if collection_date_field not in metadata.columns:
        raise ValueError(
            f"collection_date_field '{collection_date_field}' not found in "
            f"{metadata_path}. Available columns: {', '.join(metadata.columns)}"
        )
    metadata[collection_date_field] = pd.to_datetime(
        metadata[collection_date_field], errors='coerce')
    return metadata


def describe_collection_dates(metadata, collection_date_field):
    """Summarise a collection date column via :meth:`pandas.Series.describe`.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata containing the collection date column (as datetime).
    collection_date_field : str
        Name of the collection date column.

    Returns
    -------
    pandas.DataFrame
        A single column DataFrame of descriptive statistics for the collection
        dates.
    """
    dates = pd.to_datetime(metadata[collection_date_field], errors='coerce')
    description = dates.describe()
    return description.to_frame(name=collection_date_field)


def describe_collection_dates_by_xml_set(metadata_dict,
                                         collection_date_field,
                                         xml_set_label='xml_set'):
    """Summarise collection dates for several xml sets side by side.

    Parameters
    ----------
    metadata_dict : dict of {str: pandas.DataFrame}
        Mapping of xml set name to its metadata DataFrame.
    collection_date_field : str
        Name of the collection date column.
    xml_set_label : str, default 'xml_set'
        Name to use for the columns index of the returned DataFrame.

    Returns
    -------
    pandas.DataFrame
        Descriptive statistics with one column per xml set.
    """
    describes = {}
    for xml_set, metadata in metadata_dict.items():
        dates = pd.to_datetime(metadata[collection_date_field], errors='coerce')
        describes[xml_set] = dates.describe()
    summary = pd.DataFrame(describes)
    summary.columns.name = xml_set_label
    return summary


def plot_collection_date_histogram(metadata, collection_date_field, bins=30,
                                   ax=None, x_tick_freq='automatic'):
    """Plot a histogram of sample collection dates.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata containing the collection date column (as datetime).
    collection_date_field : str
        Name of the collection date column.
    bins : int, default 30
        Number of histogram bins.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None a new figure and axes are created.
    x_tick_freq : str, default 'automatic'
        Suggested x-axis tick frequency. Options are 'automatic', 'yearly',
        'quarterly', 'monthly', 'half month' or 'weekly'.

    Returns
    -------
    tuple of (matplotlib.figure.Figure, matplotlib.axes.Axes)
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        fig = ax.get_figure()
    dates = pd.to_datetime(metadata[collection_date_field], errors='coerce').dropna()
    year_decimals = dates.map(date_to_decimal)
    sns.histplot(x=year_decimals, bins=bins, ax=ax)
    ax.set_xlabel('Collection date')
    ax.set_ylabel('Number of samples')
    if len(year_decimals) > 0:
        tick_year_decimals, tick_labels = year_decimal_to_date_tick_labels(
            year_decimals, tick_freq=x_tick_freq)
        ax.xaxis.set_ticks(tick_year_decimals)
        ax.set_xticklabels(tick_labels)
        ax.tick_params(axis='x', labelrotation=45)
        plt.setp(ax.get_xticklabels(), ha='right')
    plt.tight_layout()
    return fig, ax


def plot_stacked_collection_date_histogram(metadata_dict, collection_date_field,
                                           xml_set_label='xml_set', bins=30,
                                           ax=None, x_tick_freq='automatic'):
    """Plot a stacked histogram of collection dates for several xml sets.

    Parameters
    ----------
    metadata_dict : dict of {str: pandas.DataFrame}
        Mapping of xml set name to its metadata DataFrame.
    collection_date_field : str
        Name of the collection date column.
    xml_set_label : str, default 'xml_set'
        Label used for the grouping (legend title).
    bins : int, default 30
        Number of histogram bins.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None a new figure and axes are created.
    x_tick_freq : str, default 'automatic'
        Suggested x-axis tick frequency. Options are 'automatic', 'yearly',
        'quarterly', 'monthly', 'half month' or 'weekly'.

    Returns
    -------
    tuple of (matplotlib.figure.Figure, matplotlib.axes.Axes)
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        fig = ax.get_figure()
    frames = []
    for xml_set, metadata in metadata_dict.items():
        dates = pd.to_datetime(metadata[collection_date_field],
                               errors='coerce').dropna()
        entry = pd.DataFrame({
            'year_decimal': dates.map(date_to_decimal),
            xml_set_label: xml_set,
        })
        frames.append(entry)
    combined = pd.concat(frames, ignore_index=True)
    sns.histplot(data=combined, x='year_decimal', hue=xml_set_label,
                 multiple='stack', bins=bins, ax=ax)
    ax.set_xlabel('Collection date')
    ax.set_ylabel('Number of samples')
    if len(combined) > 0:
        tick_year_decimals, tick_labels = year_decimal_to_date_tick_labels(
            combined['year_decimal'], tick_freq=x_tick_freq)
        ax.xaxis.set_ticks(tick_year_decimals)
        ax.set_xticklabels(tick_labels)
        ax.tick_params(axis='x', labelrotation=45)
        plt.setp(ax.get_xticklabels(), ha='right')
    plt.tight_layout()
    return fig, ax