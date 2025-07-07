from matplotlib import pyplot as plt
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from arviz import hdi
from beast_pype.date_utilities import decimal_to_date, date_to_decimal
from scipy.interpolate import interp1d
from copy import deepcopy
import seaborn as sns

# stop annoying panadas warnings
pd.options.mode.chained_assignment = None



def read_log_file(file_path, cut_to_first=None):
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


def read_strain_logs_for_plotting(file_path_dict,
                                  convert_become_uninfectious_rate=False,
                                  youngest_tips_dict=None):
    """

    Parameters
    ----------
    file_path_dict: dict {str: str}
        Dictionary of paths to log files. Keys are strain names.
    convert_become_uninfectious_rate:   bool, default False
        If true column 'becomeUninfectiousRate_BDSKY_Serial' will be used to
         calculate columns 'Rate of Becoming Uninfectious (per day)' and
          'Infection period (per day)'.
    youngest_tips_dict: dict {str: float}
        Dictionary of paths to log files. Keys are strain names.


    Returns
    -------

    """
    dfs_list = []
    for xml_set, file_path in file_path_dict.items():
        if file_path.endswith('.log'):
            temp_df = read_log_file(file_path)
        elif file_path.endswith('.csv'):
            temp_df = pd.read_csv(file_path)
        elif file_path.endswith('.tsv'):
            temp_df = pd.read_csv(file_path, sep='\t')
        else:
            raise ValueError('Only csv, tsv or log files are supported.')
        if youngest_tips_dict is not None:
            youngest_tip = youngest_tips_dict[xml_set]
            if not isinstance(youngest_tip, (int, float)):
                youngest_tip = date_to_decimal(youngest_tip)
            temp_df['Origin'] = youngest_tip - temp_df['origin_BDSKY_Serial']
        if convert_become_uninfectious_rate:
            temp_df['Rate of Becoming Uninfectious (per day)'] = temp_df['becomeUninfectiousRate_BDSKY_Serial'] / 365
            temp_df['Infection period (per day)'] = 1 / temp_df['Rate of Becoming Uninfectious (per day)']

        temp_df.insert(loc=0, column='xml_set', value=xml_set)
        dfs_list.append(temp_df)

    df = pd.concat(dfs_list)

    id_vars = ['xml_set']
    value_vars = [col for col in df.columns if col not in id_vars]
    df_melt = df.melt(id_vars=id_vars,
                      value_vars=value_vars,
                      value_name='Estimate')
    return df, df_melt


def percentile_5th(g):
    return np.percentile(g, 5)

def percentile_95th(g):
    return np.percentile(g, 95)

def percentile_pivot(df, column):
    """
    Create a percentile based pivot tale.
    
    """
    df_to_return = pd.pivot_table(df, values=[column], index=['xml_set'],
                                  aggfunc= [percentile_5th, np.median, percentile_95th])
    if column == 'Origin':
        df_to_return = df_to_return.map(decimal_to_date).map(lambda x: x.strftime('%Y-%m-%d'))
    df_to_return.index = df_to_return.index.str.split('_', expand=True)
    df_to_return.reset_index(inplace=True)
    df_to_return.columns = ['Type of Strain', 'Strain', 'Sample Size','5th Percentile', 'Median (50th Percentile)', '95th Percentile']
    return df_to_return


def hdi_columns_starting_with(df, starting_with, hdi_prob=0.95):
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
        lower_interval, upper_interval = hdi(selected_value, hdi_prob=hdi_prob)
        median_val = np.median(selected_value)
        records.append({'Parameter': column,
                        f'Lower {str(hdi_prob)} HDI': lower_interval,
                        'Median': median_val,
                        f'Upper {str(hdi_prob)} HDI': upper_interval})
    return pd.DataFrame.from_records(records)

def hdi_pivot(df, column, hdi_prob=0.95):
    """
    Create a hdi based pivot tale.

    Parameters
    ----------
    df: pd.DataFrame
        Data Frame to be pivoted.
    column: str
        Column to calculate HDI for.
    hdi_prob: float, default=0.95
        HDI probability to use.

    Returns
    -------
    pd.DataFrame
    """
    records = []
    for selection in df['xml_set'].unique():
        selected_value = df[column][df['xml_set']==selection].to_numpy()
        lower_interval, upper_interval = hdi(selected_value, hdi_prob=hdi_prob)
        median_val = np.median(selected_value)
        records.append(
            {'xml_set': selection,
             f'Lower {str(hdi_prob)} HDI': lower_interval,
             'Median':median_val,
             f'Upper {str(hdi_prob)} HDI': upper_interval
             })
    
    return pd.DataFrame.from_records(records)


def plot_comparative_box_violine(df_melted, parameter, prior_draws=None,
                                 include_grid=True):
    """
    Plot box violin plots comparing strains.

    Parameters
    ----------
    df_melted : pandas.DataFrame
        DataFrame melted for use in Seaborne. Contains columns 'xml_set', 'variable' and
        'Estimate'.
    parameter   :   str
        Parameter name.
    prior_draws :   numpy.ndarray,default None
        If given included in plot as 'Draws from Prior'

    Returns
    -------
    Seaborn box violin plots.
    """
    df_melted = df_melted[df_melted['variable'] == parameter]
    if prior_draws is not None:
        temp_data = {'xml_set': 'Draws from Prior',
                     'Sample': None,
                     'variable': parameter,
                     'Estimate': prior_draws}
        temp_df = pd.DataFrame(temp_data)
        df_melted = pd.concat([df_melted, temp_df])

    fig = sns.violinplot(data=df_melted, y='Estimate', x='xml_set', inner=None, hue='xml_set')
    sns.boxplot(data=df_melted, y='Estimate', x='xml_set', ax=fig, color='grey', width=0.3)
    fig.set(ylabel=parameter)
    fig.set_xticklabels(fig.get_xticklabels(), rotation=90)
    plt.tight_layout()
    if include_grid:
        plt.grid(linestyle='-.')
    return fig


def _set_dates_skyline_plotting_df(log_df, parameter_start, style, youngest_tip, hdi_prob=0.95, dates_of_change=None):
    """Generate a DataFrame for plotting HDI values from a skyline model log.

    Parameters
    ---------------
    log_df: pd.DataFrame
        BDSKY log dataframe.
    youngest_tip: float/int or dt.datetime
        Date of youngest tip. If date or datetime will be converted to decimal year float.
    hdi_prob: float
        Highest density interval to be used.


    Returns
    -----------
    rt_plotting_df: pd.DataFrame
        HDI and median R_t values with deci dates (year fractions) and dates.
    """
    r_t = log_df.loc[:, log_df.columns.str.startswith(parameter_start)]
    r_t.dropna(axis=1, how='all', inplace=True)
    r_t_hdi = pd.DataFrame({label: hdi(series.to_numpy(), hdi_prob=hdi_prob) for label, series in r_t.items()})
    r_t_median = r_t.median()
    r_t_values = pd.concat([r_t_hdi.iloc[0, :], r_t_median, r_t_hdi.iloc[-1, :]], axis=1)
    r_t_values.columns = ['lower', 'median', 'upper']
    if not isinstance(youngest_tip, (int, float)):
        youngest_tip_date = deepcopy(youngest_tip)
        youngest_tip = date_to_decimal(youngest_tip)

    if not isinstance(youngest_tip, (int, float)):
        raise ValueError('youngest_tip must be int, float, date or datetime')
    median_origin = log_df.origin_BDSKY_Serial.median()
    median_origin_corrected = youngest_tip - median_origin
    nr_t = r_t.shape[1]
    if dates_of_change is None:
        year_decimals = np.linspace(median_origin_corrected, youngest_tip, nr_t + 1)
        year_decimals = pd.Series(year_decimals)
    else:
        if isinstance(dates_of_change, pd.DatetimeIndex):
            dates_of_change = dates_of_change.to_series()
        if isinstance(dates_of_change, (tuple, list, np.ndarray)):
            dates_of_change = pd.Series(dates_of_change).sort_values()
        median_origin_corrected = pd.Series(median_origin_corrected)
        year_decimals = dates_of_change.map(date_to_decimal)
        youngest_tip = pd.Series([youngest_tip])
        year_decimals = pd.concat([median_origin_corrected,
                                   year_decimals,
                                   youngest_tip])
    mid_points = [np.mean([i, j])
                  for i, j in zip(year_decimals[:-1], year_decimals[1:])]
    if style.startswith('skyline'):
        rt_plotting_df = pd.DataFrame({
            'year_decimal': year_decimals.repeat(2)[1:-1].to_list(),
            'lower': r_t_values['lower'].repeat(2).to_list(),
            'median': r_t_values['median'].repeat(2).to_list(),
            'upper': r_t_values['upper'].repeat(2).to_list()
        })
    elif style.startswith('smooth spline'):
        kind = 'quadratic'
        lower_spline = interp1d(mid_points, r_t_values['lower'], kind=kind)
        median_spline = interp1d(mid_points, r_t_values['median'], kind=kind)
        upper_spline = interp1d(mid_points, r_t_values['upper'], kind=kind)
        new_x_vals = np.linspace(mid_points[0], mid_points[-1], 1000)
        rt_plotting_df = pd.DataFrame({
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
            'lower': r_t_values['lower'],
            'median': r_t_values['median'],
            'upper': r_t_values['upper']
        })
        return rt_plotting_df, mid_point_df
    else:
        return rt_plotting_df


def _plot_gridded_skyline(trace_df,
                 youngest_tip,
                 parameter_start='reproductiveNumber',
                 grid_size=100):
    """
    Plot skyline parameters from Birth-Death Skyline trace over time.

    Parameters
    ----------
    trace_df : DataFrame
        Data frame of BEAST2 MCMC runs.
    youngest_tip : float
        Year decimal of youngest tip.
    parameter_start : str
        Starting string of parameters.
    x_tick_freq : str, default 'monthly'
        Frequency of x ticks.
    ylabel : str, default '$R_t$'
        Label of y-axis.
    color : str, default 'blue'
        Color to use for plotting, see matplotlib.colors module
         (https://matplotlib.org/stable/users/explain/colors/colors.html#colors-def).
    grid_size: int, default 100
        Grid size for smoothing skyline interpolation.
    include_grid : bool, default True
        Whether to include grid lines.
    ax : matplotlib.axis, default None
        Matplotlib axis to plot on.

    Returns
    -------
    Matplotlib figure and axes.

    Notes
    --------
    Skyline griding interpolation creation adapted from the R function gridSkyline
    from https://github.com/laduplessis/bdskytools/tree/master
    """
    parameter_subset = trace_df.loc[:, trace_df.columns.str.startswith(parameter_start)]
    origin_column = [column for column in trace_df.columns if column.startswith('origin_BDSKY')][0]
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
                                                              hdi_prob=0.95), axis=0, result_type='expand')
    skyline_gridded_hpd = skyline_gridded_hpd.transpose()
    skyline_gridded_hpd.columns = ['lower', 'upper']
    skyline_gridded_hpd.insert(1, 'median', skyline_gridded.median())
    skyline_gridded_hpd['year_decimal'] = youngest_tip - grid_times
    return skyline_gridded_hpd

def plot_skyline(traces_df,
                 youngest_tip,
                 parameter_start='reproductiveNumber',
                 x_tick_freq='monthly',
                 style='skyline',
                 ylabel='$R_t$',
                 palette=sns.color_palette(),
                 include_grid=True,
                 grid_size=100,
                 partition_dates=None):
    """
    Plot comparative parameters skylines with set dates from Birth-Death Skyline traces over time.

    Parameters
    ----------
    traces_df : DataFrame
        Data frame of BEAST2 MCMC runs.
    partition_dates : Series or list of dates
        Dates of partitoning skyline.
    youngest_tip : float
        Year decimal of youngest tip.
    parameter_start : str
        Starting string of parameters.
    x_tick_freq : str, default 'monthly'
        Frequency of x ticks.
    style : str, default 'skyline'
        Style of plotting.
    ylabel : str, default '$R_t$'
        Label of y-axis.
    palette : str, default sns.color_palette()
        Palette of plotting.
    include_grid : bool, default True
        Whether to include grid lines.

    Returns
    -------
    Matplotlib figure and axes.
    """
    if 'xml_set' in traces_df.columns:
        strains = traces_df.xml_set.unique()
    else:
        strains = []
    if len(strains) > 1:
        if partition_dates is not None:
            for name, value in {'partition_dates': partition_dates,
                                'youngest_tip': youngest_tip}.items():
                if not isinstance(value, dict):
                    raise TypeError(name + ' must be a dict if there are multiple strains.')
                if any(stain not in strains for stain in value.keys()):
                    raise ValueError('Each key of ' + name + ' should be the name of a strain.')
                if any(strain not in value.keys() for strain in strains):
                    raise ValueError('All of the stains should be represented in ' + name + '.')
            plot_dfs = {strain: _set_dates_skyline_plotting_df(traces_df[traces_df['xml_set'] == strain],
                                                                   parameter_start=parameter_start,
                                                                   style=style,
                                                                   youngest_tip=youngest_tip[strain],
                                                                   dates_of_change=partition_dates[strain])
                            for strain in strains}
        else:
            if style != 'skyline':
                raise ValueError('Currently only "skyline" style is supported when plotting without partition dates.')
            plot_dfs = {
                strain: _plot_gridded_skyline(
                    traces_df[traces_df['xml_set'] == strain],
                    parameter_start=parameter_start,
                    youngest_tip=youngest_tip[strain],
                    grid_size=grid_size)
                for strain in strains}
    else:
        if partition_dates is not None:
            plot_dfs = {strains[0]: _set_dates_skyline_plotting_df(traces_df,
                                                                       parameter_start=parameter_start,
                                                                       style=style,
                                                                       youngest_tip=youngest_tip,
                                                                       dates_of_change=partition_dates)}
        else:
            if style != 'skyline':
                raise ValueError('Currently only "skyline" style is supported when plotting without partition dates.')
            plot_dfs = _plot_gridded_skyline(
                traces_df,
                parameter_start=parameter_start,
                youngest_tip=youngest_tip,
                grid_size=grid_size)

    if style.endswith('with mid-points'):
        x_vals = pd.concat([item[0] for item in plot_dfs.values()])['year_decimal']
    else:
        if isinstance(plot_dfs, dict):
            x_vals = pd.concat(plot_dfs.values())['year_decimal']
        else:
            x_vals = plot_dfs['year_decimal']
    tick_year_decimals, tick_labels = _yeardecimal_date_tick_lebals(x_vals, tick_freq=x_tick_freq)
    fig, ax = plt.subplots(1, 1, layout='constrained', figsize=(7.5, 5))
    if isinstance(plot_dfs, dict):
        for i, (strain, plot_df) in enumerate(plot_dfs.items()):
            if style.endswith('with mid-points'):
                mid_point_df = plot_df[1]
                plot_df = plot_df[0]
                ax.scatter(mid_point_df['year_decimal'], mid_point_df['upper'], marker='v', c=palette[i])
                ax.scatter(mid_point_df['year_decimal'], mid_point_df['median'],marker='D', c=palette[i])
                ax.scatter(mid_point_df['year_decimal'], mid_point_df['lower'], marker="^", c=palette[i])

            ax.plot(plot_df['year_decimal'], plot_df['median'], color=palette[i], label=strain)
            ax.fill_between(plot_df['year_decimal'], plot_df['lower'], plot_df['upper'], color=palette[i],
                            alpha=0.2)
    else:
        i = 0
        plot_df = plot_dfs
        if style.endswith('with mid-points'):
            mid_point_df = plot_df[1]
            plot_df = plot_df[0]
            ax.scatter(mid_point_df['year_decimal'], mid_point_df['upper'], marker='v', c=palette[i])
            ax.scatter(mid_point_df['year_decimal'], mid_point_df['median'], marker='D', c=palette[i])
            ax.scatter(mid_point_df['year_decimal'], mid_point_df['lower'], marker="^", c=palette[i])

        ax.plot(plot_df['year_decimal'], plot_df['median'], color=palette[i])
        ax.fill_between(plot_df['year_decimal'], plot_df['lower'], plot_df['upper'], color=palette[i],
                        alpha=0.2)

    ax.xaxis.set_ticks(tick_year_decimals)
    ax.set_xticklabels(tick_labels)
    ax.tick_params(axis='x', labelrotation=45)
    ax.set_xlabel('Date')
    ax.set_ylabel(ylabel)
    if len(strains) > 1:
        plt.legend()
    plt.tight_layout()
    if include_grid:
        plt.grid(linestyle='-.')

    return fig, ax


def _yeardecimal_date_tick_lebals(series, tick_freq='monthly'):
    max_date = decimal_to_date(series.max())
    min_date = decimal_to_date(series.min())
    if tick_freq == 'yearly':
        freq = 'YS'
        strftime = '%Y'
    elif tick_freq == 'quarterly':
        freq = 'QS'
        strftime = '%Y-%b'
    elif tick_freq == 'monthly':
        freq = 'MS'
        strftime = '%Y-%b'
    elif tick_freq == 'half month':
        freq = 'SME'
        strftime = '%Y-%b-%d'
    elif tick_freq == 'weekly':
        freq = 'W'
        strftime = '%Y-%b-%d'
    else:
        raise ValueError('x_tick_freq must be "yearly", "quarterly", "monthly", ''half month" or "weekly".' +
                         '\n x_tick_freq=' + str(tick_freq))

    tick_dates = pd.date_range(start=min_date.replace(day=1), end=max_date, freq=freq)
    tick_labels = tick_dates.to_series().dt.strftime(strftime).to_list()
    tick_year_decimals = tick_dates.map(date_to_decimal)
    return tick_year_decimals, tick_labels

def plot_hist_kde(trace_df,
                  parameter,
                  xlabel=None,
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
    xlabel: str, default None
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
        Matplotlib figure, axis and dictionay of 'lower' and 'upper' hdi with median.
    """

    if xlabel is None:
        xlabel = parameter
    fig, ax = plt.subplots()  # initializes figure and plots
    sns.histplot(trace_df[parameter], ax=ax, kde=True, color=color)
    if hdi_prob is not None:
        hdi_est = hdi(trace_df[parameter].to_numpy(), hdi_prob=hdi_prob)
        upper_key = f'Upper {str(hdi_prob)} HDI'
        lower_key = f'Lower {str(hdi_prob)} HDI'
        hdi_est = {lower_key: hdi_est[0],
                   'Median': trace_df[parameter].median(),
                   upper_key: hdi_est[1]}
        ax.axvline(hdi_est['Median'], color='k', lw=2)
        ax.axvline(hdi_est[lower_key], color='k', ls='--', lw=1)
        ax.axvline(hdi_est[upper_key], color='k', ls='--', lw=1)
    ax.set_xlabel(xlabel)
    if tight_layout:
        plt.tight_layout()

    if hdi_prob is None:
        return fig, ax
    else:
        return fig, ax, hdi_est


def plot_origin_or_tmrca(trace_df, parameter, x_tick_freq='monthly', hdi_prob=None, color=None):
    """Plot histograms of origin or TMRCA, with kde line.

    Parameters
    ----------
    trace_df: pandas.DataFrame
        DataFrame with column headed with parameters.
    parameter: str
        Parameter name.
    x_tick_freq : str, default 'monthly'
        Frequency of x ticks.
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
        Matplotlib figure, axis and dictionay of 'lower' and 'upper' hdi with median.
    """
    if parameter not in ['Origin', 'TMRCA']:
        raise ValueError('parameters must be "Origin" or "TMRCA".')
    tick_year_decimals, tick_labels = _yeardecimal_date_tick_lebals(trace_df[parameter], tick_freq=x_tick_freq)
    output = plot_hist_kde(trace_df,
                            parameter,
                            xlabel=f'{parameter} (date)',
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
    plt.tight_layout()
    if hdi_prob is None:
        return fig, ax
    else:
        hdi_est = {key: decimal_to_date(value) for key, value in hdi_est.items()}
        return fig, ax, hdi_est


def plot_comparative_origin(df_melted, tick_freq='monthly', one_figure=False, palette=None):
    """Plot histograms of origins for each strain, with kde line..

    Parameters
    ----------
    df_melted : pandas.DataFrame
        DataFrame melted for use in Seaborne. Contains columns 'xml_set', 'variable' and
        'Estimate'.
    one_figure: bool, default False
        If True, then one figure will be created with the histograms for each strain
        overlayed on top of each other.
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
    df = df_melted[df_melted.variable == 'Origin']
    df.rename(columns={'Estimate': 'Origin'}, inplace=True)
    tick_year_decimals, tick_labels = _yeardecimal_date_tick_lebals(df.Origin, tick_freq=tick_freq)
    if one_figure:
        ax = sns.histplot(data=df, x='Origin', hue='xml_set', kde=True, palette=palette)
        if palette is None:
            palette = sns.color_palette()
        for colour, strain in zip(palette, df.xml_set.unique()):
            strain_df = df[df.xml_set == strain]
            ax.axvline(strain_df['Origin'].median(), color=colour)
            ax.axvline(strain_df['Origin'].quantile(0.05), color=colour, ls='--', lw=1)
            ax.axvline(strain_df['Origin'].quantile(0.95), color=colour, ls='--', lw=1)
        fig = ax
    else:
        fig = sns.FacetGrid(df, row="xml_set", hue='xml_set', margin_titles=True, aspect=4)
        fig.map_dataframe(sns.histplot, x='Origin', kde=True)
        fig.set_titles(row_template='')
        fig.add_legend()
        for ax, strain in zip(fig.axes.flat, df.xml_set.unique()):
            strain_df = df[df.xml_set == strain]
            ax.axvline(strain_df['Origin'].median(), color='k', lw=2)
            ax.axvline(strain_df['Origin'].quantile(0.05), color='k', ls='--', lw=1)
            ax.axvline(strain_df['Origin'].quantile(0.95), color='k', ls='--', lw=1)

    ax.xaxis.set_ticks(tick_year_decimals)
    ax.set_xticklabels(tick_labels)
    ax.tick_params(axis='x', labelrotation=45)
    plt.tight_layout()
    return fig
