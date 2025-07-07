from src.beast_pype.outputs import (read_strain_logs_for_plotting, plot_box_violin,
                                    plot_comparative_skyline, plot_comparative_origin, percentile_pivot)
from src.beast_pype.outputs import _set_dates_skyline_plotting_df
from src.beast_pype.date_utilities import decimal_to_date
import scipy
import os
import pandas as pd
from datetime import datetime
import json
import warnings

# stop annoying matplotlib warnings
warnings.filterwarnings("ignore", module="matplotlib\*")


# def is_date(string):
#     """
#     Return whether the string can be interpreted as a date.
#
#     :param string: str, string to check for date
#     :param fuzzy: bool, ignore unknown tokens in string if True
#
#     Source
#     ---------
#     https://stackoverflow.com/questions/25341945/check-if-string-has-date-any-format
#
#    """
#     try:
#         parser(string)
#         return True
#
#     except ValueError:
#         return False

def _last_period_rt_ratio(rt_type, r_t_plot_dfs_dict, chosen_samples):
    ratio_of_rts = []
    for sample in chosen_samples:
        denominator = r_t_plot_dfs_dict[sample][rt_type].iloc[-1]
        ratio_of_rts.append({column:r_t_plot_dfs_dict[column][rt_type].iloc[-1]/denominator for column in chosen_samples})

    ratio_of_rts = pd.DataFrame.from_records(ratio_of_rts)
    ratio_of_rts.index = chosen_samples
    return ratio_of_rts

class BP2CW:
    """
    Attributes
    ---------------
    latest_run_dir: str
        Path to the directory containing the results of the latest run of the beast_pype .
    date_of_run: str
        Date the latest run of pipeline was started, formatted as YYYY-MM-DD.
    sample_info: Pandas DataFrame
        Information on the samples used in the run.
    strain_sub_varients: {str: [str]}
        Subvariants of variants listed in sample_info.
    infection_period_box_violin_with_prior: matplotlib.axes._axes.Axes
        Box violin plot of infection period with prior.
    infection_period_box_violin: matplotlib.axes._axes.Axes
        Box violin plot of infection period.
    infection_period_percentile_df: Pandas DataFrame
        Percentiles of infection period.
    clock_rate_box_violin_with_prior: matplotlib.axes._axes.Axes
        Box violin plot of clock rate with prior.
    clock_rate_box_violin: matplotlib.axes._axes.Axes
        Box violin plot of clock rate.
    clock_rate_percentile_df: Pandas DataFrame
        Percentiles of clock rate.
    sampling_proportion_box_violin_with_prior: matplotlib.axes._axes.Axes
        Box violin plot of sampling proportion with prior.
    sampling_proportion_box_violin: matplotlib.axes._axes.Axes
        Box violin plot of sampling proportion.
    sampling_proportion_percentile_df: Pandas DataFrame
        Percentiles of sampling.
    rt_true_skyline_fig: matplotlib.figure.Figure
        Rt true skyline plot.
    rt_smoothed_fig: matplotlib.figure.Figure
        Rt skyline smoothed over srid mid-points plot.
    rt_smoothed_with_points_fig: matplotlib.figure.Figure
        Rt skyline smoothed over srid mid-points plot. Mid-points plotted.
    rt_df: Pandas DataFrame
        Rt information
    last_period_rt_ratio_lower_df: Pandas DataFrame
        Ratios of lower HPD Rt value in last time period.
    last_period_rt_ratio_median_df: Pandas DataFrame
        Ratios of median Rt value in last time period.
    last_period_rt_ratio_upper_df: Pandas DataFrame
        Ratios of upper HPD Rt value in last time period.
    orign_overlayed_fig: matplotlib.axes._axes.Axes
        Origin histogram/kde plots overlayed into one figure.
    orign_stacked_fig: seaborn.axisgrid.FacetGrid
        Origin histogram/kde plots stacked vertically (several plots).
    orign_percentile_df: Pandas DataFrame
        Percentiles of origin estimates.
    """

    def __init__(self, pipeline_runs_dir):
        """

        Parameters
        ---------------
        pipeline_runs_dir: str
            Directory where beast_pype pipeline runs are located. Each subdirectory
            must:
             * Be a run of the beast_pype pipeline.
             * A date of format YYYY-MM-DD.
        """
        self.pipeline_runs_dir = pipeline_runs_dir
        availleble_runs = [f'{pipeline_runs_dir}/{path}' for path in os.listdir(pipeline_runs_dir)
                           if os.path.isdir(f'{pipeline_runs_dir}/{path}')]
        self.latest_run_dir = max(availleble_runs, key=os.path.getctime)
        self.date_of_run = self.latest_run_dir.rsplit('/', 1)[-1]
        with open(self.latest_run_dir + "/pipeline_run_info.json", "r") as file:
            data = file.read()
        file.close()
        pipeline_run_info = json.loads(data)
        sample_dirs = pipeline_run_info["sample directories"]
        chosen_samples = list(sample_dirs.keys())

        voi_strain_sample_dirs = {}
        dr_strain_sample_dirs = {}
        for sample, path in sample_dirs.items():
            if sample.split('_')[0]:
                dr_strain_sample_dirs[sample] = path
            else:
                voi_strain_sample_dirs[sample] = path

        records = [
            {'Type of Strain': sample.split('_')[0], 'Strain': sample.split('_')[1], 'Sample Size': sample.split('_')[2]}
            for sample in dr_strain_sample_dirs.keys()]
        records += [
            {'Type of Strain': sample.split('_')[0], 'Strain': sample.split('_')[1], 'Sample Size': sample.split('_')[2]}
            for sample in voi_strain_sample_dirs.keys()]

        self.sample_info = pd.DataFrame.from_records(records)
        self.strain_sub_varients = pipeline_run_info["strain sub-varients"]
        log_path_dict = {sample: directory + '/merged.log' for sample, directory in sample_dirs.items()}
        youngest_tips_dict = {sample: datetime.strptime(date_str, '%Y-%m-%dir')
                              for sample, date_str in pipeline_run_info['youngest tip dates'].items()
                              if sample in chosen_samples}
        df, df_melted_for_seaborn = read_strain_logs_for_plotting(
            file_path_dict=log_path_dict,
            convert_become_uninfectious_rate=True,
            youngest_tips_dict=youngest_tips_dict)

        ## Infection period suff
        gamma_prior_params = {'a': 26.809699121992104, 'loc': 0, 'scale': 1.4548995665894278}
        prior = scipy.stats.gamma(**gamma_prior_params)
        yearly_rate_prior_draw = prior.rvs(size=int(1e5))
        daily_rate_prior_draw = yearly_rate_prior_draw / 365.25
        inf_period_draw = 1 / daily_rate_prior_draw
        self.infection_period_box_violin_with_prior = plot_box_violin(df_melted_for_seaborn, 'Infection period (per day)', prior_draws=inf_period_draw)
        self.infection_period_box_violin = plot_box_violin(df_melted_for_seaborn, 'Infection period (per day)')
        self.infection_period_percentile_df = percentile_pivot(df, 'Infection period (per day)')

        ### The evolutionary substitution rate (i.e. clock rate) is estimated. The prior given here is based on a strict clock estimate for SARS-CoV-2. Note clock rate's unit is 'number of nucleotide substitutions per site per year'.
        gamma_prior_params = {'a': 1.7777777777777781, 'scale': 0.00022499999999999994}
        prior = scipy.stats.gamma(**gamma_prior_params)
        subs_rate_prior_draw = prior.rvs(size=int(1e5))
        self.clock_rate_box_violin_with_prior = plot_box_violin(df_melted_for_seaborn, 'clockRate', prior_draws=subs_rate_prior_draw)
        self.clock_rate_box_violin = plot_box_violin(df_melted_for_seaborn, 'clockRate')
        self.clock_rate_percentile_df = percentile_pivot(df, 'clockRate')

        ## Sampling Proportion
        beta_prior_params = {'a': 1, 'b': 999}
        prior = scipy.stats.beta(**beta_prior_params)
        sampling_prior_draw = prior.rvs(size=int(1e5))
        self.sampling_proportion_box_violin_with_prior = plot_box_violin(df_melted_for_seaborn,
                                                                           'samplingProportion_BDSKY_Serial',
                                                                           prior_draws=sampling_prior_draw)
        self.sampling_proportion_box_violin = plot_box_violin(df_melted_for_seaborn, 'samplingProportion_BDSKY_Serial')
        self.sampling_proportion_percentile_df = percentile_pivot(df, 'samplingProportion_BDSKY_Serial')

        ## $R_T$
        type_change_date_dicts = {sample: [datetime.strptime(date_str, '%Y-%m-%dir') for date_str in date_list]
                                  for sample, date_list in pipeline_run_info['Re change dates'].items()}
        sample_change_date_dicts = {sample: type_change_date_dicts[sample.split('_')[0]] for sample in chosen_samples}
        self.rt_true_skyline_fig, axs = plot_comparative_skyline(df,
                                                                 sample_change_date_dicts,
                                                                 youngest_tips_dict)
        self.rt_smoothed_fig, axs = plot_comparative_skyline(df,
                                                             sample_change_date_dicts,
                                                             youngest_tips_dict,
                                                             style='smooth spline')
        self.rt_smoothed_with_points_fig, axs = plot_comparative_skyline(df,
                                                                         sample_change_date_dicts,
                                                                         youngest_tips_dict,
                                                                         style='smooth spline with mid-points')


        ### RT tables
        parameter_start = 'reproductiveNumber'
        style = 'skyline'
        rt_dfs_dict = {}
        for sample in df.Strain_and_Sample_Size.unique():
            temp_df = _set_dates_skyline_plotting_df(df[df['Strain_and_Sample_Size'] == sample],
                                                     parameter_start=parameter_start,
                                                     style=style,
                                                     youngest_tip=youngest_tips_dict[sample],
                                                     dates_of_change=sample_change_date_dicts[sample])

            temp_df['Strain_and_Sample_Size'] = sample
            start_df = temp_df[temp_df.index % 2 == 0]
            end_df = temp_df[temp_df.index % 2 != 0]
            start_df['Start of Period'] = start_df.year_decimal.map(decimal_to_date).map(
                lambda x: x.strftime('%Y-%m-%dir'))
            start_df['End of Period'] = end_df.year_decimal.map(decimal_to_date).map(
                lambda x: x.strftime('%Y-%m-%dir')).to_list()
            rt_dfs_dict[sample] = start_df

        rt_df = pd.concat(rt_dfs_dict.values())
        rt_df = rt_df[
            ['Strain_and_Sample_Size', 'Start of Period', 'End of Period', 'lower', 'median', 'upper']]
        self.rt_df = rt_df.set_index('Strain_and_Sample_Size')
        self.rt_df.index = self.rt_df.index.str.split('_', expand=True)
        self.rt_df.reset_index(inplace=True)
        self.rt_df.columns = ['Type of Strain', 'Strain', 'Sample Size', 'Start of Period', 'End of Period', 'lower',
                               'median', 'upper']

        ## Ratios of Rt in last time period
        self.last_period_rt_ratio_lower_df = _last_period_rt_ratio('lower', rt_dfs_dict, chosen_samples)
        # %% md
        # ### Median
        self.last_period_rt_ratio_median_df = _last_period_rt_ratio('median', rt_dfs_dict, chosen_samples)
        # %% md
        # Rows are denominators and columns are numerators.
        self.last_period_rt_ratio_upper_df = _last_period_rt_ratio('upper', rt_dfs_dict, chosen_samples)

        ## Origin
        self.orign_overlayed_fig = plot_comparative_origin(df_melted_for_seaborn, one_figure=True)
        self.orign_stacked_fig = plot_comparative_origin(df_melted_for_seaborn)
        self.orign_percentile_df = percentile_pivot(df, 'Origin')
