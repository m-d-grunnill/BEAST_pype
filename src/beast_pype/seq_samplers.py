import pandas as pd
from random import sample


def split_df_by_weeks(df, date_field='date'):
    """ 
    Splits a dataframe into chunks of a week based on Monday being the first day of the week.
    """
    return [g for n, g in df.groupby(pd.Grouper(key=date_field, freq='W'))]
    
    
def df_weekly_sampler(df, n, draws=1, date_field='date'):
    """
    Draw a random sample from a dataframe spread evenly accross weekly chunks.
    
    If a week has fewer rows than the sample needed then the remainder is spread
     across all other rows.
    Any remainder after this is distributed by drawing an extra row from other
    weeks selected at random.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Pandas DataFrame from which the random samples are drawn.
    n : int
        Total number of samles to draw.
    draws : int
        Number of draws of size N you wish to make. Replacement occurs between draws.
    date_field : str, default 'date'
        Field of date variable to base chopping up selection into weeks.

    Returns
    -------
    draws_data : pandas.DataFrame or Mast of pandas.DataFrames
        If draws > 1 a list of pandas.DataFrames is returned, otherwise a single dataframe.
    
    """
    if len(df) < n:
        raise ValueError('Dataframe has less rows than the sample size selected.')
    if len(df) == n:
        raise ValueError('Dataframe has the same number of rows as the sample size selected.')

    week_split_dfs = split_df_by_weeks(df, date_field=date_field)
    number_of_weeks = len(week_split_dfs)
    modulus = n % number_of_weeks
    typical_weekly_sample = int(n/number_of_weeks)

    available_per_week = [len(week_df) for week_df in week_split_dfs]

    less_than_equal_typical_sample = [available <= typical_weekly_sample for available in available_per_week]
    less_than_needed = typical_weekly_sample - pd.Series(available_per_week)[less_than_equal_typical_sample]
    extras_needed = less_than_needed.sum()+modulus

    if extras_needed > 0:
        n_weeks_more_than_needed = less_than_equal_typical_sample.count(False)
        extras_needed_per_available = int(extras_needed / n_weeks_more_than_needed)
        typical_weekly_sample += extras_needed_per_available
        less_than_equal_typical_sample = [available <= typical_weekly_sample for available in available_per_week]
        index_greater_than_typical_sample = [
            week_i for week_i, available in
            enumerate(available_per_week)
            if available > typical_weekly_sample
        ]
        n_weeks_more_than_needed = len(index_greater_than_typical_sample)
        final_modulus = extras_needed % n_weeks_more_than_needed
        weeks_to_add_one_extra = sample(index_greater_than_typical_sample, final_modulus)
    else:
        weeks_to_add_one_extra = []

    draws_dfs = []
    for draw in range(draws):
        sampled_data = []
        for week_i, weeK_df in enumerate(week_split_dfs):
            if less_than_equal_typical_sample[week_i]:
                sampled_data.append(weeK_df)
            else:
                this_weeks_sample = typical_weekly_sample + int(week_i in weeks_to_add_one_extra)
                sampled_data.append(weeK_df.sample(this_weeks_sample))
        sampled_data = pd.concat(sampled_data)
        draws_dfs.append(sampled_data)
    
    if draw == 0:
        return draws_dfs[0]
    else:
        return draws_dfs