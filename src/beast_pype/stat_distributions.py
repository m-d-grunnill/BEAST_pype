"""Functions for creating dataframes of gamma distribution and making gama distrbution related calcualations."""

from scipy.stats import gamma
import seaborn as sns
import numpy as np
import pandas as pd


## Gamma

def gamma_from_mean(mu, sd=None, variance=None):
    """Convert mean and either sd or variance to shape (a) and scale of gama distribution.

    See:
        Bolker, Benjamin M. 2008. “Gamma.” In Ecological Models and Data in R, 131–133. Princeton University Press.
        https://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation
        https://stats.stackexchange.com/questions/342639/how-to-find-alpha-and-beta-from-a-gamma-distribution


    Parameters
    ----------
    mu : int or float
        Mean
    sd : int or float, optional
        Standard Deviation, by default None. If variance is None value is assumed to be 1.
    variance : int or float, optional
        Variance, by default None. If sd is None value is assumed to be 1.

    Returns
    -------
    dict {'a' : float, 'scale': float}
       Diction containing shape (a) and scale of gama distribution.
    """
    if variance is None:
        if sd is None:
            variance = 1
        else:
            variance = sd ** 2
    return {'a': (mu ** 2) / variance, 'scale': variance / mu}


def mean_rate_years_gamma_df(mu, sd=None, variance=None, loc=0, min_range=0.1, max_range=500, step=0.1,  return_gamma_dict=False):
    """Generate dataframe outlining pdf of gamma distribution for mean yearly rate.

    Parameters
    ----------
    mu : float or int
        Mean in years
    sd : float or int, optional
        Standard deviation in years., by default None
    variance : float or int, optional
        Variance in years, by default None
    loc : float or int, optional
        Location, by default 0.
    min_range : int, optional
        Min range value in years, by default 0.1
    max_range : int, optional
        Max range value in years, by default 100
    step : int, optional
        Step in years, by default 0.1


    Returns
    -------
    gamma_dict: dict {'a' : float, 'scale': float}
        Diction containing shape (a) and scale of gama distribution.
    traces_df: pandas.DataFrame
        Dataframe outlining pdf of gamma distribution for mean yearly rate.
    """

    gamma_dict = gamma_from_mean(mu=mu, sd=sd, variance=variance)
    gamma_dict['loc'] = loc
    df = gamma_df(**gamma_dict, min_range=min_range, max_range=max_range, step=step)
    if return_gamma_dict:
        return gamma_dict, df
    else:
        return df

def gamma_df(a, loc=0, scale=1, min_range=0.1, max_range=500, step=0.1):
    """Generate dataframe outlining pdf of gamma distribution showing yearly rate converted to daily rate and day period.


    Parameters
    ----------
    a : float or int
        Shape
    loc : float or int, optional
        Location, by default 0.
    scale : float
        Scale, by default 1.
    min_range : int, optional
        Min range value in years, by default 0.1
    max_range : int, optional
        Max range value in years, by default 100
    step : int, optional
        Step in years, by default 0.1

    Returns
    -------
    pandas.DataFrame


    """
    rate_years = np.arange(min_range, max_range + step, step=step)
    rate_days = rate_years / 365
    period_days = 1 / rate_days
    probabilities = gamma.pdf(x=rate_years, a=a, loc=loc, scale=scale)
    return pd.DataFrame({'Rate per Year': rate_years, 'Rate per Day':  rate_days,
                         'Period in Days': period_days,
                         'Probability': probabilities})

def plot_gamma_df(gamma_df):
    """Plot dataframe resulting from gamma_df function.

    Parameters
    ----------
    gamma_df : pandas.DataFrame

    Returns
    -------
    seaborn.lineplot

    """
    df_melt = pd.melt(gamma_df, id_vars=['Probability'])
    df_melt_mod = df_melt[df_melt['Probability'] > 0.0001]
    fig = sns.relplot(data=df_melt_mod, x="value", y="Probability",
                      col="variable", col_order=['Period in Days', 'Rate per Day', 'Rate per Year'],
                      kind='line',
                      facet_kws={'sharex': False})
    fig.set_titles(col_template="")
    fig.axes[0, 0].set_xlabel('Infectious Period in Days')
    fig.axes[0, 1].set_xlabel('Rate of Becoming Uninfectious per Day')
    return fig

## LogNormal


## Beta


