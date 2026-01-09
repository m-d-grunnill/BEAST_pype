from copy import deepcopy
from pandas.tseries.offsets import DateOffset

def gen_partition_dates(partition_freq, date_of_youngest_tip, voi_youngest_oldest_tip, voi_strains, dr_strain):
    """

    Parameters
    ----------
    partition_freq
    date_of_youngest_tip
    voi_youngest_oldest_tip
    voi_strains
    dr_strain

    Returns
    -------

    """
    voi_partition_dates = []
    offsets = {partition_freq['unit']: partition_freq['every']}
    date_to_append = date_of_youngest_tip - DateOffset(**offsets)
    if not date_to_append > voi_youngest_oldest_tip:
        raise ValueError(
            f'Change dates dictionary for partition_freq makes no sense.\n' +
            f"First partition date generated is not after the most recent tip date out of oldest tips for each VOI lineage  ({voi_youngest_oldest_tip.strftime('%Y-%m-%d')}).\n" +
            'Suggest using smaller "units" or "every" values or pushing "end" date back.'
        )
    while date_to_append > voi_youngest_oldest_tip:
        voi_partition_dates.append(date_to_append.strftime('%Y-%m-%d'))
        date_to_append = date_to_append - DateOffset(**offsets)

    dr_partition_dates = deepcopy(voi_partition_dates)
    dr_partition_dates.append(date_to_append.strftime('%Y-%m-%d'))
    partition_dates = {f"VOI_{voi_strain}": voi_partition_dates for voi_strain in voi_strains}
    partition_dates[f'DR_{dr_strain}'] = dr_partition_dates
    return partition_dates