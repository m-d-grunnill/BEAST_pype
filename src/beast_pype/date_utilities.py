from datetime import timedelta, datetime
from datetime import datetime as dt
import calendar
from dateutil.relativedelta import relativedelta, SU, MO
import time


def _since_epoch(date):  #
    """
    Returns seconds since epoch

    Parameters
    ----------
    date: datetime.date or datetime.datetime
        Date to be converted.

    Returns
    -------
    float
    """
    return time.mktime(date.timetuple())


def decimal_to_date(decimal):
    """
    Converts a decimal year to a date object.

    Parameters
    ----------
    decimal: float
        The decimal value to convert.

    Returns
    -------
    date: datetime.date
    """
    ''' Convert year decimal to date.'''
    year = int(decimal)
    d = timedelta(days=(decimal - year) * (365 + calendar.isleap(year)))
    day_one = datetime(year, 1, 1)
    date = d + day_one
    return date


def date_to_decimal(date):
    """ Convert date to year decimal (year fraction).

    Parameters
    ----------
    date: datetime.date or datetime.datetime
        Date to be converted.

    Returns
    -------
    date_as_year_decimal : float
    """
    year = date.year
    start_of_this_year = dt(year=year, month=1, day=1)
    start_of_next_year = dt(year=year + 1, month=1, day=1)

    year_elapsed = _since_epoch(date) - _since_epoch(start_of_this_year)
    year_duration = _since_epoch(start_of_next_year) - _since_epoch(start_of_this_year)
    fraction = year_elapsed / year_duration

    return date.year + fraction


def the_next_sunday(date, dont_if_sunday=False):
    """
    Find the next Sunday after a given date.

    Parameters
    ----------
    date: datetime.date or datetime.datetime
    dont_if_sunday: bool, default False
        If already a Sunday, return date.

    Returns
    -------
    datetime.datetime: The next Sunday after date.
    """
    if dont_if_sunday and date.weekday()==6:
        return date
    else:
        return date + relativedelta(weekday=SU(1))

def the_previous_monday(date, dont_if_monday=False):
    """
    Find the previous Monday before a given date.

    Parameters
    ----------
    date: datetime.date or datetime.datetime
    dont_if_monday: bool, default False
        If already a Monday, return date.

    Returns
    -------
    datetime.datetime: The previous Monday before date.
    """
    if dont_if_monday and date.weekday()==0:
        return date
    else:
        return date + relativedelta(weekday=MO(-1))