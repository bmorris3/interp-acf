"""
Calculate the autocorrelation function for an evenly sampled time-series with
missing data.
"""
from __future__ import print_function, absolute_import, division

from scipy.ndimage import gaussian_filter
from scipy import signal
import numpy as np

__all__ = ["interpolated_acf", "autocorrelation", "interpolate_missing_data",
           "dominant_period"]
__version__ = "0.1"
__author__ = "Brett Morris (bmmorris@uw.edu)"


def interpolate_missing_data(times, fluxes):
    """
    Assuming ``times`` are uniformly spaced with missing cadences,
    fill in the missing cadences with linear interpolation.

    Parameters
    ----------
    times : numpy.ndarray
        Incomplete but otherwise uniformly sampled times
    fluxes : numpy.ndarray
        Flux for each time in ``times``

    Returns
    -------
    interpolated_times : numpy.ndarray
        ``times`` with filled-in missing cadences
    interpolated_fluxes : numpy.ndarray
        ``fluxes`` with filled-in missing cadences
    """
    # Find typical time between cadences:
    dt = np.median(np.diff(times))
    first_time = times[0]

    # Approximate the patchy grid of integer cadence indices,
    # i.e.: (0, 1, 3, 4, 5, 8, ...)
    cadence_indices = np.rint((times - first_time)/dt)
    # Find missing cadence indices if that grid were complete
    expected_cadence_indices = set(np.arange(cadence_indices.min(),
                                             cadence_indices.max()))
    missing_cadence_indices = expected_cadence_indices.difference(set(cadence_indices))
    # Convert the missing cadences to times
    missing_times = first_time + np.array(list(missing_cadence_indices))*dt

    # Interpolate to find fluxes at missing times
    interp_fluxes = np.interp(missing_times, times, fluxes)

    # Combine the interpolated and input times, fluxes
    interpolated_fluxes = np.concatenate([fluxes, interp_fluxes])
    interpolated_times = np.concatenate([times, missing_times])

    # Sort the times, fluxes, so that you can compute the ACF on them:
    sort_by_time = np.argsort(interpolated_times)
    interpolated_fluxes = interpolated_fluxes[sort_by_time]
    interpolated_times = interpolated_times[sort_by_time]
    return interpolated_times, interpolated_fluxes


def autocorrelation(x):
    """
    Calculate the autocorrelation function of array ``x``.
    """
    result = np.correlate(x, x, mode='full')
    return result[result.size//2:]


def interpolated_acf(times, fluxes):
    """
    Calculate the autocorrelation function after interpolating over
    missing times and fluxes.

    Parameters
    ----------
    times : numpy.ndarray
        Incomplete but otherwise uniformly sampled times
    fluxes : numpy.ndarray
        Flux for each time in ``times``

    Return
    ------
    lag : numpy.ndarray
        Time lag
    acf : numpy.ndarray
        Autocorrelation function
    """
    if not np.all(np.sort(times) == times):
        raise ValueError("Arrays must be in chronological order to compute ACF")

    # Interpolate over missing times, fluxes
    interpolated_times, interpolated_fluxes = interpolate_missing_data(times,
                                                                       fluxes)
    # Calculate the grid of "lags" in units of ``times``
    dt = np.median(np.diff(interpolated_times))
    lag = dt*np.arange(len(interpolated_fluxes))

    # Compute the autocorrelation function on interpolated fluxes
    acf = autocorrelation(interpolated_fluxes)
    return lag, acf


def dominant_period(lag, acf, min=None, max=None, plot=False):
    """
    Find the dominant period in the autocorrelation function.

    Parameters
    ----------
    lag : numpy.ndarray
        Time lags
    acf : numpy.ndarray
        Autocorrelation function
    min : float (optional)
        Return dominant period greater than or equal to ``min``
    max : float (optional)
        Return dominant period less than or equal to ``max``
    plot : bool (optional)
        Plot the autocorrelation function, peak detected

    Return
    ------
    acf_period : float
        Dominant period detected via the autocorrelation function
    """
    lag_limited = np.copy(lag)
    acf_limited = np.copy(acf)

    # Apply limits if any are input:
    if min is not None and max is None:
        acf_limited = acf_limited[lag_limited > min]
        lag_limited = lag_limited[lag_limited > min]
    elif max is not None and min is None:
        acf_limited = acf_limited[lag_limited < max]
        lag_limited = lag_limited[lag_limited < max]
    elif max is not None and min is not None:
        acf_limited = acf_limited[(lag_limited < max) & (lag_limited > min)]
        lag_limited = lag_limited[(lag_limited < max) & (lag_limited > min)]

    smooth_acf = gaussian_filter(acf_limited, 10)
    relative_maxes = signal.argrelmax(smooth_acf)[0]
    if len(relative_maxes) == 0:
        raise ValueError("No period found. Be sure to check limits and to "
                         "median-subtract your fluxes.")

    absolute_max_index = relative_maxes[np.argmax(smooth_acf[relative_maxes])]
    acf_period = lag_limited[absolute_max_index]

    if plot:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(lag, acf/np.max(acf), label='ACF')
        plt.plot(lag, smooth_acf/np.max(smooth_acf), label='Smoothed ACF')
        plt.axvline(acf_period, ls='--', color='r', label='Primary period')
        plt.legend()
        plt.title('ACF')
        plt.xlabel('Lag')

    return acf_period
