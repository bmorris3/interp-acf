import numpy as np
import pytest
from ..interpacf import interpolated_acf, dominant_period

@pytest.mark.parametrize("n_points_missing",
                         [(i) for i in np.linspace(100, 300, 10, dtype=int)])
def test_periods(n_points_missing):
    # Make flux time-series with random noise, and
    # two periodic oscillations, one 70% the amplitude
    # of the other:
    np.random.seed(42)
    n_points = 1000
    primary_period = 2.5*np.pi
    secondary_period = 1.3*np.pi
    all_times = np.linspace(0, 6*np.pi, n_points)
    all_fluxes = 10 + (0.1*np.random.randn(len(all_times)) +
                       np.sin(2*np.pi/primary_period * all_times) +
                       0.7*np.cos(2*np.pi/secondary_period * (all_times - 2.5)))

    # Remove some fluxes, times from those data:
    #n_points_missing = 200  # This number is approximate
    missing_indices = np.unique(np.random.randint(0, n_points,
                                                  size=n_points_missing))
    mask = list(set(np.arange(len(all_times))).difference(set(missing_indices)))
    times_incomplete = all_times[mask]
    fluxes_incomplete = all_fluxes[mask]

    # Need zero-mean fluxes:
    fluxes_incomplete -= np.mean(fluxes_incomplete)

    # Compute autocorrelation function
    lag, acf = interpolated_acf(times_incomplete, fluxes_incomplete)

    # Find dominant period in autocorrelation function
    detected_period = dominant_period(lag, acf)

    assert ((primary_period - detected_period)/primary_period) < 0.02