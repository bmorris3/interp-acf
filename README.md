# interp-acf
Calculate the autocorrelation function for an evenly sampled time-series with 
missing data.

### Install
Clone this repository, change directories into it, and run this command to 
install the package: 
```
python setup.py install
```

### Example
**Check out the a longer example with plots 
[here](https://github.com/bmorris3/interp-acf/blob/master/example.ipynb)**. 

For time series (`times_incomplete`, `fluxes_incomplete`) with missing data, 
compute the autocorrelation function after interpolating over missing data, and
identify the dominant period in that autocorrelation function: 
```
from interpacf import interpolated_acf, dominant_period
lag, acf = interpolated_acf(times_incomplete, fluxes_incomplete)
period = dominant_period(lag, acf)
```
