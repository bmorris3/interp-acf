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
identify the dominant period in the smoothed autocorrelation function (after 
smoothing according to the recommendations of 
[McQuillan, Aigrain & Mazeh (2013)](http://adsabs.harvard.edu/abs/2013MNRAS.432.1203M)): 
```python
from interpacf import interpolated_acf, dominant_period

lag, acf = interpolated_acf(times_incomplete, fluxes_incomplete)
period = dominant_period(lag, acf, plot=True)
```
![](http://staff.washington.edu/bmmorris/images/acf.png)

#### Dependencies
* python 2/3
* numpy
* scipy
* matplotlib (optional)
