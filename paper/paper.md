---
title: 'interpacf: compute autocorrelation functions in light curves with missing data.'
tags:
  - astronomy
  - python
authors:
 - name: Brett M. Morris
   orcid: 0000-0003-2528-3409
   affiliation: 1
 - name: John C. Lurie
   orcid: 0000-0002-8114-0835
   affiliation: 1
affiliations:
 - name: University of Washington, Seattle, WA
   index: 1
date: 28 January 2018
bibliography: paper.bib
---

# Summary

This is a simple package for measuring autcorrelation functions from light 
curves to search for periodicities. The major feature included in `interpacf` 
which is not included in numpy's autocorrelation function is the ability to
interpolate over missing data, as is commonly the case in astronomical light 
curves. 

# References

