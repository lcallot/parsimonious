parsimonious
=====

A R package for the estimation of models with parsimoniously time varying-parameters, by Laurent Callot and Johannes Tang Kristensen. 

[![Build Status](https://travis-ci.org/lcallot/pcvar.png?branch=master)](https://travis-ci.org/lcallot/parsimonious)



Content
-----------

This package permits the estimation of parsimoniously time-varying parameter models as in *Vector Autoregressions with Parsimoniously Time-Varying Parameters and an Application to Monetary Policy.*

 * [Link to the paper](http://lcallot.github.io/papers/ptv-var/)
 * [Link to the replication repo](https://github.com/lcallot/ptv-var)
 
 
Usage
----------

The function `ptvfit` estimates the model given by a formula and returns a _parsimonious_ object for which a `plot` and a `print` method are available.

