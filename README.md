## time-series

Code for time series analysis.

### `cmefilt.m` Common-mode error filter

Code for removing common-mode error in GPS time series. Common mode is error that is coherent across an entire network of GPS receivers.

### `fitTS.m` GPS time series trend estimation

Code for fitting horizontal GPS displacement time series for secular (linear + quadratic) trend, seasonal (annual + semi-annual), coseismic offsets, and postseismic transients. The time series model is given by `tsmodel.m`. Code uses MATLAB's built-in `lsqcurvefit` nonlinear least-squares solver.
