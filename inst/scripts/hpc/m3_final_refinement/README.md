# Final refinement of the M3 maximum likelihood fit

These scripts completed the fit of M3 (`taildown_covariates.R`) after the
equal-step finite-difference Hessian of the package refinement was found to be
wrong (negative eigenvalue, zero standard errors): the log-likelihood has a long
flat ridge, and its curvature differs by orders of magnitude between directions.
They were run in this order, with the paths of the original analysis:

1. `hess2_M3.R`: Hessian at the first estimate along the eigenvectors of the
   equal-step Hessian (steps of half a standard error); showed a predicted gain of 1.1.
2. `newton_M3.R`: Newton steps with whitened Hessians (two iterations, -1193.25 to -1192.96).
3. `bfgs_M3.R`: BFGS in whitened coordinates (to -1192.31, last steps gaining 0.005).
4. `final_M3.R`: whitened Hessian and Wald intervals at the estimate, cross-validation
   and grid predictions; writes `fit_M3.rds` and `intervals_M3.csv`.

The same whitened Hessian is now used by the refinement in `odisp_fit_zil_ml`.
