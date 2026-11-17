# DISQ 0.2.1

- Organize the fixed-projection implementation into local roots, global
  projection, scores, communication, and spatial inference.
- Add a simulated example from local SQBiT fitting through communication.
- Revise the README and function descriptions around the analysis workflow.
- Keep the existing function arguments, estimator defaults, and numerical
  calculations unchanged.

# DISQ 0.2.0

- Add `disq_from_local_fits()` for decorrelated initialization, global
  smoothness correction, fixed-projection communication, and spatial covariance.
- Add `disq_bandwidths()` with the manuscript cube-root rules and default
  common multiplier 0.20.
- Retain the score and coefficient traces for all requested communication rounds.
