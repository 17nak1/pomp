## pomp:
An R package for statistical inference on partially observed Markov processes

#mif2:
Maximum likelihood by iterated filtering

An iterated filtering algorithm for estimating the parameters of a partially-observed Markov process.
Running `mif2` causes the algorithm to perform a specified number of particle-filter iterations.
At each iteration, the particle filter is performed on a perturbed version of the model, in which the parameters to be estimated are subjected to random perturbations at each observation.
This extra variability effectively smooths the likelihood surface and combats particle depletion by introducing diversity into particle population.
