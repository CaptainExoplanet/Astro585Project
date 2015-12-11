# Astro585Project

This repository contains all the files neccessary to run an intermediate version of the adaptive anealed importance sampling algorithm described in Liu (2014). This version of the algorithm includes only the expectation and maximaization steps (leaving out deletion, addition, and merging).

Serial and parallel versions of the algorithm are provided in serial_version.jl and parallel_version.jl, respectively. The parallel version is set up to run on a single multi-core work station or a distributed cluster.

Testing functions are provided for serial and parallel code in testing_utilities.jl and testing_utilities_parallel.jl, respectively. For a simple, timed test, run run_serial_tests.jl or run_parallel_tests.jl.

Earlier versions of the code are provided in the EarlyAttempts subdirectory and can be ignored.
