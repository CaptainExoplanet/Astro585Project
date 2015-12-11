#This file tests and benchmarks the parallel version of the algorthm at varying numbers of iterations and samples.

#addprocs(3);

include("testing_utilities_parallel.jl");

#run once with small N to compile
test_run_algorithm_2D(10,100,2.);
test_run_algorithm_3D(10,100,2.);

println("Running test in 2 dimensions with 10 iterations and 100000 samples");
@time test_run_algorithm_2D(10,100000,2.);

println("Running test in 3 dimensions with 10 iterations and 100000 samples");
@time test_run_algorithm_3D(10,100000,2.);

println(“all tests pass”);


