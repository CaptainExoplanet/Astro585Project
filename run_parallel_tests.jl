#This file tests and benchmarks the parallel version of the algorthm at varying numbers of iterations and samples.

#addprocs(3);

include("testing_utilities_parallel.jl");

#run once with small N to compile
@time test_run_algorithm_2D(10,100000,2.);
#test_run_algorithm_3D(10,100,2.);

### Test algorith at varying N in 2 dimensions ###
#=
println("Running test in 2 dimensions with 10 iterations and 100 samples");
@time test_run_algorithm_2D(10,100,2.);

println("Running test in 2 dimensions with 10 iterations and 1000 samples");
@time test_run_algorithm_2D(10,1000,2.);

println("Running test in 2 dimensions with 10 iterations and 10000 samples");
@time test_run_algorithm_2D(10,10000,2.);

println("Running test in 2 dimensions with 10 iterations and 100000 samples");
@time test_run_algorithm_2D(10,100000,2.);

println("Running test in 2 dimensions with 100 iterations and 1000 samples");
@time test_run_algorithm_2D(100,1000,1.);

println("Running test in 2 dimensions with 1000 iterations and 1000 samples");
@time test_run_algorithm_2D(1000,1000,1.);


### Test algorithm at varying N in 3 dimentions ###

println("Running test in 3 dimensions with 10 iterations and 100 samples");
@time test_run_algorithm_3D(10,100000,2.);

println("Running test in 3 dimensions with 10 iterations and 1000 samples");
@time test_run_algorithm_3D(10,1000,2.);

println("Running test in 3 dimensions with 10 iterations and 10000 samples");
@time test_run_algorithm_3D(10,10000,2.);

println("Running test in 3 dimensions with 10 iterations and 100000 samples");
@time test_run_algorithm_3D(10,100000,1.);

println("Running test in 3 dimensions with 100 iterations and 1000 samples");
@time test_run_algorithm_3D(100,1000,1.);

println("Running test in 3 dimensions with 1000 iterations and 1000 samples");
@time test_run_algorithm_3D(1000,1000,1.);

println("All tests have been passed!");
=#
