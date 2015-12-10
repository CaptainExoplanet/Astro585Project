cd("$(homedir())/Astro585Project");

include("testing_utilities.jl");
#include("testing_utilities_parallel.jl");

# Test serial accuracy as N increases

test_accuracy_2D(10,100);
test_accuracy_2D(100,100);
test_accuracy_2D(1000,100);
test_accuracy_2D(10000,100);

test_accuracy_2D(10,1000);
test_accuracy_2D(10,10000);

test_accuracy_2D(100,1000);
test_accuracy_2D(1000,1000);