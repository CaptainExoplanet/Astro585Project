println("Hello World")

using Distributions
using Base.Test
include("$(homedir())/Astro585Project/serial_version.jl")

function test_e_step()
  psi = build_psi();
  theta = [1. 3. 5. 7.;2. 4. 6. 8.];
  a_prime = expectation(psi,theta);
  @test_approx_eq_eps(a_prime[5],0.0,0.1);
end

function test_m_step(w)
  @assert length(w)==4;
  psi = build_psi();
  theta = [1. 3. 5. 7.;2. 4. 6. 8.];
  x_prime,s_prime = maximization(psi,theta,w);
  @test x_prime[5][1]<psi[5,3][1];
end

function test_run_algorithm_2D(iterations::Int64=100,samples::Int64=1000,epsilon::Float64=0.1)
  @assert iterations>0;
  @assert samples>0;
  @assert epsilon>0;
  test_index=1;
  psi=run_algorithm(build_psi(),MvNormal([3.,4.],[2. 2.; 3. 3.]),iterations,samples);
  for i in 1:length(psi[:,1])
    if psi[i,1]>psi[test_index,1]
      test_index = i;
    end
  end
  @test_approx_eq_eps(psi[test_index,3][1],3.0,epsilon);
end

function test_run_algorithm_3D(iterations::Int64=100,samples::Int64=1000,epsilon::Float64=0.1)
  @assert iterations>0;
  @assert samples>0;
  @assert epsilon>0;
  test_index=1;
  psi=build_psi(ones(5)./5,ones(Float64,5).*2,Array[[3.,4.,5.],[3.,4.,5.],[3.,4.,5.],[3.,4.,5.],[3.,4.,5.]],Matrix[[3. 2. 1.; 3. 2. 1.; 3. 2. 1.],[3. 2. 1.; 3. 2. 1.; 3. 2. 1.],[3. 2. 1.; 3. 2. 1.; 3. 2. 1.],[3. 2. 1.; 3. 2. 1.; 3. 2. 1.],[3. 2. 1.; 3. 2. 1.; 3. 2. 1.]]);
  psi=run_algorithm(psi,MvNormal([3.,4.,5.],[2. 2. 2.; 3. 3. 3.; 4. 4. 4.]),iterations,samples);
  for i in 1:length(psi[:,1])
    if psi[i,1]>psi[test_index,1]
      test_index = i;
    end
  end
  @test_approx_eq_eps(psi[test_index,3][1],3.0,epsilon);
end


srand(123)
@time test_run_algorithm_2D(1000,10000,0.1)

