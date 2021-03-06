using Distributions
using Base.Test
srand(123)

function iterate(p,psi::Matrix,N::Int64=1000)
  @assert N>0;
  w = zeros(N);
  w_norm = zeros(N);
  sum = 0.0;
  theta = [];
  density_func = 0.0;
  lamda=linspace(0,1,10)
  for i in 1:N
    for j in 1:length(psi[:,1])
      if j==1
        if i==1
          theta = vcat(theta,psi[j,1]*rand(psi[j,5]));
        else
          theta = hcat(theta,psi[j,1]*rand(psi[j,5]));
        end
      else
        theta[:,i] = theta[:,i] .+ psi[j,1]*rand(psi[j,5]);
      end
    end
    for j in 1:length(psi[:,1])
      density_func += psi[j,1]*pdf(psi[j,5],vec(theta[:,i]));
    end
    for k in 1:length(lamda)
      w[i] += exp(lamda[k]*logpdf(p,theta[:,i])+(1-lamda[k])*log(pdf(psi[1,5],vec(theta[:,i])))-log(density_func));
    end
    sum += w[i];
  end
  for i in 1:N
    w_norm[i] = w[i]/sum;
  end
  new_psi = update_comp(theta,w_norm,psi);
  return new_psi;
end

# This function initializes psi as a 2-dimensional array of varying types
function build_psi(alpha=ones(5)./5,df=ones(Float64,5).*2,x=Array[[1.,2.],[3.,4.],[5.,6.],[7.,8.],[9.,10.]],sigma=Matrix[[2. 1.; 2. 1.],[2. 1.; 2. 1.],[2. 1.; 2. 1.],[2. 1.; 2. 1.],[2. 1.; 2. 1.]])
  q = MvTDist[];
  for i in 1:length(df)
    q=vcat(q,MvTDist(df[i],vec(x[i]),sigma[i]));
  end
  psi = [alpha df x sigma q];
end

function calc_epsilon(psi,theta)
  epsilon = zeros(length(psi[:,1]),length(theta[1,:]));
  for i in 1:length(theta[1,:])
    for j in 1:length(psi[:,1])
      epsilon[j,i] = psi[j,1]*pdf(psi[j,5],vec(theta[:,i]));
    end
    epsilon[:,i] = epsilon[:,i]/sum(epsilon[:,i])
  end
  return epsilon
end

calc_u_m(psi,theta,j) = (psi[j,2]+length(psi[:,1]))/(psi[j,2]+reshape((theta.-psi[j,3])'*psi[j,4]*(theta.-psi[j,3]),1)[1])

calc_C_n(psi,theta,j) = (theta.-psi[j,3])*(theta.-psi[j,3])'

function expectation(psi,theta)
  epsilon = calc_epsilon(psi,theta);
  alpha_prime = zeros(length(psi[:,1]));
  for j in 1:length(psi[:,1])
    for i in 1:length(theta[1,:])
      alpha_prime[j] += pdf(psi[j,5],vec(theta[:,i]))*epsilon[j,i];
    end
  end
  return alpha_prime/sum(alpha_prime);
end

function maximization(psi,theta,w)
  epsilon = calc_epsilon(psi,theta);
  x_prime = Array(Array,length(psi[:,1]));
  sig_prime = Array(Matrix,length(psi[:,1]));
  for j in 1:length(psi[:,1])
    top_x = 0.0;
    bottom_x = 0.0;
    top_sig = Array(Float64,length(theta[:,1]),length(theta[:,1]));
    top_sig = top_sig.*0.0;
    bottom_sig = 0.0;
    for i in 1:length(theta[1,:])
      top_x += w[i]*epsilon[j,i]*calc_u_m(psi,theta[:,i],j).*theta[:,i];
      bottom_x += w[i]*epsilon[j,i]*calc_u_m(psi,theta[:,i],j);
      top_sig += w[i]*epsilon[j,i]*calc_u_m(psi,theta[:,i],j).*calc_C_n(psi,theta[:,i],j);
      bottom_sig += w[i]*epsilon[j,i];
    end
    x_prime[j]=top_x./bottom_x;
    sig_prime[j]=top_sig./bottom_sig;
  end
  return x_prime,sig_prime
end

function update_comp(theta,w,psi::Array{Any,2}) #add types
  #Expectation
  alpha_prime = expectation(psi,theta);
  #Maximization
  x_prime,sig_prime = maximization(psi,theta,w)
  return build_psi(alpha_prime,psi[:,2],x_prime,sig_prime);
end

function run_algorithm(psi=build_psi(),p=MvNormal([3.,4.],[2. 2.;3. 3.]),iterations::Int64=100,samples::Int64=1000)
  @assert iterations>0;
  @assert samples>0;
  for i in 1:iterations
    psi=iterate(p,psi,samples);
  end
  return psi;
end

#### Testing my functions ####

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

psi=build_psi()
bs=linspace(2,3,2)
bss=[bs bs]
p=MvNormal([9.,10.],bss)
@time for i in 1:10
  psi=iterate(p,psi,100);
end
psi

@time test_run_algorithm_3D(10,100,1.0)

test_e_step()
test_m_step([.2,.2,.2,.2])


#Tests with kepler planet are below. These don't work yet.
#cd("$(homedir())/Astro585Project/TTVFaster/Julia/benchmark")
#include("TTVFaster/Julia/test_ttv.jl")
#data=readdlm("$(homedir())/Astro585Project/TTVFaster/Julia/kepler62ef_planets.txt",',',Float64)
#data=vec(data)
#ttv1,ttv2=test_ttv(5,40,20,data)
#include("TTVFaster/Julia/benchmark/benchmark_grad_ttv.jl")

