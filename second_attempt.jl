using Distributions
srand(123)

#= This function iterates psi as desribed in Liu (2014) section 3
     - Inputs:
       -a density funtion for the target population
       -an array of psi containing the mixture distribution
       -The number of iterations
     - Outputs:
       -An updated version of psi
     - Functions called:
       -update_comp()
=#
function iterate(p_i,psi::Matrix,N::Int64=1000)
  w = zeros(N);
  w_norm = zeros(N);
  sum = 0.0;
  theta = zeros(N);
  density_func = 0.0
  for i in 1:N
    for j in 1:length(psi[:,1])
      theta[i] += psi[j,1]*rand(psi[j,5]);
    end
    for j in 1:length(psi[:,1])
      density_func += psi[j,1]*pdf(psi[j,5],theta[i]);
    end
    w[i] = exp(logpdf(p_i,theta[i])-log(density_func));
    sum += w[i];
  end
  for i in 1:N
    w_norm[i] = w[i]/sum;
  end
  #delete_comp(q);
  #merge_comp(q);
  #add_comp(q);
  new_psi = update_comp(theta,w_norm,psi);
  return new_psi;
end

# This function initializes psi as a 2-dimensional array of varying types
function build_psi(alpha::Array{Float64,1}=ones(5)./10,df::Array{Float64,1}=ones(Float64,5).*2,x=Array[[1.,2.],[3.,4.],[5.,6.],[7.,8.],[9.,10.]],sigma=Matrix[[2. 1.; 2. 1.],[2. 1.; 2. 1.],[2. 1.; 2. 1.],[2. 1.; 2. 1.],[2. 1.; 2. 1.]])
  q = MvTDist[];
  for i in 1:length(df)
    q=vcat(q,MvTDist(df[i],vec(x[i]),sigma[i]));
  end
  psi = [alpha df x sigma q];
end

psi=build_psi()
psi[:,1]

""" This function updates the components of psi as desribed in Liu (2014) section 3.4
     - Inputs:
       -a random sample from the mixture distribution
       -an array of psi containing the mixture distribution
       -The dimensions of the problem
     - Outputs:
       -An updated version of psi
     - Functions called:
       -build_psi
"""
function update_comp(theta,w,psi::Array{Any,2},dim::Int=1) #add types
  epsilon = zeros(length(psi[:,1]),length(theta));
  alpha_prime = zeros(length(psi[:,1]));
  for i in 1:length(psi[:,1])
    for j in 1:length(theta)
      epsilon[i,j] = psi[i,1]*pdf(psi[i,4],theta[j])*length(theta)/theta[j];
      alpha_prime[i] += pdf(psi[i,4],theta[j])*epsilon[i,j]
      #should have an x_prime and sigma calculation here
    end
  end
  return build_psi(alpha_prime/sum(alpha_prime),psi[:,2],psi[:,3])
end

#Testing my functions
psi=build_psi()
p=Normal(5,3)
for i in 1:10
  psi=iterate(p,psi,1000);
end
psi


#Tests with kepler planet are below. These don't work yet.
#cd("$(homedir())/Astro585Project/TTVFaster/Julia/benchmark")
#include("TTVFaster/Julia/test_ttv.jl")
#data=readdlm("$(homedir())/Astro585Project/TTVFaster/Julia/kepler62ef_planets.txt",',',Float64)
#data=vec(data)
#ttv1,ttv2=test_ttv(5,40,20,data)
#include("TTVFaster/Julia/benchmark/benchmark_grad_ttv.jl")

