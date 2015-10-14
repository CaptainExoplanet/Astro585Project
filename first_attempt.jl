# In one dimension for now

using Distributions
srand(123)

function iterate(p_i,psi::Array{4},N::Int=1000)
  w = zeros(N);
  w_norm = zeros(N);
  sum = 0.0;
  theta = zeros(N);
  for i in 1:N
    theta[i] = cumsum(psi[:,1].*rand(psi[:,4])); #This might be wrong
    density_func = cumsum(psi[:,1].*pdf(psi[:,4],theta[i]));
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

# setup psi
function build_psi(alpha::Array{Float64,1}=ones(10)./10,df::Array{Int8,1}=ones(Int8,10).*2,x::Array{Float64,1}=[1.,2.,3.,4.,5.,6.,7.,8.,9.,10.])
  q = Normal[];
  for i in 1:length(df)
    q=vcat(q,Normal(x[i],df[i]));
  end
  psi = [alpha df x q];
end

function update_comp(theta,w,psi::Array{4},dim::Int=1) #add types
  epsilon = zeros(length(psi[:,1]),length(theta));
  alpha_prime = zeros(length(psi[:,1]));
  for i in 1:length(psi[:,1])
    for j in 1:length(theta)
      epsilon[i,j] = psi[i,1]*pdf(psi[i,4],theta[j])*length(theta)/theta[j]; #also could be wrong
      alpha_prime[i] += pdf(psi[i,4],theta[j])*epsilon[i,j]
      #should have an x_prime and sigma calculation here
    end
  end
  return build_psi(alpha_prime,psi[:,2],psi[:,3])
end

function test_simple(N::Int=100)
  psi=build_psi();
  p=Normal(5,3)
  for i in 1:N
    psi=iterate(p,psi);
  end
  return psi;
end

test_simple()



#Tests with kepler planet are below. These don't work yet.
#cd("$(homedir())/Astro585Project/TTVFaster/Julia/benchmark")
#include("TTVFaster/Julia/test_ttv.jl")
#data=readdlm("$(homedir())/Astro585Project/TTVFaster/Julia/kepler62ef_planets.txt",',',Float64)
#data=vec(data)
#ttv1,ttv2=test_ttv(5,40,20,data)
#include("TTVFaster/Julia/benchmark/benchmark_grad_ttv.jl")

