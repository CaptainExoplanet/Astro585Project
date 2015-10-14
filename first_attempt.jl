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

function main_func(p::Array,q,big_i::Int=1000)
  i=2
  psi = build_psi();
  while i<big_i
    psi = iterate(p,psi)
    i++
  end
end

# setup psi
function build_psi(alpha::Array{Float64,1}=ones(10)./10,df::Array{Int8,1}=ones(Int8,10).*2,x::Array{Float64,1}=[1.,2.,3.,4.,5.,6.,7.,8.,9.,10.])
  q = Normal[];
  for i in 1:length(df)
    q=vcat(q,Normal(x[i],df[i]));
  end
  psi = [alpha df x q];
end

function update_comp(theta,w,psi,dim::Int=1) #add types
  epsilon = zeros(length(psi[:,1]),length(theta));
  alpha_prime = zeros(length(psi[:,1]));
  for i in 1:length(psi[:,1])
    for j in 1:length(theta)
      epsilon[i,j] = psi[i,1]*pdf(psi[i,4],theta[j])*length(theta)/theta[j]; #also could be wrong
      alpha_prime[i] += pdf(psi[i,4],theta[j])*epsilon[i,j]
      #should have an x_prime calculation here
    end
  end
  return build_psi(alpha_prime,psi[:,2],psi[:,3])
end

cd("$(homedir())/Astro585Project/")
