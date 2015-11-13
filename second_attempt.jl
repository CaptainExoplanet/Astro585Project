using Distributions
srand(123)

function iterate(p,psi::Matrix,N::Int64=1000)
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
bs=[]
bs=vcat(bs,psi[2,1]*rand(psi[2,5]))
bs=hcat(bs,psi[2,1]*rand(psi[2,5]))
bs=hcat(bs,psi[2,1]*rand(psi[2,5]))
bs[1,2]

function calc_epsilon(psi,theta)
  epsilon = zeros(length(psi[:,1]),length(theta[1,:]));
  density_func = 0.0;
  for i in 1:length(theta[1,:])
    for j in 1:length(psi[:,1])
      density_func += psi[j,1]*pdf(psi[j,5],vec(theta[:,i]));
      epsilon[j,i] = psi[j,1]*pdf(psi[j,5],vec(theta[:,i]));
    end
  end
  return epsilon/density_func;
end

#not finished product
function calc_u_m(psi,theta)
  u_m = zeros(length(psi[:,1]));
  for j in 1:length(psi[:,1])
    u_m[j] = (psi[j,2]+length(psi[:,1]))/(psi[j,2]+(theta[j,:].-psi[j,3])*psi[j,4]*(theta[j,:].-psi[j,3]));
  end
  return u_m;
end

#not finished product
function calc_C_n(psi,theta)
  return ones(vec(psi[:,4])); #needs fixing
end

function update_comp(theta,w,psi::Array{Any,2},dim::Int=1) #add types
  epsilon = calc_epsilon(psi,theta);
  u_m = calc_u_m(psi,theta);
  C_n = calc_C_n(psi,theta);
  alpha_prime = zeros(length(psi[:,1]));
  x_prime = zeros(length(psi[:,3]));
  #Expectation
  for j in 1:length(psi[:,1])
    for i in 1:length(theta[1,:])
      alpha_prime[j] += pdf(psi[j,5],vec(theta[:,i]))*epsilon[j,i];
    end
  end
  #Maximization
  for j in 1:length(psi[:,1])
    top_x = 0.0;
    bottom_x = 0.0;
    top_sig = 0.0;
    bottom_sig = 0.0;
    for i in 1:length(theta[1,:])
      top_x += w[i]*epsilon[j,i]*u_m[i]*theta[j,i];
      bottom_x += w[i]*epsilon[j,i]*u_m[i];
      top_sig += w[i]*epsilon[j,i]*u_m*C_n[i];
      bottom_sig += w[i]*epsilon[j,i];
    end
    x_prime[j]=top_x/bottom_x;
    sig_prime[j]=top_sig/bottom_sig;
  end
  return build_psi(alpha_prime/sum(alpha_prime),psi[:,2],x_prime,sig_prime);
end

#Testing my functions
psi=build_psi()
bs=linspace(1,2,2)
bss=[bs bs]
p=MvNormal([3.,4.],bss)
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

