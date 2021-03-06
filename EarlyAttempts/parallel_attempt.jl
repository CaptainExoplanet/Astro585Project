nprocs()
addprocs(3)
proclist = workers()         # Get list of workers

# test that can use tasks on all workers
for proc in proclist
  println("Julia says welcome from proc ",proc," w/ hostname ",@fetchfrom(proc,gethostname() ))
end

@everywhere using Distributions
srand(123)

function iterate(p,psi::Matrix,N::Int64=1000)
  w = zeros(N);
  w_norm = zeros(N);
  sum = 0.0;
  theta = [];
  density_func = 0.0;
  lamda=linspace(0,1,10)
  sum = @parallel (+) for i in 1:N
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
    w[i];
  end
  w_norm = w/sum;
  #delete_comp(q);
  #merge_comp(q);
  #add_comp(q);
  new_psi = update_comp(theta,w_norm,psi);
  return new_psi;
end

# This function initializes psi as a 2-dimensional array of varying types
function build_psi(alpha=ones(5)./10,df=ones(Float64,5).*2,x=Array[[1.,2.],[3.,4.],[5.,6.],[7.,8.],[9.,10.]],sigma=Matrix[[2. 1.; 2. 1.],[2. 1.; 2. 1.],[2. 1.; 2. 1.],[2. 1.; 2. 1.],[2. 1.; 2. 1.]])
  q = MvTDist[];
  for i in 1:length(df)
    q=vcat(q,MvTDist(df[i],vec(x[i]),sigma[i]));
  end
  psi = [alpha df x sigma q];
end

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

@everywhere function calc_u_m(psi,theta)
  u_m = zeros(length(psi[:,1]));
  for j in 1:length(psi[:,1])
    u_m[j] = (psi[j,2]+length(psi[:,1]))/(psi[j,2]+reshape((theta.-psi[j,3])'*psi[j,4]*(theta.-psi[j,3]),1)[1]);
  end
  return u_m;
end

@everywhere calc_C_n(psi,theta,j) = (theta.-psi[j,3])*(theta.-psi[j,3])'

function update_comp(theta,w,psi::Array{Any,2},dim::Int=1) #add types
  epsilon = calc_epsilon(psi,theta);
  alpha_prime = zeros(length(psi[:,1]));
  x_prime = Array(Array,length(psi[:,1]));
  sig_prime = Array(Matrix,length(psi[:,1]));
  #Expectation
  for j in 1:length(psi[:,1])
    alpha_prime[j] = @parallel (+) for i in 1:length(theta[1,:])
      pdf(psi[j,5],vec(theta[:,i]))*epsilon[j,i];
    end
  end
  #Maximization
  for j in 1:length(psi[:,1])
    top_x = 0.0;
    bottom_x = 0.0;
    top_sig = Array(Float64,2,2);
    top_sig = top_sig.*0.0;
    bottom_sig = 0.0;
    top_x = @parallel (+) for i in 1:length(theta[1,:])
      w[i]*epsilon[j,i]*calc_u_m(psi,theta[:,i])[j].*theta[:,i];
    end
    bottom_x = @parallel (+) for i in 1:length(theta[1,:])
      w[i]*epsilon[j,i]*calc_u_m(psi,theta[:,i])[j];
    end
    top_sig = @parallel (+) for i in 1:length(theta[1,:])
      w[i]*epsilon[j,i]*calc_u_m(psi,theta[:,i])[j].*calc_C_n(psi,theta[:,i],j);
    end
    bottom_sig = @parallel (+) for i in 1:length(theta[1,:])
      w[i]*epsilon[j,i];
    end
    x_prime[j]=top_x./bottom_x;
    sig_prime[j]=top_sig./bottom_sig;
  end
  return build_psi(alpha_prime/sum(alpha_prime),psi[:,2],x_prime,sig_prime);
end

#Testing my functions
psi=build_psi()
bs=linspace(1,2,2)
bss=[bs bs]
p=MvNormal([3.,4.],bss)
@time for i in 1:10
  psi=iterate(p,psi,1000);
end


#Tests with kepler planet are below. These don't work yet.
#cd("$(homedir())/Astro585Project/TTVFaster/Julia/benchmark")
#include("TTVFaster/Julia/test_ttv.jl")
#data=readdlm("$(homedir())/Astro585Project/TTVFaster/Julia/kepler62ef_planets.txt",',',Float64)
#data=vec(data)
#ttv1,ttv2=test_ttv(5,40,20,data)
#include("TTVFaster/Julia/benchmark/benchmark_grad_ttv.jl")
