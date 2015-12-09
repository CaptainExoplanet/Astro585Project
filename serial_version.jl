#= This module provides functions that complete a portion of the algorithm described in Liu (2014). This version
of the algorithm runs through the expectation and maximization steps to update the components of the proposal
distribution. For testing purposes, see the testing_utilities.jl module.
=#

using Distributions

#= This function calculates the posterior mass as described in Lui (2014) section 3.2
    - Inputs:
      -an array of psi containing the mixture distribution
      -a random sample from the mixture distribution
    - Ouputs:
      -an array containing posterior masses
=#
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

#= This function calculates the u_m as described in Lui (2014) section 3.4
    - Inputs:
      -an array of psi containing the mixture distribution
      -a random sample from the mixture distribution
      -an index identifying the component of psi
    - Ouputs:
      -u_m (a scalar)
=#
calc_u_m(psi,theta,j) = (psi[j,2]+length(psi[:,1]))/(psi[j,2]+reshape((theta.-psi[j,3])'*psi[j,4]*(theta.-psi[j,3]),1)[1])

#= This function calculates the C_n as described in Lui (2014) section 3.4
    - Inputs:
      -an array of psi containing the mixture distribution
      -a random sample from the mixture distribution
      -an index identifying the component of psi
    - Ouputs:
      -C_n (a matrix)
=#
calc_C_n(psi,theta,j) = (theta.-psi[j,3])*(theta.-psi[j,3])'

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

#= This function initializes psi as a 2-dimensional array of varying types
    - Inputs:
      -an array of probability masses
      -an array of degrees of freedom
      -a two-dimensional array of component means
      -an array of covariance matrices
    - Outputs:
      -an array of psi containing the mixture distribution
=#
function build_psi(alpha=ones(5)./5,df=ones(Float64,5).*2,x=Array[[1.,2.],[3.,4.],[5.,6.],[7.,8.],[9.,10.]],sigma=Matrix[[1. 0.5; 0.5 1.],[1. 0.5; 0.5 1.],[1. 0.5; 0.5 1.],[1. 0.5; 0.5 1.],[1. 0.5; 0.5 1.]])
  q = MvTDist[];
  for i in 1:length(df)
    q=vcat(q,MvTDist(df[i],vec(x[i]),sigma[i]));
  end
  psi = [alpha df x sigma q];
end

#= This function updates the parameters described in the expectation section of Lui (2014) section 3.4
    - Inputs:
      -an array of psi containing the mixture distribution
      -a random sample from the mixture distribution
    - Outputs:
      -an updated probability mass array
    - Functions called
      -calc_epsilon()
=#
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

#= This function updates the parameters described in the maximization section of Lui (2014) section 3.4
    - Inputs:
      -an array of psi containing the mixture distribution
      -a random sample from the mixture distribution
      -an array of weights representing the importance of each component of the mixture distribution
    - Outputs:
      -an updated mean array
      -an updated covariance matrix array
    - Functions called
      -calc_epsilon()
      -calc_u_m()
      -calc_C_n()
=#
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

#= This function updates the components of psi as desribed in Liu (2014) section 3.4
     - Inputs:
       -a random sample from the mixture distribution
       -an array of weights representing the importance of each component of the mixture distribution
       -an array of psi containing the mixture distribution
     - Outputs:
       -an updated version of psi
     - Functions called:
       -expectation()
       -maximization()
       -build_psi()
=#
function update_comp(theta,w,psi::Array{Any,2}) #add types
  #Expectation
  alpha_prime = expectation(psi,theta);
  #Maximization
  x_prime,sig_prime = maximization(psi,theta,w)
  return build_psi(alpha_prime,psi[:,2],x_prime,sig_prime);
end

#= This function runs through the algorithm a number of times specified by iterations
    - Inputs:
      -an array of psi containing the mixture distribution
      -a density function for the target population
      -the number of iterations
      -the number of samples to be taken at each iteration
    - Outputs:
      -a final version of psi
    - Functions called:
      -iterate()
=#
function run_algorithm(psi=build_psi(),p=MvNormal([3.,4.],[2. 2.;3. 3.]),iterations::Int64=100,samples::Int64=1000)
  @assert iterations>0;
  @assert samples>0;
  for i in 1:iterations
    psi=iterate(p,psi,samples);
  end
  return psi;
end
