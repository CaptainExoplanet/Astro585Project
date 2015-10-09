using Distributions
srand(123)

function iterate(p_i,psi::Array{4},N::Int=1000)
  w = zeros(N);
  sum = 0.0;
  theta = zeros(N);
  for i in 1:N
    theta[i] = cumsum(psi[:,1].*rand(psi[:,4])); #This might be wrong
    density_func = cumsum(psi[:,1].*pdf(psi[:,4],theta[i]));
    w[i] = exp(logpdf(p_i,theta[i])-log(density_func));
    sum += w[i];
  end
  w_norm = w[i]/sum;
  #delete_comp(q);
  #merge_comp(q);
  #add_comp(q);
  new_psi = update_comp(w_norm,psi);
  return new_psi;
end

function main_func(p::Array,q,big_i::Int=1000)
  i=1
  while i<big_i
    psi[i] = iterate(p[i],q,psi[i-1])
    i++
  end
end

# setup psi
function build_psi(alpha::Array{Float64,1}=ones(10)./10,df::Array{Int8,1}=ones(Int8,10).*2,x::Array{Float64,1}=[1.,2.,3.,4.,5.,6.,7.,8.,9.,10.])
  q = Normal[]
  for i in 1:length(df)
    q=vcat(q,Normal(x[i],df[i]))
  end
  psi = [alpha df x q]
end

