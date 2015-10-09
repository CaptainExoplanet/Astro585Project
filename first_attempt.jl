using Distributions
srand(123)

function iterate(p_i,q,psi,N::Int=1000)
  w = zeros(N);
  sum = 0.0;
  theta = rand(q,N);
  for i in 1:N
    w[i] = exp(logpdf(p_i,theta[i])-logpdf(q,theta[i]));
    sum += w[i];
  end
  w_norm = w[i]/sum;
  delete_comp(q);
  merge_comp(q);
  add_comp(q);
  new_psi = update_comp(q);
  return new_psi;
end

function main_func(p::Array,q,big_i::Int=1000)
  i=1
  while i<big_i
    psi[i] = iterate(p[i],q,psi[i-1])
    i++
  end
end

# setup q
M=10
df=2
N=1000
q = TDist[]
for i in 1:N
  q=vcat(q,TDist(df))
end
alpha = ones(N)./N
#for i in 1:M
#big_q = alpha[1] +
