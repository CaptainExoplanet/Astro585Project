function iterate(p_i,q,psi,N::Int=1000)
  w_i = zeros(N);
  sum = 0.0;
  for i in 1:N
    theta_n = calc_theta(i);
    w_i = p_i(theta_n)/q(theta_n,psi);
    sum += w_i;
  end
  w = w_i/sum;
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
