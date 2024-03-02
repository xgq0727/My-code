
function T=gen_ZadoffChu_Seq(N,r,K)
 T =zeros(N,K);
 w=exp(-1j*2*pi*r/N);
for q = 1 : K
  k=[0:1:N-1];
 
  if mod(N,2)==0
      Z=w.^(k.*k/2+q*k);
  else
      Z=w.^(k.*(k+1)/2+q*k);
  end
  T(:,q) = Z.';
end