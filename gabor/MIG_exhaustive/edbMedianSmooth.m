function y = edbMedianSmooth(x,n);
Npnt = length(x);
for k = 1:Npnt
  l = n;
  r = n;
  if k <= n
    l = k-1;
  end
  if k > Npnt - n
    r = Npnt-k;
  end
  v = x(k-l:k+r);
  y(k) = median(v);
end
