function jd = jointDistr(x12, Nbin);
jd = [];
jd = zeros(Nbin,Nbin);
M = size(x12,1);

x = x12(:,1);
y = x12(:,2);

x_max = unique(max(x));
x_min = unique(min(x));
dx = (x_max - x_min)/(Nbin-1);
if dx == 0
  return;
end

y_max = unique(max(y));
y_min = unique(min(y));
dy = (y_max - y_min)/(Nbin-1);
if dy == 0
  return;
end

x = x - x_min;
y = y - y_min;

for k = 1:M
  ix = fix(x(k)/dx)+1;
  iy = fix(y(k)/dy)+1;
  jd(iy,ix) = jd(iy,ix) + 1;
end

jd = jd/sum(jd(:));
