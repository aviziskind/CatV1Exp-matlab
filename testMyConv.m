% function testMyConv
clear Y1 Y2;
n = 16;
x = rand(n,1);

c = gaussian(-3:3, 0, .6)'; c = c/sum(c);
m = length(c);

y1 = conv(x,c);
y2 = myConv(x,c);

% a = ceil(m/2)-1;
% figure(1); clf; hold on;
% plot(a+[1:n], x, 'bo-'); 
% plot(1:m+n-1, y1, 'g.:'); 

3;
% sum(abs(y1-y2))
assert( sum(abs(y1-y2)) < 1e-10)

ny = 5;
nz = 10;

% X = rand(n,ny);
% for iy = 1:ny
%      Y1(:,iy) = conv(X(:,iy), c);    
% end

X = rand(n,ny, nz);
for iy = 1:ny
    for iz = 1:nz
        y = conv(X(:,iy, iz), c);
        Y1(:,iy, iz) = y; %(a+[1:n]);    
    end
end
Y2 = myConv(X, c, 1);    
assert ( sum(abs(Y1(:)-Y2(:))) < 1e-10 )


clear Y1 Y2;
w = 2;
X = rand(n,ny, nz);
tic;
for iy = 1:ny
    for iz = 1:nz
        y = gaussSmooth(X(:,iy, iz), w);
        Y1(:,iy, iz) = y; %(a+[1:n]);    
    end
end
t1 = toc;

tic;
Y2 = gaussSmooth(X, w, 1);
t2 = toc;
t1/t2

assert( sum(abs(Y1(:)-Y2(:))) < 1e-10)
3;


% test gaussian smoothing: padding with 1, end, vs no padding
% 
% x = gaussSmooth(rand(n,1)*3, 3)+randn(n,1)/2;
% c = gaussian(-4:4, 0, 2)'; c = c/sum(c);
% m = length(c);
% 
% y_func = gaussSmooth(x, 2); % with padding
% y_conv = myConv(x,c, 'same');  % no padding
% 
% % figure(1); clf; hold on;
% plot([1:n], x, 'bo-'); 
% plot([1:n], y_func, 'gs-'); 
% plot([1:n], y_conv, 'r.:'); 






