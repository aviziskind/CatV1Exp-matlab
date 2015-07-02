% function temp1

figure(1);
h1 = subplot(1,2,1); axis([-1, 1, -1, 1]); zeroaxes;
h2 = subplot(1,2,2); 

first = true;
xs = [0:.1:1]; n = length(xs);
ys = [0:.1:1]; 
w0 = 1;
[xs_g, ys_g] = meshgrid(xs, ys);
X = [xs_g(:)'; ys_g(:)'];
w = [1; 1];




while true
   
    z = w' * X;
    z = reshape(z, [n,n]);
    if first
        subplot(h2);
        h = mesh(xs,ys, z);
        axis(axis);
        zlim([-1, 1]);
        first = false;
    else
        set(h, 'zdata', z)
    end
    subplot(h1);
    [x,y] = ginput(1);
    w = [x;y];
end






% end
% 