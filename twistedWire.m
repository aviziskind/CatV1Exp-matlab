
nPerTwist = 100;
nTwists = 1;
z = linspace(0, nTwists, nTwists*nPerTwist);
theta = linspace(0, 2*pi*nTwists, nPerTwist*nTwists);

w = .5;

nWires = 4;
phis = [0:nWires-1]*pi/2;

x = cell(1,nWires);
y = cell(1,nWires);
for i = 1:nWires
    x{i} = w*cos(theta+phis(i));  
    y{i} = w*sin(theta+phis(i));
end
 

figure(101); clf; hold on;
for i = 1:nWires
    h(i) = plot3(x{i}, y{i}, z, color_s(i));
end
axis([-1, 1, -1, 1, 0 nTwists]);
set(h, 'linewidth', 5, 'marker', 'o', 'linestyle', '-')
view(3);
%%

