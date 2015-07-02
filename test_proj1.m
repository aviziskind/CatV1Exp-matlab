function test_proj1



x = [-5:.2:5];
y = [-5:.2:5];

[xg, yg] = meshgrid(x,y);
zg = d_actual(xg, yg);


figure(41);
surf(x, y, zg');

figure(42);

% f = @(x,y)  sqrt( (x- (x+y)/2).^2 + (y-(x+y)/2).^2 );
f = @(x,y)  sqrt( .5*(x-y)^2 );
fmesh(f, [-5, 5], [-5, 5])

3;

end

function z = d_actual(x,y)

    z = zeros(size(x));
    
    
    for i = 1:numel(x)

        if x(i) == y(i)
           3; 
        end
%     x = [3; 4];
        H = 1/2*[1 1; 1 1];
        xp = H*[x(i);y(i)];
        3;

    z(i) =   norm([x(i)-xp(1), y(i)-xp(2)]);

    end
end