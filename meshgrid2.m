function [x_g, y_g] = meshgrid2(x, y)

    nx = length(x);
    ny = length(y);

%     x_g = zeros(1, nx*ny);
%     y_g = zeros(1, nx*ny);
    x_g = zeros(ny, nx);   %[x=1, y=1], [x=1, y=2], ...[x=1], [y=n]
    y_g = zeros(ny, nx);

    for xi = 1:nx
        for yi = 1:ny
            x_g(yi + (xi-1)*ny ) = x(xi);
            y_g(yi + (xi-1)*ny ) = y(yi);            
        end
    end
    
    [x_g1, y_g1] = meshgrid(x,y);
    [x_g1 - x_g,  y_g1 - y_g]

end