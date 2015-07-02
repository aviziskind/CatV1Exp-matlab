function [x_s, y_s, z_s] = sort3(x,y,z)
    [x_s, idx] = sort(x);
    y_s = y(idx);
    z_s = z(idx);
end