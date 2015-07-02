function [a,b,c] = divrem(x, y)
    % x / y =>  a*y + b
    a = floor( (x-1)./y) ;
    b = rem( (x-1), y) + 1;
    c = (a)*y + b;
end
