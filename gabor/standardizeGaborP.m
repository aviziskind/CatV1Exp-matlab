function params = standardizeGaborP(params)
    % ensure theta is between 0 and 180
    % ensure phi is between 0 and 360;

    idx_theta = 6;
    idx_phi = 8;    
    
    theta = params(idx_theta);
    phi   = params(idx_phi);
    
    phi = mod(phi, 2*pi);
    theta = mod(theta, 2*pi);
    if theta > pi
        theta = mod(theta, pi);  % ensure theta is between 0 and 180. (if between 180 and 360, 
%         phi = mod(phi+pi, 2*pi); % is the same as if we subtract 180, but with opposite phase sign).
        phi = 2*pi-phi; % is the same as if we subtract 180, but with opposite phase sign).
    end
    
    assert(ibetween(theta, 0, pi));
    assert(ibetween(phi, 0, 2*pi));
    
    params(idx_theta) = theta;
    params(idx_phi) = phi;
    
end

%{
    xs = [1:64];
    ys = [1:80];



%}