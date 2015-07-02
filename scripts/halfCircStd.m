function cs = halfCircStd(r, phi, phi_mean)
    if nargin < 2
        n = length(r);
        phi = linspace(0, pi, n+1);
        phi(end) = [];
        phi = reshape(phi, size(r));
    end
    if nargin < 3
        phi_mean = .5 *(circMean(r, 2*phi));
    end
    
    phi_diff = abs(phi-phi_mean);
    phi_diff = min(phi_diff, pi-phi_diff);
    
    cs = sqrt( sum((phi_diff .^2) .* r) / sum(r) );        
end