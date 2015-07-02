function y = gabor1D(A, mu, sigma, k, phi, X) 
    y = A*exp( -(X-mu).^2./(2*sigma.^2)) .* cos(abs(k)*(X-mu) + phi);    
end