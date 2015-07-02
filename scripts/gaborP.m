function [Z, J] = gaborP(P, X) 
    Z = gabor(P(1), P(2), P(3), P(4), P(5), P(6), P(7), P(8), P(9), X);
    
    if nargout > 1
        J = dgabor(P(1), P(2), P(3), P(4), P(5), P(6), P(7), P(8), P(9), X);
    end    
    
end