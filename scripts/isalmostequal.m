function tf = isalmostequal(a, b, eps)    
    tf = (length(a) == length(b)) && (max( abs(a-b) ) < eps);
end