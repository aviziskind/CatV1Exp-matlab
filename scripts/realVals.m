function y = realVals(x) 
    y = x(~isinf(x) & ~isnan(x));
end