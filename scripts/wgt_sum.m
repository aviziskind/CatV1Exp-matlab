function y = wgt_sum(x,w, normFlag)
    x = x(:); w = w(:);           
    
    if exist('normFlag', 'var') && ~isempty(normFlag)
        w = w/sum(w);
    end        
    y = sum( x .* w );
end

% formula for weighted variance:
%     w = w/sum(w);
%     x0 = w.*x;
%     x = x - sum(x0);
%     var = sum(w .* x^2)
