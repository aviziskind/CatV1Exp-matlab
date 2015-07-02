function [X_norm, stdev] = normalizeToNStd(X, stdev) 
    if isvector(X)
        X = X(:)';
    end
    X_ms = bsxfun(@minus, X, mean(X,2));
    if (nargin < 2) || isempty(stdev)
        stdev = std_med(X, 2);
    end
    X_norm = bsxfun(@rdivide, X_ms, stdev);
end