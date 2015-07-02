function [h, p, z] = testSignificanceOfMagOfSum(X, sumDim, shuffleDim, alpha, showWorkingFlag)
    if nargin < 4
        alpha = .05;
    end
    if ~strcmp(class(X), 'double');
        X = double(X); % for 'mag' function
    end
    showWorking = exist('showWorkingFlag', 'var') && ~isempty(showWorkingFlag);
    
    n_shuffles = 100;
    
    r_actual = norm( sum(X, sumDim) );
    rs_shuffled = zeros(1, n_shuffles);
    for rep_i = 1:n_shuffles                
        rs_shuffled(rep_i) = norm( sum( shuffle(X, shuffleDim), sumDim));
    end
    
    mu = mean(rs_shuffled);
    sigma = std(rs_shuffled);    
    
    [h, p] = ztest(r_actual, mu, sigma, alpha, 'right');
    
%     if nargout > 2
        z = (r_actual - mu)/sigma;
        z = max(0, z);
%     end
 
    if showWorking
        hist(rs_shuffled);
        drawVerticalLine(r_actual, 'color', 'r');
        drawVerticalLine(mu, 'color', 'g');
        drawVerticalLine([mu + sigma,mu - sigma], 'color', 'g','linestyle', ':');
        drawVerticalLine(median(rs_shuffled), 'color', 'b');
        title( sprintf('h = %d, p = %2.4g (z = %2.3g)', h, p, z));
    end    
    
end
 