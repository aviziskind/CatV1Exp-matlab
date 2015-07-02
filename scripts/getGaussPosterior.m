function   [post, log_lh,mahalaD] = getGaussPosterior(X, mu, Sigma, p)
%WDENSITY Weighted conditional density and mahalanobis distance.
%   LOG_LH = WDENSITY(...) returns log of component conditional density
%   (weighted by the component probability) of X. LOG_LH is a N-by-K matrix
%   LOG_LH, where K is the number of Gaussian components. LOG_LH(I,J) is
%   log (Pr(point I|component J) * Prob( component J))
%
%   [LOG_LH, MAHALAD]=WDENSITY(...) returns the Mahalanobis distance in
%   the N-by-K matrix MAHALAD. MAHALAD(I,J) is the Mahalanobis distance of
%   point I from the mean of component J.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:20:40 $

    log_prior = log(p);
    [n,d]=size(X);
    k=size(mu,1);
    log_lh = zeros(n,k);
    mahalaD = zeros(n,k);
    
    for j = 1:k
                
        % compute the log determinant of covariance
        [L,f] = chol(Sigma(:,:,j) );
        diagL = diag(L);
        if (f ~= 0) || any(abs(diagL) < eps(max(abs(diagL)))*size(L,1))
            error('stats:gmdistribution:wdensity:IllCondCov', ...
                'Ill-conditioned covariance created');
        end
        logDetSigma = 2*sum(log(diagL));
        

        Xcentered = bsxfun(@minus, X, mu(j,:));
       
        xRinv = Xcentered /L ;
        mahalaD(:,j) = sum(xRinv.^2, 2);

        log_lh(:,j) = -0.5 * mahalaD(:,j) +...
            (-0.5 *logDetSigma + log_prior(j)) - d*log(2*pi)/2;
        %get the loglikelihood for each point with each component
        %log_lh is a N by K matrix, log_lh(i,j) is log \alpha_j(x_i|\theta_j)
    end

    [ll, post] = estep(log_lh);
    
end
   



function  [ll, post, logpdf]=estep(log_lh)
%ESTEP E-STEP for Gaussian mixture distribution
%   LL = ESTEP(LOG_LH) returns the loglikelihood of data in LL.  LOG_LH
%   is the log of component conditional density weighted by the component
%   probability.
%
%   [LL, POST] = ESTEP(LOG_LH) returns the posterior probability in the
%   matrix POST. POST(i,j) is the posterior  probability of point i
%   belonging to cluster j.
%
%   [LL, POST, DENSITY] = ESTEP(LOG_LH) returns the logs of the pdf values
%   of data in the vector density.
%
%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:20:37 $

maxll = max (log_lh,[],2);
%minus maxll to avoid underflow
post = exp(bsxfun(@minus, log_lh, maxll));
%density(i) is \sum_j \alpha_j P(x_i| \theta_j)/ exp(maxll(i))
density = sum(post,2);
%normalize posteriors
post = bsxfun(@rdivide, post, density);
logpdf = log(density) + maxll;
ll = sum(logpdf) ;


end