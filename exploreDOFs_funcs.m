function exploreDOFs_funcs


%     y =  exp( -(((x-a).^2)./(2*b^2)+((x-c).^4)./(2*d^4)));  

    N = 20;
    nT = 50000;
    
    theta = linspace(0, 2*pi, N+1); theta = theta(1:end-1);
    
%     freq1 = 1;
%     freq2 = 2;
%     freq3 = 3;
%     f1 = sin(freq1*theta)';
%     f2 = sin(freq2*theta)';
%     f3 = sin(freq3*theta)';
% 
%     ampFreq1s = randn(2, nT);
%     ampFreq2s = randn(2, nT);
%     ampFreq3s = randn(2, nT);
%             
%     v1 = f1 * ampFreq1s(1,:) + f2 * ampFreq2s(1,:) + f3 * ampFreq3s(1,:);
%     v2 = f1 * ampFreq1s(2,:) + f2 * ampFreq2s(2,:) + f3 * ampFreq2s(2,:);



    freq1 = 1;
    freq2 = 2;
    freq3 = 3;
    f1 = sin(freq1*theta)';
    f2 = sin(freq2*theta)';
    f3 = sin(freq3*theta)';

    ampFreq1s = randn(2, nT);
    ampFreq2s = randn(2, nT);
    ampFreq3s = randn(2, nT);
            
    v1 = f1 * ampFreq1s(1,:) + f2 * ampFreq2s(1,:) + f3 * ampFreq3s(1,:);
    v2 = f1 * ampFreq1s(2,:) + f2 * ampFreq2s(2,:) + f3 * ampFreq2s(2,:);


    ccs = pearsonR_v(v1, v2);
    meanCCs = mean(ccs);    
    stdCCs = std(ccs);    
    dof = 1/(stdCCs^2);

    figure(1); clf;
    hist(ccs, 50);
    title(sprintf('%.2f \\pm %.2f (%.2f)', meanCCs, stdCCs, dof));

end

%     muRange = [0, N];
%     sigRange = [1, 5];

%         mu1 = rand(1,nT)*diff(muRange)+muRange(1);
%         mu2 = rand(1,nT)*diff(muRange)+muRange(1);
%         sig = 1;%rand(1,nT)*diff(sigRange)+sigRange(1);
%                 
%         v1 = getGaussianVecs(N,mu1,sig);
%         v2 = getGaussianVecs(N,mu2,sig);



        
%{
                bias_prop = bias_props(bi);
                progressBar;    

                abs_bias_prop = abs(bias_prop);
                sgn_bias_prop = real(bias_prop);

                x1w = bsxfun(@plus, (1-abs_bias_prop)*x1, abs_bias_prop*(bias_vec));
                x2w = bsxfun(@plus, (1-abs_bias_prop)*x2, abs_bias_prop*(bias_vec)*sgn_bias_prop);            

                ccs = pearsonR_v(x1w, x2w);                
                meanCCs(ni, vi, bi) = mean(ccs);    
                stdCCs(ni, vi, bi) = std(ccs);    




    % ga

    Ns = [4, 5, 6, 7, 8, 9, 10];    
%     Ns = [3, 5, 10, 15, 20, 30, 50];    
    nNs = length(Ns);
    nT = 10000;

    bias_props = [0:.01:1];
%     bias_props = [0:.01:.4, 0.405:.005:.8 ,.81:.01:1];
    
    nB = length(bias_props);
%     bias_vec = randn(N, 1);    
    
    
    nVs = 10;    
    bias_vecs = cell(1,nNs);
    for ni = 1:nNs
        bias_vecs{ni} = randn(Ns(ni), nVs);
        bias_vecs{ni} = bsxfun(@minus, bias_vecs{ni}, mean(bias_vecs{ni},1));
    end

    
    figure(102); clf; hold on;
    
    meanCCs = zeros(nNs,nVs,nB);
    stdCCs = zeros(nNs,nVs,nB);
    
    progressBar('init', nNs*nVs*nB, 40);
    for ni = 1:nNs
        N = Ns(ni);
    
        for vi = 1:nVs
            bias_vec = bias_vecs{ni}(:,vi);        
            x1 = randn(N,nT); x1 = bsxfun(@minus, x1, mean(x1, 1));
            x2 = randn(N,nT); x2 = bsxfun(@minus, x2, mean(x2, 1));

            for bi = 1:nB
    %             ccs = zeros(nT,1);
                bias_prop = bias_props(bi);
                progressBar;    

                abs_bias_prop = abs(bias_prop);
                sgn_bias_prop = real(bias_prop);

                x1w = bsxfun(@plus, (1-abs_bias_prop)*x1, abs_bias_prop*(bias_vec));
                x2w = bsxfun(@plus, (1-abs_bias_prop)*x2, abs_bias_prop*(bias_vec)*sgn_bias_prop);            

                ccs = pearsonR_v(x1w, x2w);                
                meanCCs(ni, vi, bi) = mean(ccs);    
                stdCCs(ni, vi, bi) = std(ccs);    
            end
        
        end
        
    end
        

    figure(102); clf; hold on;
    hs = zeros(nNs, nVs);
    for ni = 1:nNs
        for vi = 1:nVs        
            x = meanCCs(ni,vi,:);
            y = stdCCs(ni,vi,:);
            hs(ni, vi) = plot(x(:), y(:), ['.', color_s(ni)]);
        end
    end

%     legend(hs(:,1), legendarray('N = ', Ns));
    3;
    
    hh(1) = plot(wvfm_raw1(1), wvfm_raw1(2), 'ko', 'markersize', 3, 'markerfacecolor', 'k');
    hh(2) = plot(wvfm_raw2(1), wvfm_raw2(2), 'ko', 'markersize', 5, 'markerfacecolor', 'k');
    hh(3) = plot(wvfm_raw3(1), wvfm_raw3(2), 'ko', 'markersize', 7, 'markerfacecolor', 'k');
    hh(4) = plot(wvfm_raw(1), wvfm_raw(2), 'ks', 'linewidth', 1, 'markerfacecolor', 'k');
    hh(5) = plot(wvfm_ccw(1), wvfm_ccw(2), 'rs', 'linewidth', 1, 'markerfacecolor', 'r');
    

    leg_strs = [legendarray('N = ', Ns); {'raw(1)', 'raw(2)', 'raw(3)', 'raw', 'ccw'}'];
    legend([hs(:,1);hh'], leg_strs, 'fontsize', 9 );
    xlabel('mean'); ylabel('std');
    box on;
    axis([-.01 1.1, -.01 .7]);
    
    3;
    
%     figure(100); hold on;
%     hold on; 
%     plot(bias_props, meanCCs, 'k.-', bias_props, stdCCs, 'b.-'); legend('mean', 'std')


end

%}
function v = getGaussianVecs(N,mu,sigma, nvec)
    if (nargin < 4) 
        nvec = length(mu);
    end
    x = [1:N]';
    
    if length(mu) < nvec
        mu = mu(ones(1,nvec));
    end
    if length(sigma) < nvec
        sigma = sigma(ones(1,nvec));
    end
        
    v = zeros(N,nvec);
    for i = 1:nvec
        v(:,i) = gaussian(x, mu(i), sigma(i));
    end
end
