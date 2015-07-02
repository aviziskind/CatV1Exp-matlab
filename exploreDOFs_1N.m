function exploreDOFs_1N

%     wvfm_ccw = [.459, .300];
%     wvfm_raw = [.723, .176];
% 
%     wvfm_raw1 = [.847, .134];
%     wvfm_raw2 = [.782, .144];
%     wvfm_raw3 = [.746, .161];


    N = 6;
    nT = 1000;
    cc_mean = 0;

%     bias_props = [0:.1:1];
%     bias_props = [0:.01:1];
    bias_props = [0:.01:.4, 0.405:.005:.8 ,.81:.01:1];
    nB = length(bias_props);
%     bias_vec = randn(N, 1);    
    
    
    
    nVs = 8;
    bias_vecs = randn(N, nVs);
    bias_vecs = bsxfun(@minus, bias_vecs, mean(bias_vecs,1));

    cc_halfs_mean = zeros(1,nVs);
    cc_pows_mean   = zeros(1,nVs);
    cc_off_mean   = zeros(1,nVs);
    
    cc_halfs_std = zeros(1,nVs);
    cc_pows_std   = zeros(1,nVs);
    cc_off_std    = zeros(1,nVs);
    
    std_zero_bias = 1./ sqrt(N-1);
    sig_func_pos = @(beta, x) sigmoid(x-beta(3), 1, beta(1), beta(2));
    sig_func_neg = @(beta, x) std_zero_bias*(1-sigmoid(x-beta(3), 1, beta(1), beta(2)));

    beta0 = [.5, 1, 0];    
    
    progressBar('init-', nB*nVs)
    figure(100); clf; hold on;
    figure(101); clf; hold on;

    
    
    meanCCs = zeros(nVs,nB);
    stdCCs = zeros(nVs,nB);
    skewCCs = zeros(nVs,nB);
    kurtCCs = zeros(nVs,nB);
    
    for vi = 1:nVs
        bias_vec = bias_vecs(:,vi);        

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
            if nB < 5
                figure(bi); clf; hist(ccs, 30);
            end
            meanCCs(vi, bi) = mean(ccs);    
            stdCCs(vi, bi) = std(ccs);    
            skewCCs(vi,bi) = skewness(ccs);
            kurtCCs(vi,bi) = kurtosis(ccs);
        end
        
        beta_mean = nlinfit(bias_props, meanCCs(vi,:), sig_func_pos, beta0);
        cc_halfs_mean(vi) = beta_mean(1);
        cc_pows_mean(vi) = beta_mean(2);
        beta_std = nlinfit(bias_props, stdCCs(vi,:), sig_func_neg, beta0);
        cc_halfs_mean(vi) = beta_std(1);
        cc_pows_mean(vi) = beta_std(2);
                
    end
        

    figure(100); clf; hold on;
    figure(101); clf; hold on;
    figure(102); clf; hold on;
    figure(103); clf; hold on;
    figure(104); clf; hold on;
    for vi = 1:nVs
        beta_mean = nlinfit(bias_props, meanCCs(vi,:), sig_func_pos, beta0);
        cc_halfs_mean(vi) = beta_mean(1);
        cc_pows_mean(vi) = beta_mean(2);
        cc_off_mean(vi) = beta_mean(3);
        
        beta_std = nlinfit(bias_props, stdCCs(vi,:), sig_func_neg, beta0);
        cc_halfs_std(vi) = beta_std(1);
        cc_pows_std(vi) = beta_std(2);
        cc_off_std(vi) = beta_std(3);

        figure(100); 
        plot(bias_props, meanCCs(vi,:), ['.', color_s(vi)]);
        fplot(@(x) sig_func_pos(beta_mean, x), lims(bias_props), [color_s(vi)]);        
        
        figure(101); 
        plot(bias_props, stdCCs(vi,:), ['.', color_s(vi)]);
        fplot(@(x) sig_func_neg(beta_std, x), lims(bias_props), [color_s(vi)]);     
        
        figure(102);
        plot(meanCCs(vi,:), stdCCs(vi,:), ['.', color_s(vi)])
        
        figure(103);
        plot(bias_props, skewCCs(vi,:), ['.', color_s(vi)]);        

        figure(104);
        plot(bias_props, kurtCCs(vi,:), ['.', color_s(vi)]);                
    end
    
%     figure(100); box on; xlabel('bias'); ylabel('mean')
%     title(sprintf('Mean cc vs bias for N = %d, for %d vectors', N, nVs));
% 
%     figure(101); box on; xlabel('bias'); ylabel('std')
%     title(sprintf('Std cc vs bias for N = %d, for %d vectors', N, nVs));
% 
%     figure(103); box on; xlabel('bias'); ylabel('skewness')
%     title(sprintf('Skewness of cc vs bias for N = %d, for %d vectors', N, nVs));
% 
%     figure(104); box on; xlabel('bias'); ylabel('kurtosis')
%     title(sprintf('Kurtosis of cc vs bias for N = %d, for %d vectors', N, nVs));
%     
%     figure(105); plot(meanCCs', stdCCs', '.'); xlabel('mean'); ylabel('std');
%     figure(106); plot(meanCCs', skewCCs', '.'); xlabel('mean'); ylabel('skewness');
%     figure(107); plot(stdCCs', skewCCs', '.'); xlabel('std'); ylabel('skewness');
%     
    3;
    bias_means = mean(bias_vecs,1);
    bias_stds = std(bias_vecs,[], 1)
    
    3;
    
%     figure(100); hold on;
%     hold on; 
%     plot(bias_props, meanCCs, 'k.-', bias_props, stdCCs, 'b.-'); legend('mean', 'std')


end



