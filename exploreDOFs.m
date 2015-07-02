function exploreDOFs

    fig_id = 110;

    doWvfmPts = 1;
    
    if doWvfmPts
    
    wvfm_ccw = [.459, .300];
    wvfm_raw = [.723, .176];
    wvfm_raw1 = [.847, .134];
    wvfm_raw2 = [.782, .144];
    wvfm_raw3 = [.746, .161];

    
    pca_raws_2_20 = [0.0555    0.9985
        0.0444    0.7177
        0.0432  0.5999;
        0.0426  0.5467;
        0.0438  0.5263;
        0.0433  0.5038;
        0.0428  0.4889;
        0.0421  0.4803;
        0.0422  0.4723
        0.0417    0.4651
        0.0415    0.4624
        0.0414    0.4605
        0.0414    0.4587
        0.0412    0.4576
        0.0411    0.4566
        0.0410    0.4557
        0.0411    0.4550
        0.0409    0.4541
        0.0410    0.4535];

    
    
    pca_ccws_2_20 = [    0.0271    0.9996
        0.0344    0.7259
        0.0382    0.6042;
        0.0316    0.5313;
        0.0250    0.4872;
        0.0232    0.4653;
        0.0225    0.4526;
        0.0230    0.4410;
        0.0232    0.4355
        0.0233    0.4321
        0.0231    0.4295
        0.0229    0.4255
        0.0226    0.4226
        0.0225    0.4203
        0.0221    0.4180
        0.0219    0.4166
        0.0217    0.4149
        0.0215    0.4135
        0.0214    0.4124    ];

     pca_rawSep_1_5 = [0.0000546  0.6044
    0.0452    0.4732
    0.0407    0.4597
    0.0403    0.4531
    0.0392    0.4478];

    glf_raw_2_20 = [0.8756    0.4830
        0.8863    0.1588            
        0.8450    0.1656
        0.8203    0.1707
        0.7966    0.1725
        0.7840    0.1714
        0.7745    0.1712
        0.7675    0.1708
        0.7616    0.1724
        0.7501    0.1747
        0.7437    0.1759
        0.7391    0.1759
        0.7337    0.1750
        0.7304    0.1752
        0.7275    0.1752
        0.7253    0.1751
        0.7223    0.1751
        0.7205    0.1751
        0.7166    0.1757];

        
    glf_ccw_2_20 = [-0.0026    1.0000
        0.0421    0.7077;
        0.3602    0.4818
        0.4352    0.3866
        0.4876    0.3369
        0.4799    0.3338
        0.4657    0.3274
        0.4607    0.3244
        0.4572    0.3229
        0.4522    0.3210
        0.4482    0.3176
        0.4440    0.3158
        0.4412    0.3141
        0.4395    0.3129
        0.4358    0.3107
        0.4339    0.3098
        0.4331    0.3090
        0.4341    0.3042
        0.4318    0.3031];

    end

    Ns = [2, 3, 4, 5, 6, 7, 8, 9, 10];    
%     Ns = [2, 3, 4, 5, 10, 20];    
%     Ns = [2:10, 100, 200, 400, 800, 1600];
%     Ns = [2, 3, 4];    
%     Ns = [3, 5, 10, 15, 20, 30, 50];    
    nNs = length(Ns);
    nT = 1000;

%     bias_props = [0:.01:1];
    bias_props = [0:.01:.4, 0.405:.005:.8 ,.81:.01:1];
    
    nB = length(bias_props);
%     bias_vec = randn(N, 1);    
    
    
    nVs = 10;    
    bias_vecs = cell(1,nNs);
    for ni = 1:nNs
%         bias_vecs{ni} = exprnd(1, Ns(ni), nVs);
%         bias_vecs{ni} = bsxfun(@minus, bias_vecs{ni}, mean(bias_vecs{ni},1));
        
        bias_vecs_i = exprnd(1, Ns(ni), nVs);
        bias_vecs_i = bsxfun(@minus, bias_vecs_i, mean(bias_vecs_i,1) );
        bias_vecs_i = bsxfun(@rdivide, bias_vecs_i, normV(bias_vecs_i, 1) );

        bias_vecs{ni} = bias_vecs_i;        
        
    end

        
    
    meanCCs = zeros(nNs,nVs,nB);
    stdCCs = zeros(nNs,nVs,nB);
    
    progressBar('init', nNs*nVs*nB, 40);
    for ni = 1:nNs
        N = Ns(ni);
    
        for vi = 1:nVs
            bias_vec = bias_vecs{ni}(:,vi);        

            for bi = 1:nB
    %             ccs = zeros(nT,1);
                bias_prop = bias_props(bi);
                progressBar;    

                abs_bias_prop = abs(bias_prop);
                x1 = exprnd(1, N,nT); x1 = bsxfun(@minus, x1, mean(x1, 1));
                x2 = exprnd(1, N,nT); x2 = bsxfun(@minus, x2, mean(x2, 1));                
                
                x1w = bsxfun(@plus, (1-abs_bias_prop)*x1, abs_bias_prop*(bias_vec));
                x2w = bsxfun(@plus, (1-abs_bias_prop)*x2, abs_bias_prop*(bias_vec));            

                ccs = pearsonR_v(x1w, x2w);                
                meanCCs(ni, vi, bi) = mean(ccs);    
                stdCCs(ni, vi, bi) = std(ccs);    
            end
        
        end
        
    end
        

    figure(fig_id); clf; hold on;
    hs = zeros(nNs, nVs);
    for ni = 1:nNs
        for vi = 1:nVs        
            x = meanCCs(ni,vi,:);
            y = stdCCs(ni,vi,:);
            hs(ni, vi) = plot(x(:), y(:), ['o', color_s(ni)], 'markersize', 3);
        end
    end

    
%     figure(fig_id+1); clf; hold on;
%     hs = zeros(nNs, nVs);
%     for ni = 1:nNs
%         y0 = 1./((Ns(ni)-1));
%         for vi = 1:nVs        
%             x = meanCCs(ni,vi,:);
%             y = (stdCCs(ni,vi,:).^2)/y0;
%             hs(ni, vi) = plot(x(:), y(:), ['o', color_s(ni)], 'markersize', 3);
%         end
%     end
    
    
%     legend(hs(:,1), legendarray('N = ', Ns));
    3;

    if doWvfmPts
        figure(fig_id);
        hh(1) = plot(wvfm_raw1(1), wvfm_raw1(2), 'ko', 'markersize', 3, 'markerfacecolor', 'k');
        hh(2) = plot(wvfm_raw2(1), wvfm_raw2(2), 'ko', 'markersize', 5, 'markerfacecolor', 'k');
        hh(3) = plot(wvfm_raw3(1), wvfm_raw3(2), 'ko', 'markersize', 7, 'markerfacecolor', 'k');
        hh(4) = plot(wvfm_raw(1), wvfm_raw(2), 'ks', 'linewidth', 1, 'markerfacecolor', 'k');
        hh(5) = plot(wvfm_ccw(1), wvfm_ccw(2), 'rs', 'linewidth', 1, 'markerfacecolor', 'r');
        
        h_pca_raw = plot(pca_raws_2_20(:,1), pca_raws_2_20(:,2), 'ko-', 'markerfacecolor', 'k', 'markersize', 3);
        h_pca_ccw = plot(pca_ccws_2_20(:,1), pca_ccws_2_20(:,2), 'ko-', 'markerfacecolor', 'r', 'markersize', 3);
        h_pca_rawS = plot(pca_rawSep_1_5(:,1), pca_rawSep_1_5(:,2), 'v-', 'markerfacecolor', 'k', 'markersize', 3);
        
        h_glf_raw = plot(glf_raw_2_20(:,1), glf_raw_2_20(:,2), 'ks-', 'markerfacecolor', 'k', 'markersize', 3);
        h_glf_ccw = plot(glf_ccw_2_20(:,1), glf_ccw_2_20(:,2), 'ks-', 'markerfacecolor', 'r', 'markersize', 3);
        
        set(h_pca_raw, 'markerfacecolor', 'k', 'markersize', 3);
        set(h_pca_ccw, 'markerfacecolor', 'r', 'markersize', 3);
        set(h_pca_rawS, 'markerfacecolor', 'r', 'markersize', 3);
        
        leg_strs = [legendarray('N = ', Ns-1); {'raw(1)', 'raw(2)', 'raw(3)', 'raw', 'ccw'}'];
        legend([hs(:,1);hh'], leg_strs, 'fontsize', 9 );
        
    else
        leg_strs = legendarray('N = ', Ns-1); 
        legend([hs(:,1)], leg_strs, 'fontsize', 9 );
    end
    
    xlabel('mean'); ylabel('std');
    box on;
    axis([-.05 1.02, -.01 1.02]);
    
    3;
    
%     figure(100); hold on;
%     hold on; 
%     plot(bias_props, meanCCs, 'k.-', bias_props, stdCCs, 'b.-'); legend('mean', 'std')


end



