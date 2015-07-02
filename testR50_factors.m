function testR50_factors
%%
        exp_r50_factor_1D = 0.6931;
        exp_r50_factor_2D = 1.6785;
        exp_r50_factor_3D = 2.6741;
        gauss_r50_factor_1D = 0.6745;
        gauss_r50_factor_2D = 1.17741; %sqrt(2*log(2));
        gauss_r50_factor_3D = 1.538172;

        exp_r50_factors   = [exp_r50_factor_1D, exp_r50_factor_2D, exp_r50_factor_3D];
        gauss_r50_factors = [gauss_r50_factor_1D, gauss_r50_factor_2D, gauss_r50_factor_3D];

        exp_r50_factors   = [0.6931, 1.6785, 2.6741];
        gauss_r50_factors = [0.6745, 1.17741, 1.538172];
        
        %%
    doD1_test_analytical = 0;
    if doD1_test_analytical
        D = 1;
        N = 100000;
        mu = 1;
        r = exprnd(mu, D, N);
        %%
        mn = mean(r, 2);

        med = median(r, 2);
        med_pred = mu * exp_r50_factor_1D;

        %%
        sig = 1;
        g = abs(normrnd(0, sig, D, N));
    %     st = std(g);
        med = median(g, 2);
        med_pred = sig * gauss_r50_factor_1D;
    end
        
    %%
    exp_mus = [.3:.1:.8];
%     exp_mus = [.5:.1:1]; 
    gauss_sig = exp_mus; 
    n = length(exp_mus);
    
    exp_means = zeros(1,n); 
    exp_medians = zeros(1,n); 
    gauss_means = zeros(1,n); 
    gauss_medians = zeros(1,n);

    D = 3;
    dx = 0.01;
    
    if D == 1
        x = [0:dx:10];
    elseif D == 2;
        x = [0:.02:10]; y = x;
        [x_grid, y_grid] = meshgrid(x,y);
        r_grid = sqrt(x_grid.^2 + y_grid.^2);    
    elseif D == 3;
        %%
        x = [0:.03:5]; y = x; z = x;
        [x_grid, y_grid, z_grid] = meshgrid(x,y,z);
        r_grid = sqrt(x_grid.^2 + y_grid.^2 + z_grid.^2);    
    end
    
    for i = 1:n
        
        if D == 1
            F_exp = exppdf(x * exp_r50_factors(D) ,exp_mus(i));
            F_gauss = normpdf(x * gauss_r50_factors(D), 0, gauss_sig(i)); 
        elseif D == 2 || D == 3
            
            F_exp = exppdf(r_grid * exp_r50_factors(D) ,exp_mus(i));
            F_gauss = normpdf(r_grid * gauss_r50_factors(D), 0, gauss_sig(i));             
            
        end
        
        
        F_exp_norm = F_exp/sum(F_exp(:));                
        F_gauss_norm = F_gauss(:)/sum(F_gauss(:));
        
        if D == 1
            exp_medians(i) = x( find(cumsum(F_exp_norm) > 0.5, 1) );
            gauss_medians(i) = x( find(cumsum(F_gauss_norm) > 0.5, 1) );            
            
        elseif D == 2 || D == 3
        
            [r_sorted, r_srt_idx] = sort(r_grid(:));
            F_exp_norm_srt = F_exp_norm(r_srt_idx);
            F_gauss_nrm_sort = F_gauss_norm(r_srt_idx);
            
            exp_medians(i) = r_sorted( find(cumsum(F_exp_norm_srt) > 0.5, 1) );
            gauss_medians(i) = r_sorted( find(cumsum(F_gauss_nrm_sort) > 0.5, 1) );
            
        end

%%        
%         exp_means(i) = wgt_sum(r_grid(:), F_exp_norm(:));
%%
%         gauss_means(i) = wgt_sum(x, F_gauss);

        
        
        
        
%         
%         figure(55);  clf;
%         plot(x, F_exp);
%         drawVerticalLine(exp_means(i), 'color', 'b')
%         drawVerticalLine(exp_medians(i), 'color', 'r')
%         drawVerticalLine(exp_mus(i), 'color', 'g', 'linestyle', ':')
%         3;
%         figure(56);  clf;
%         plot(x, F_gauss);
%         drawVerticalLine(gauss_means(i), 'color', 'b')
%         drawVerticalLine(gauss_medians(i), 'color', 'r')
%         drawVerticalLine(gauss_sig(i), 'color', 'g', 'linestyle', ':')
%         3;
        
        
        
    end    
    %%
    figure(1 + (10*D)); clf;  
    subplot(1,2,1); hold on;    
    fplot(@(x) x, lims(exp_mus), 'k-');
%     plot(exp_mus, exp_means, 'bo-');
    
%     plot(exp_mus, exp_medians / exp_r50_factors(D), 'r.-');
    plot(exp_mus, exp_medians, 'r.-');
    xlabel('mu'); ylabel('median / factor');
    
    title('Exponential');

    subplot(1,2,2); hold on;    
    fplot(@(x) x, lims(exp_mus), 'k-');
%     plot(gauss_sig, gauss_means, 'bo-');

%     plot(gauss_sig, gauss_medians / gauss_r50_factors(D), 'r.-');
    plot(gauss_sig, gauss_medians , 'r.-');
    xlabel('sigma'); ylabel('median / factor');
%     
    title('Gaussian');    
    3;
    
    %%
    
    D = 2;
    x = [-3:.002:3];  y = x; z = x;
    exp_mus = [.3:.1:.8];
    [x_grid, y_grid] = meshgrid(x,y);
    r_grid = sqrt(x_grid.^2 + y_grid.^2);
    for i = 1:length(exp_mus);    
        F_exp = exppdf(r_grid,exp_mus(i));
        F_exp = F_exp(:);        
        exp_means(i) = mean(F_exp);
        exp_medians(i) = median(F_exp);
    end
    figure(5); clf; hold on;
    plot(exp_mus, exp_means, 'bo-');
    plot(exp_mus, exp_medians, 'r.-');
%     hist(F_exp, 100);
    3;
    %%
    
    
    x_pos = x(x>=0); y_pos = x_pos; 
    
    
3;



  


end

