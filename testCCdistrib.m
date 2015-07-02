function testCCdistrib
% Investigating the factors that affect the distribution of cc/rho between two vectors
% A.
%    1. dependence on length of vector (n), or parameters of distribution (mean/std/skew)
% B. dependence on #degrees of freedom
%    1  compare gaussian profile with shuffled gaussian profile:
%    2  look at correlation function of actual tuning curves
%    3  with real data, compare distribution before/after shuffling tuning curves


% 4. is distribution of means centered at zero?

cmpFunc = @pearsonR;
% cmpFunc = @spearmanRho;

tic;
fmt = '%.2f';
str_mn_std = @(x) [num2str(nanmean(x), fmt) ' \pm ' num2str(nanstd(x)/sqrt(nnz(~isnan(x))), fmt)];

doTests = [1];


if any(doTests == 1)   % 1. test dependence on 1st, 2nd, 3rd momnts
    % test dependence on local structure.
    
    n = 8;  % length of vector
    B = 10000;
    
%     elementDistrib = 'gaussian';
    elementDistrib = 'exponential';
    
    % use gaussian.    
    ccs1 = zeros(1,B);
    ccs2 = zeros(1,B);
    sig = 2;
    % progressBar('init', B, 30);
    tic;    
    
    % control mean / var / skewness / n:
    param_name = 'n';

    plotDistribs = true;
    plotMoments = false;
    
%     p = [1:20:80]; % parameter values
%     p = [2 4 8, 16, 32, 64, 128]; % parameter values
%     p = 50;
%     p = [2 3 4 5 6 8 16 32 64]; % parameter values
    p = [4 8 16]; % parameter values
    nP = length(p);
    ccs = zeros(B, nP);
    progressBar('init-', nP, nP);
    
    doIndiv = n*B > 1e7;
    sizeVec = iff(doIndiv, [n,1], [n,B]);
    
    switch elementDistrib
        case 'gaussian',    genfun = @(sizeVec) randn(sizeVec);
        case 'exponential', genfun = @(sizeVec) exprnd(1, sizeVec);
        case 'uniform',     genfun = @(sizeVec) rand(sizeVec);
    end
    
    for p_i = 1:nP
        progressBar(p_i);
        switch param_name
            case 'mean'
                mu = p(p_i);
                sig = 1;
%                 genfun = @() randn(n,B)*sig + mu;
                
                genfun = @() exprnd(mu, sizeVec);
            case 'var'
                mu = 0;
                sig = p(p_i);
                if sig == 0
                    sig = .1;
                end
                genfun = @() randn(sizeVec)*sig + mu;
            case 'skew'
                %want to find mu, sigma such that var = V, skewness = S.
                V = p(p_i);
                S = 5;
                
                sigma = fzero(@(s) skw(s)-S, 1);
                mu = fzero(@(mu) vr(mu, sigma)-V, 1);                                
                genfun = @() exp(randn(sizeVec)*sig + mu);
                
            case 'n'
%                 mu = 1;
%                 sig = 4;
                n = p(p_i);
                
                doIndiv = n*B > 1e7;
                sizeVec = iff(doIndiv, [n,1], [n,B]);
                
%                 genfun = @() randn(sizeVec)*sig + mu;
%                 genfun = @() exp(randn(sizeVec)*sig + mu);
%                 genfun = @() exprnd(mu, sizeVec);
        
                genfun2 = @() genfun(sizeVec);
                
            otherwise
                error('wrong paramname');
        end
        
        if ~doIndiv % generate all random numbers at once
            X = genfun2();
            Y = genfun2();
        end
        
        for b = 1:B
            if doIndiv
                x = genfun2();
                y = genfun2();
            else
                x = X(:,b);
                y = Y(:,b);
            end
                        
            ccs(b, p_i) = cmpFunc(x, y);
            
        end
    end
    progressBar('done');
    
    fprintf('Computed all cc''s'); toc; 
    
    if plotDistribs, nBins = 31;
        binEdges = linspace(-1, 1, nBins+1);
        binCenters = mean([binEdges(1:nBins); binEdges(2:nBins+1)],1);
        for p_i = 1:nP
            nBins = 50;
            
%             X = /(sum(ccs(:,p_i))/50);
            figure(p_i);
            n = histcnt(ccs(:,p_i), binEdges); 
            bar(binCenters, n/(sum(n)), 'hist');
            xlim([-1 1]); title(str_mean_median(ccs(:,p_i)));
            xlabel([ param_name ' = ' num2str(p(p_i), '%d')]);
        end
    end
     return;
    if plotMoments
        funs = {@mean, @std, @skewness, @kurtosis};
        nF = length(funs);
        M = zeros(nF, nP);
        E = zeros(nF, nP);
        figure(11); clf;
        for fi = 1:nF
            tic;
            fprintf([' computing '  func2str(funs{fi}) ' : ' ]); 
            nboot = 1000;
            h(fi) = subplot(1,nF, fi);
            for p_i = 1:nP
                fprintf('.');
                b = bootstrp(nboot, funs{fi}, ccs(:, p_i));
                M(fi, p_i) = mean(b);
                E(fi, p_i) = std(b); 
            end
            errorbar(p, M(fi,:), E(fi,:), 'bo-'); hold on        

            title([func2str(funs{fi})]);
            xlabel(param_name);
            toc; 
        end
    end
end

return;
if any(doTests == 4)
    for b = 1:B
    %     progressBar;
        mu1 = n/2;
        mu2 = rand*n;
    %     x = gaussian(0:n-1, mu1, sig)';
    %     y = gaussian(0:n-1, mu2, sig)';


        ccs1(b) = cmpFunc(x,y);
        xp = x(randperm(n));
        yp = y(randperm(n));

        ccs2(b) = cmpFunc(x,y);


    end
    figure(1);
    hist(ccs1, 20);
    title(str_moments(ccs1, 4))
    figure(2);
    title(str_moments(ccs2, 4))

end

if any(doTests == 5)
    
    % test dependence on local structure.
    n = 8;
    B = 1000;
    BB = 1000;
    % use gaussian.

    meanCCs = zeros(1, BB);
    for bb = 1:BB
    
    
%     p = [1:20:80]; % parameter values
%     p = [2 4 8, 16, 32, 64, 128]; % parameter values
        ccs = zeros(1, B);
%     p = 50;
%         progressBar('init-', nP, nP);
    
        doIndiv = n*B > 1e7;
        sizeVec = iff(doIndiv, [n,1], [n,B]);
    
        mu = 1;
        sig = 1;
        genfun = @() exprnd(mu, sizeVec);
                
        
        X = genfun();
        Y = genfun();
    end
        
    for b = 1:B
        if doIndiv
            x = genfun();
            y = genfun();
        else
            x = X(:,b);
            y = Y(:,b);
        end

        ccs(b, p_i) = cmpFunc(x, y);
            
        end
    end
    progressBar('done');
    
    fprintf('Computed all cc''s'); toc; 
    
    if plotDistribs, nBins = 31;
        binEdges = linspace(-1, 1, nBins+1);
        binCenters = mean([binEdges(1:nBins); binEdges(2:nBins+1)],1);
        for p_i = 1:nP
            nBins = 50;
            
%             X = /(sum(ccs(:,p_i))/50);
            figure(p_i);
            n = histcnt(ccs(:,p_i), binEdges); 
            bar(binCenters, n/(sum(n)), 'hist');
%             xlim([-1 1]); title(str_mean_median(ccs(:,p_i)));
            xlabel([ param_name ' = ' num2str(p(p_i), '%d')]);
        end
    end

    
    toc;
end





function s = str_mean_median(x)

    fmt = '%.3f';
        
    x = x(~isnan(x));
    N = length(x);
    mean_x = mean(x);
    stderr_x = stderr(x); 
    
    median_x = median(x);
%     med_samples = bootstrp(100, @median, x);
%     std_median1 = stderr(med_samples);
    
    med_ci = bootci(100, @median, x );
%     std_median2 = diff(std_med)/2;
    
    [h,p_t] = ttest(x);
    p_W = signrank(x);    
    
    s1 = sprintf([fmt ' \\pm ' fmt ' || ' fmt ' [ ' fmt ', ' fmt ' ] (%d)'], mean_x, stderr_x, median_x, med_ci(1), med_ci(2), N);
    s2 = sprintf(['p_t = ' fmt '. p_W = ' fmt], p_t, p_W);
    s = {s1, s2};
    
end