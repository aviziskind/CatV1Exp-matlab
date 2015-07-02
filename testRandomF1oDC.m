allNBins = [5, 11];
allNSamples = [5, 10, 15, 20, 25, 30, 40, 50:10:500];

nBoot = 1000;

nnbins = length(allNBins);
nnsamples = length(allNSamples);

means = zeros(nnbins, nnsamples);
stds = zeros(nnbins, nnsamples);
just_stds = zeros(nnbins, nnsamples);
h= zeros(nnbins, nnsamples);

showHistograms = false;

if showHistograms
    figure(243); clf;
    plot_j = 1;
end

progressBar('init-', nnbins*nnsamples, 30);
for bi = 1:nnbins
    
    nbins = allNBins(bi);
    
    for si = 1:nnsamples
        progressBar;
        nsamples = allNSamples(si);
        
        binEdges = linspace(0,360, nbins+1);
        binCenters = binEdge2cent(binEdges);
        
        F1oDC_control = zeros(1, nBoot);
        for i = 1:nBoot
            randBinSamples = randi(nbins, 1, nsamples);
            [uVals_tmp, uBinCount] = unique(randBinSamples);
            binCount = zeros(1,nbins);
            binCount(uVals_tmp) = uBinCount;
%             valsControl_refl = binCount(idx);
            f1odc = getF1oDC(binCenters, binCount, 360, 1, 'cos' );
            F1oDC_control(i) = f1odc;
        end
        [n,xout] = normhist(F1oDC_control, 30);
            m = mean(F1oDC_control);
            s = std(F1oDC_control);
            just_stds(bi,si) = s;
        
        beta = nlinfit(xout, n, @(beta, x) beta(3)*gaussian(x, beta(1), beta(2)), [m, s, n(1)]);
        if showHistograms
            h(bi, si) = subplot(nnbins,nnsamples, plot_j);
            bar(xout, n, 1); hold on;            
            fplot(@(x) beta(3)*gaussian(x, beta(1), beta(2)), [xout(1), xout(end)], 'r:')
            s1 = sprintf('[nBin = %d. nSamp = %d.]', nbins, nsamples);
            s2 = sprintf('m = %.2f, s = %.2f', beta(1), beta(2));
            title({s1, s2});
            plot_j = plot_j +1;
%             xlim([-1, 1])
%             ylim([0 12]);
        end
        
        
        means(bi, si) = beta(1);
        stds(bi, si) = abs(beta(2));
        
    end    
    
    
end
matchAxes('XY', h);
progressBar('done');
   
% figure(244); clf;
% L = legendarray('N = ', allNSamples);
% subplot(1,2,1); plot(allNBins, means, 'o-'); title('mean'); xlabel('nbins'); legend(L, 'location', 'SE')as
% subplot(1,2,2); plot(allNBins, stds, 'o-'); title('std'); xlabel('nbins'); legend(L, 'location', 'SE'); hold on;
% plot(allNBins, just_stds, 's:'); 

figure(245); clf;
L = legendarray('nbin = ', allNBins);
subplot(1,2,1); plot(allNSamples, means', 'o-'); title('mean'); xlabel('nsamples'); legend(L, 'location', 'NE')
subplot(1,2,2); plot(allNSamples, stds', 'o'); title('std'); xlabel('nsamples'); legend(L, 'location', 'NE'); hold on;
plot(allNSamples, just_stds', 's'); 
for bi = 1:nnbins
    expfun = @(beta, X) beta(1).*exp(-X./beta(2))+beta(3);
    beta = nlinfit(allNSamples, just_stds(bi,:), expfun, [1 1]);
    fplot(@(x) expfun(beta, x), xlim, [color_s(bi) ':']);    
end
