function [samp_mean, samp_std] = calcNullF1oDCdistribution(allNBins, allNSamples, Fn, f1Type, reflect, nBoot)
    
    nnbins = length(allNBins);
    nnsamples = length(allNSamples);
    
    samp_mean = zeros(nnbins,nnsamples);
    samp_std = zeros(nnbins,nnsamples);
    
    if nargin < 6
        nBoot = 500;    
    end
    
    doBar = nnbins * nnsamples > 20;
    if doBar, progressBar('init-', nnbins * nnsamples, 30); end;
    
    for bi = 1:nnbins
        nbins = allNBins(bi);
        
        if reflect
            idx = [1:nbins, nbins-1:-1:2];
            nbins_ext = length(idx);            
            binEdges = linspace(0,360, nbins_ext+1);
            binCenters_refl = binEdge2cent(binEdges);            
        else
            binEdges = linspace(0,360, nbins+1);
            binCenters = binEdge2cent(binEdges);
        end            
        
        for si = 1:nnsamples
            if doBar, progressBar; end
            nsamples = allNSamples(si);

            F1oDC_control = zeros(1, nBoot);
            for i = 1:nBoot
                randBinSamples = randi(nbins, 1, nsamples);
                [uVals_tmp, uBinCount] = unique(randBinSamples);
                binCount = zeros(1,nbins);
                binCount(uVals_tmp) = uBinCount;
                if reflect
                    binCount_refl = binCount(idx);
                    f1odc = getF1oDC(binCenters_refl, binCount_refl, 360, Fn*2, f1Type);
                else                    
                    f1odc = getF1oDC(binCenters, binCount, 360, Fn, f1Type);                    
                end
                F1oDC_control(i) = f1odc;
            end
            [n,xout] = normhist(F1oDC_control, 30);

            m = mean(F1oDC_control);
            s = std(F1oDC_control);
            switch f1Type
                case 'cos', 
                    se = stderr(F1oDC_control);
                    assert( ibetween(m, -5*se, 5*se) );
                    samp_mean(bi,si) = m;
                    samp_std(bi,si) = s;
                case 'F1',
                    if nsamples > 10
                        beta = nlinfit(xout, n, @(beta, x) beta(3)*gaussian(x, beta(1), beta(2)), [m, s, n(1)]);
                        samp_mean(bi,si) = beta(1);
                        samp_std(bi,si) = abs(beta(2));
                    else
                        samp_mean(bi,si) = m;
                        samp_std(bi,si) = s;
                    end
                    
            end

        end    
    end
    
end
