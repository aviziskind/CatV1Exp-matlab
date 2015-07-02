function p = testSignificanceOfHistF1oDC(binCent, binN_data, numVals, Fn)
%     nVals1 = length(binN_data);
%     nControl = length(vals_control);
%     assert(all(ibetween([binN_data(:); vals_control(:)], 0, 180)));

    nBoot = 1000;    
    
    %{
    For the dF1 distribution, calculating the cosine component of the
    F1/DC is equivalent to calculating the F2/DC of the reflected    
    distribution, because everything is distributed.
    
    For the dPhi distribution, they are not the same, because the 0 and 180 
    are not reflected.     
    
    %}
    reflectFlashedDphi = 1;
    isFlashedDphi = binCent(1) == 0;
    
    returnVal = 'F1oDC';
%     returnVal = 'logPval';
    tail = 'both';
    nbins = length(binCent);
    
    if ~isFlashedDphi  % ie. either dF1, or drifting grating dPhi.
                    
        idx = [1:nbins];%, nbins:-1:1];
        binC_refl = binCent;
        vals1_refl = binN_data;
        T = 180;
            
        F1oDC_data = getF1oDC(binC_refl, vals1_refl, T, Fn, 'cos');
            
        [control_mean, control_std] = getNullF1oDCdistribution(nbins, numVals, Fn, 'cos', false);

%         p = getValueProbFromDist(F1oDC_data, control_mean, control_std, 'both');
%         p = F1oDC_data;
            
    elseif isFlashedDphi  % flashed grating  dphi distribution
        
        if reflectFlashedDphi

            binW = diff(binCent(1:2));
            idx2 = [1:nbins, nbins-1:-1:2];
            binC_refl2 = binW*([0:2*nbins-3]);
            vals1_refl2 = binN_data(idx2);
            T = 360;                        
            F1oDC_data = getF1oDC(binC_refl2, vals1_refl2, T, 2*Fn, 'F1');            
            [control_mean, control_std] = getNullF1oDCdistribution(nbins, numVals, Fn, 'F1', true);
%             [control_mean_2, control_std_2] = calcNullF1oDCdistribution(nbins, numVals, Fn, 'F1', true);
%             p = F1oDC_data;
        else
            
            T = 180;            
            F1oDC_data = getF1oDC(binCent, binN_data, T, 1, 'cos');            
            [control_mean, control_std] = getNullF1oDCdistribution(nbins, numVals, Fn, 'cos', false);
%             [control_mean2, control_std2] = calcNullF1oDCdistribution(nbins, numVals, Fn, 'cos', false);
%             p = F1oDC_data;
        end
            
            
    end

    
    switch returnVal
        case 'F1oDC', p = F1oDC_data;
        case 'pval',  p = getValueProbFromDist(F1oDC_data, control_mean, control_std, tail);
        case 'logPval', pval = getValueProbFromDist(F1oDC_data, control_mean, control_std, tail);
            if pval == 0
                p = 99;
            else
                p = -log10(pval);
                if p > 20
                    p = round(p);
                end
            end
            
    end
    if ~(imag(p) == 0)
        3;
    end
        
end




%             binW = diff(binCent(1:2));
%             idx2 = [1:nbins, nbins:-1:1];
%             binC_refl2 = binW*([1:2*nbins]);
%             vals1_refl2 = binN_data(idx2);
%             T2 = 360;                        
%             F1oDC_data2 = getF1oDC(binC_refl2, vals1_refl2, T2, 2*Fn, 'F1');
            
%             assert( abs(F1oDC_data2 - abs(F1oDC_data)) < 1e-3 );            
%             F1oDC_control = zeros(1,nBoot);
%                         
%             for i = 1:nBoot
%                 randBinSamples = randi(nbins, 1, numVals);
%                 [uVals_tmp, uBinCount] = unique(randBinSamples);
%                 binCount = zeros(1,nbins);
%                 binCount(uVals_tmp) = uBinCount;
%                 valsControl_refl = binCount(idx);
%                 f1odc = getF1oDC(binC_refl, valsControl_refl, T, Fn, 'cos');
%                 F1oDC_control(i) = f1odc;        
%             end
                        
%             m = mean(F1oDC_control);
%             s = std(F1oDC_control);
%             [m-s, m+s, F1oDC_data]
%             p = F1oDC_data;
%             p = getValueProbFromDist(F1oDC_data, F1oDC_control, 'both');
%             p = getValueProbFromDist(F1oDC_data, beta(1), beta(2));



%             assert( abs(F1oDC_data2 - abs(F1oDC_data)) < 1e-3 );            

%         for i = 1:nBoot
%             randBinSamples = randi(nbins, 1, numVals);
%             [uVals_tmp, uBinCount] = unique(randBinSamples);
%             binCount = zeros(1,nbins);
%             binCount(uVals_tmp) = uBinCount;
%             valsControl_refl = binCount(idx);
%             f1odc = getF1oDC(binC_refl, valsControl_refl, T, Fn, 'cos');
%             F1oDC_control(i) = f1odc;        
%         end
            
%         m = mean(F1oDC_control);
%         s = std(F1oDC_control);
%             [m-s, m+s, F1oDC_data]
%             p = F1oDC_data;
%         p = getValueProbFromDist(F1oDC_data, control_mean, control_std);
%             p = getValueProbFromDist(F1oDC_data, beta(1), beta(2));

            
%             if discreteValues
%                 uVals = unique([binN_data(:); vals_control(:)]);
%                 [deltaPhis, deltaPhiNullCount] = dphi_details{:};
%                 assert( all(  arrayfun(@(x) any(x == deltaPhis), uVals) ));
%                 if length(uVals) <= 5  % flashed gratings
%                     nbins = 5;
%                     binCenters = linspace(0, 180, nbins);    
%                     binEdges = binCent2edge(binCenters);    
%                 else
%                     nbins = 11;
%                     binCenters = linspace(0, 180, nbins);    
%                     binEdges = binCent2edge(binCenters);    
%                 end
%             
%             else   
%                 nbins = 11;
%                 binEdges = linspace(0, 180, nbins+1);    
%                 binCent = binEdge2cent(binEdges);
%                 binEdges(1) = -.1; binEdges(end) = 180+.1;
%             end