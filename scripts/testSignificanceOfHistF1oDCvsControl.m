function p = testSignificanceOfHistF1oDCvsControl(vals1, vals_control, F, dphi_details)
    nVals1 = length(vals1);
    nControl = length(vals_control);
    assert(all(ibetween([vals1(:); vals_control(:)], 0, 180)));

    nBoot = 100;
    
    discreteValues = nargin >= 4;
    
    if discreteValues
        uVals = unique([vals1(:); vals_control(:)]);
        [deltaPhis, deltaPhiNullCount] = dphi_details{:};
        assert( all(  arrayfun(@(x) any(x == deltaPhis), uVals) ));
        if length(uVals) <= 5  % flashed gratings
            nbins = 5;
            binCenters = linspace(0, 180, nbins);    
            binEdges = binCent2edge(binCenters);    
        else
            nbins = 11;
            binCenters = linspace(0, 180, nbins);    
            binEdges = binCent2edge(binCenters);    
        end
            
    else   
        nbins = 11;
        binEdges = linspace(0, 180, nbins+1);    
        binCent = binEdge2cent(binEdges);
        binEdges(1) = -.1; binEdges(end) = 180+.1;
    end
        
    N_data = histcnt(vals1, binEdges);
    if discreteValues
        N_data = N_data ./ deltaPhiNullCount;
    end
    
    
    if F == 1
        getFnoDC_func = @getF1oDC;
    else
        getFnoDC_func = @getFxoDC;
    end
        
    F1oDC_data = getFnoDC_func(binCent, N_data, 180);
    F1oDC_control = zeros(1,nBoot);
    
    for i = 1:nBoot
        r = randi(nControl, 1, nVals1);        
        N_control = histcnt(  vals_control(r), binEdges  );
        if discreteValues
            N_control = N_control ./ deltaPhiNullCount;
        end
        f1odc = getFnoDC_func(binCent, N_control, 180);
        F1oDC_control(i) = f1odc;        
    end
    figure(523); 
    m = mean(F1oDC_control);
    s = std(F1oDC_control);
    [m-s, m+s, F1oDC_data]
    
    p = getValueProbFromDist(F1oDC_data, F1oDC_control);
    
%     hist(F1oDC_control); 

%     DC = fourierTransform(0, dts, t, f, T);
%     3;
%     mean(DC./DC2)
    % F1 component
%     fourierTransform2 = @(w)   (1/sqrt(2))* sum( dts .* exp(1i * w * t) .* f, 2) /T;
% 
%     f1 = fourierTransform2(1);
%     f_1 = fourierTransform2(-1);
    
%     f1 = fourierTransform(1, dts, t, f, T);
%     f_1 = fourierTransform(-1, dts, t, f, T);


end


function y = fourierTransform(w, dts, t, f, T)
    % 1/sqrt(2) = 0.7071;
    y =  0.7071 * sum( dts .* exp(1i * w * t) .* f, 2) /T;
end
