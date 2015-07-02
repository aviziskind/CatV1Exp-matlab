function Xcenter = estimateGaborCenter(xs, ys, zs)
    [Ny, Nx] = elements(size(zs));
    threshold = 0.3;
    
%     zs = smooth2(zs, 3);
    % Often with gabors with many oscillations, there will be many peaks. 
    % So we findFinding all the peaks can we need to find the

    % rough estimate:
    [tmp, roughXc] = maxElement(abs(zs));
    
%     [theta] = findThetaLambdaWithStrongestOscillations(xs, ys, zs, X_center);
    
%     centerAreaX = [roughXc(1) - lambda*2
        
    
    [Xs, vals] = findPeaksAndTroughs(1:Nx, 1:Ny, zs - mean(zs(:)), threshold);
    vals = abs(vals);
    
    nPeaks = length(vals);
    if nPeaks > 2
        keyboard;
    end
    switch nPeaks
        case 1
            Xcenter = Xs;

        case 2
       
            meanval = mean(vals);
            Xcenter = [interp1(vals, Xs(1,:), meanval);
                       interp1(vals, Xs(2,:), meanval)];            
                 
        otherwise
            meanval = mean(vals);
            Xcenter = [interp1(vals, Xs(1,:), meanval);
                     interp1(vals, Xs(2,:), meanval)];            
            
    end
    Xcenter = round(Xcenter);
    
end