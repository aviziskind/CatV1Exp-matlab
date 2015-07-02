function [bestOSPs, bestOSPinds] = findBestOriSpPhs(OSP, varargin)
    % either call with 2D OSP (oris & sps):
    %       bestOSPs = findBestOriSpPhs(OSP, oris, sps)
    % or with 3D OSP (oris, sps, phs):
    %       bestOSPs = findBestOriSpPhs(OSP, oris, sps, phs)

    option = 2;
    
    % option 1:  look at the response to the top 5% of OSPs.
    if (option == 1)
        p = .15;
        n = round( numel(OSP) * p );
        [m, bestOSPinds] = maxElements(OSP, n);

    % option 2:  look at the response to the bins with above mean firing rate.
    %  -- that have firing rate over 1 std above mean.
    elseif (option == 2)
        nstd = 1;
        bestOSPinds = findV(OSP > mean(OSP(:)) + ( std(OSP(:))  * nstd ) );
        
    end
        
    bestOSPs = zeros(size(bestOSPinds));
    for i = 1:size(bestOSPinds,2)
        bestOSPs(:,i) = varargin{i}(bestOSPinds(:,i));
    end
%     bestOSPs(:,1) = oris( inds(:,1) );
%     bestOSPs(:,2) = sps( inds(:,2) );
%     bestOSPs(:,3) = phs ( inds(:,3) );

end