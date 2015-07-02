function [dPhi, allShiftCCs] = bestCircShift_Matlab(x, y1, y2, returnNanForAmbiguous, allowNegativeValues)
    
%     global count
    if nargin < 4 || isempty(returnNanForAmbiguous)
        returnNanForAmbiguous = 1;
    end
    if nargin < 5 
        allowNegativeValues = 0;
    end
        
    if isvector(y1) && size(y1,2) > 1
        y1 = y1';
    end
    if isvector(y2) && size(y2,2) > 1
        y2 = y2';
    end

     x = x-min(x);  % for some phase sets, don't start at 0
     x_max = x(end)+diff(x(1:2));
%     phimax = phis(end) + phis(2)-phis(1);  % eg. phis will be [0, 90, 180, 270]. then phimax = 360; (unless nph = 44).     

    maxFrac = .9999;
    [nPhases, nTuningCurves] = size(y1);

     if any( [all(y1 == 0), all(y2 == 0) ])
         dPhi = NaN(1, nTuningCurves);
         return;
     end            
       

    doCC_quick = 1;
    if nTuningCurves == 1
        
        if doCC_quick            
            y1_nrm = (y1-mean(y1)); y1_nrm = y1_nrm / norm(y1_nrm);
            y2_nrm = (y2-mean(y2)); y2_nrm = y2_nrm / norm(y2_nrm);
            y1 = y1_nrm; y2 = y2_nrm;            
        end
        assert( abs( dot(y1, y2) - pearsonR(y1, y2) ) < 1e-5);

        Y2 = shiftMtx(double(y2));
        Y1 = y1(:, ones(1, nPhases)); 
        dotprods = sum(Y1 .* Y2, 1);
        allShiftCCs = dotprods;
        [max_dot, ind_maxdot] = max(dotprods);
        tf_aboveTh = (dotprods > max_dot*maxFrac);
        if (nnz(tf_aboveTh) > 1) && returnNanForAmbiguous   % two (or more) peaks of the same height.
            % if just two adjacent similar dot products, allow, but use
            % half-even rounding
%             idx_aboveTh = find(tf_aboveTh);
%             if idx_aboveTh                
%             [sortedCCs, idx_srt] = sort(dotprods, 'descend');
%             assert(idx_srt(1) == ind_maxdot);
%             idx_aboveTh_srt = find(sortedCCs > max_dot*maxFrac);
            
            
            idx_aboveTh = find(tf_aboveTh);
%             if circDist(idx_srt(1), idx_srt(idx_aboveTh_srt), n) > 2

            dists = circDist(ind_maxdot(1), idx_aboveTh, nPhases);
            if any(dists >= 2) || (nPhases <= 4)
                dPhi = nan;
            else
                
                if ~( any(idx_aboveTh == 1) && any(idx_aboveTh == nPhases) )
                     dphi_idx = roundHalfEven( mean(idx_aboveTh) );
                else
                     dphi_idx = idx_aboveTh(1);
%                      keyboard;
                end
                dPhi = x( dphi_idx );                
            end
            dPhi = nan;
           
%             cc_atMaxDphi = nan;
%             if nnz(dotprods == max_dot) > 1
%                 count = count + 1;
%             end
% 
%             dPhi = x(ind_maxdot(1));
%             if dPhi > x_max/2
%                 dPhi = x_max-dPhi; % == mod(dPhi, phimax/2)???
%             end


        else
            dPhi = x(ind_maxdot);
        end
        
        if ~isnan(dPhi) && ( dPhi > x_max/2 )
            if ~allowNegativeValues
                dPhi = x_max-dPhi; % == mod(dPhi, phimax/2)???
            else
                dPhi = dPhi-x_max; % == mod(dPhi, phimax/2)???
            end
            
%             if nargout >= 2
%                 cc_atMaxDphi = pearsonR(y1, Y2(:,ind_maxdot));
%             end
        end

    elseif nTuningCurves > 1

        shift_idxs = reshape(shiftMtx(1:nPhases), [nPhases, 1, nPhases]);
        shift_idx_3d = bsxfun(@plus, [0:nTuningCurves-1]*nPhases, shift_idxs) ;
        Y1 = y1(:,:, ones(1,nPhases));
        Y2 = y2(shift_idx_3d);

        dotprods = sum(Y1 .* Y2, 1);
        [max_dot, ind_maxdot] = max(dotprods, [], 3);        
        
        dPhi = x(ind_maxdot);
        
        if returnNanForAmbiguous
            numEqToMax = sum( bsxfun(@ge, dotprods, reshape(max_dot, [1, 1, nPhases])*maxFrac), 3);        
            idx_ok = numEqToMax == 1;
            dPhi(~idx_ok) = nan;
        end
        
        idx_overHalfMax = dPhi > x_max/2;
        if ~allowNegativeValues
            dPhi(idx_overHalfMax) = x_max-dPhi(idx_overHalfMax);
        else
            dPhi(idx_overHalfMax) = dPhi(idx_overHalfMax)-x_max;
        end

%         cc_atMaxDphi = dPhi; % initialize with the same nans;
%         y2_idx = bsxfun(@plus, bsxfun(@plus, [1:n]', [find(idx_ok)-1]*n),  (ind_maxdot(idx_ok)-1)*(n*nSamples));
%         if ~isempty(y2_idx)
%             cc_atMaxDphi(idx_ok) = pearsonR_v(Y1(:,idx_ok, 1), Y2(y2_idx) );
%         end

    end
        

    
end


function y = roundHalfEven(x)
    if abs(x - floor(x) - 0.5) < 1e-5
        if ~odd(floor(x))
            y = floor(x);
        else
            y = ceil(x);
        end
    else        
        y = round(x);
    end
    
end

%                 doTest = 0;
%                 if doTest
%                     dPhi2 = zeros(1,nSamples);
%                     for i = 1:nSamples
%                         dPhi2(i) = deltaPhi(x, y1(:,i), y2(:,i), dPhiMode);
%                     end
%                     assert(isequalwithequalnans(dPhi, dPhi2));
%                 end

