function [dPhi, cc_atMaxDphi, phi1, phi2] = deltaPhi(phis, y1, y2, dPhiMode, nonCircFlag)

% when the first 2 inputs are scalars 
%     dPhi = deltaPhi(ph1, ph2, [phimax], dPhiMode)
%     	-- phi1 and phi2 are two phases, and dPhi is the shortest distance
%     between them, (mod phimax), 

% when the inputs are vectors, this function is defined as
%     dPhi = deltaPhi(phis, y1, y2, dPhiMode)
%     where phis is the set of phases, and y1 and y2 are sets of data that
%     correspond to the phase values. dPhi is deltaPhi between (a) the peaks of
%     the two datasets, or (b) the phases from computing the angle of the
%     f1 components of each data set.
%
% dPhiMode can be any of the following:
% 	  'dist between maxes'
%     'cross-correlation'
%     'F1 phase'
%     'dist between COMs'
%
%
% NOTE: ANGLES MUST BE IN DEGREES
     showWorking = false; 
     cc_atMaxDphi = 0;
     if ~exist('dPhiMode', 'var') || isempty(dPhiMode)
         dPhiMode = 'F1phase';
     end
     okModes = {'dist between maxes', 'cross-correlation', 'F1 phase', 'dist between COMs'};
     needCirc = {[], 1, 1, []};
     mode_id = find(strcmp(dPhiMode, okModes),1);     
     if isempty(mode_id)
         error('Unknown deltaPhi mode');
     end
     
     circular = ~exist('nonCircFlag', 'var') || isempty(nonCircFlag) || ~nonCircFlag;
     if ~isempty(needCirc{mode_id}) && (needCirc{mode_id} ~= circular)
         error(['(non)circular deltaPhi not available for mode == ' okModes{mode_id}]);
     end
     
     if circular
         phimax = 360;
         distFunc = @(x1, x2) circDist(x1, x2, phimax);
     else
         distFunc = @(x1, x2) abs(x2-x1);
     end
     
    
    if any( [all(y1 == 0), all(y2 == 0) ])
        dPhi = NaN;
        return;
    end
    
    ds = abs(diff(phis,2));
    if any(ds > 1e-2) %
        error('phis are not evenly spaced');
    end
    
    phis = phis - min(phis);  % for some phase sets, don't start at 0
%     phimax = phis(end) + phis(2)-phis(1);  % eg. phis will be [0, 90, 180, 270]. then phimax = 360; (unless nph = 44).     
       
    switch dPhiMode
        case 'dist between maxes'
%                 idx_maxPhi1 = indmax( y1 );
%                 idx_maxPhi2 = indmax( y2 );
            phi1 = phis( indmax( y1 ) );
            phi2 = phis( indmax( y2 ) );            
            dPhi = distFunc(phi1, phi2);

        case 'cross-correlation'
    
            % create correlation matrix
            maxFrac = .99;
            [n, nSamples] = size(y1);
            if nSamples == 1
                            
                Y2 = shiftMtx(double(y2));
                Y1 = y1(:, ones(1, n));            
                dotprods = sum(Y1 .* Y2, 1);
                [max_dot, ind_maxdot] = max(dotprods);
                if nnz(dotprods >= max_dot*maxFrac) > 1   % two (or more) peaks of the same height.
                    dPhi = nan;
                    cc_atMaxDphi = nan;
                else                
                    dPhi = phis(ind_maxdot);
                    if dPhi > phimax/2
                        dPhi = phimax-dPhi; % == mod(dPhi, phimax/2)???
                    end
                    if nargout >= 2
                        cc_atMaxDphi = pearsonR(y1, Y2(:,ind_maxdot));                        
                    end
                end
            
            elseif nSamples > 1
                
                shift_idxs = reshape(shiftMtx(1:n), [n, 1, n]);
                shift_idx_3d = bsxfun(@plus, [0:nSamples-1]*n, shift_idxs) ;
                Y1 = y1(:,:, ones(1,n));
                Y2 = y2(shift_idx_3d);
                
                dotprods = sum(Y1 .* Y2, 1);
                [max_dot, ind_maxdot] = max(dotprods, [], 3);
                numEqToMax = sum( bsxfun(@ge, dotprods, max_dot*maxFrac), 3);
                
                idx_ok = numEqToMax == 1;
                dPhi = phis(ind_maxdot);
                dPhi(~idx_ok) = nan;                
                idx_overHalfMax = dPhi > phimax/2;
                dPhi(idx_overHalfMax) = phimax-dPhi(idx_overHalfMax);
                                                
                cc_atMaxDphi = dPhi; % initialize with the same nans;                
                y2_idx = bsxfun(@plus, bsxfun(@plus, [1:n]', [find(idx_ok)-1]*n),  (ind_maxdot(idx_ok)-1)*(n*nSamples));
                if ~isempty(y2_idx)
                    cc_atMaxDphi(idx_ok) = pearsonR_v(Y1(:,idx_ok, 1), Y2(y2_idx) );
                end
                
            end
        
        
        case 'F1 phase'
    %         phi1 = rad2deg( getF1phase( phis, y1, phimax ));
    %         phi2 = rad2deg( getF1phase( phis, y2, phimax ));
    %         dPhi = distFunc(phi1, phi2);
    %         phis_rad = deg2rad(phis);
    %         phimax_rad = deg2rad(phimax);
            
            y1 = y1/max(y1);
            y2 = y2/max(y2);
            
            phi1 = rad2deg(getF1phase( phis, y1, phimax ));
            phi2 = rad2deg(getF1phase( phis, y2, phimax ));        
            dPhi = distFunc(phi1, phi2);
    
            if showWorking
    %         if rand < .01;%all(y1>0) && all(y2 >0)
                [phi1_rad, f_cos1, t_f1] = getF1phase( phis, y1, phimax );
                [phi2_rad, f_cos2, t_f2] = getF1phase( phis, y2, phimax );        
                fprintf('%.1f, %.1f\n', rad2deg(phi1_rad), rad2deg(phi2_rad))
                f_cos1 = f_cos1/max(f_cos1)*max(y1);
                f_cos2 = f_cos2/max(f_cos2)*max(y2);
                figure(213); clf;
                phis_ext = [phis(:); phis(end)+diff(phis(1:2))]; y1_ext = [y1(:); y1(1)]; y2_ext = [y2(:); y2(1)];
                plot(phis_ext, y1_ext, 'bo-', phis_ext, y2_ext, 'go-', t_f1, f_cos1, 'b:', t_f2, f_cos2, 'g:');
                title(num2str(dPhi));
                3;
            end



        case 'dist between COMs'
            x = phis(:);
            mx = mean(x);
            y1 = y1(:)/sum(y1);
            y2 = y2(:)/sum(y2);
            
            phi1 = sum( x .* y1); % = com1
            phi2 = sum( x .* y2); % = com2           
            dPhi = distFunc(phi1, phi2);
            
            if showWorking
                figure(148); clf;
                extt = @(x) [x, x(end)];
                binE = binCent2edge(x);
                
                stairs(binE, extt(binV)); hold on; stairs(binE+.03, extt(binV2), 'g');
                xlim([binE(1), binE(end)]);
            end
            
    end

    
end


%                 doTest = 0;
%                 if doTest
%                     dPhi2 = zeros(1,nSamples);
%                     for i = 1:nSamples
%                         dPhi2(i) = deltaPhi(phis, y1(:,i), y2(:,i), dPhiMode);
%                     end
%                     assert(isequalwithequalnans(dPhi, dPhi2));
%                 end

