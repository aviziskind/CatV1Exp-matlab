function w = findHowWellTunedOSP(OSP, params, show)

    [window_size, pow] = elements(params);

%     tunedMethod = 'cross-covariance';
    tunedMethod = 'pval';
    
%     showOSP = false;
    
%     function D = distMatrixWrapX(nx, ny)
%         [xs, ys] = meshgrid(1:nx, 1:ny);
% 
%         distWrapX = @(x1, x2)  min(abs(x1-x2), abs(nx-(x1-x2)) );
%         distReg   = @(a1, a2)  abs(a1-a2);        
%         
%         dxM = crossOp(xs(:), distWrapX, xs(:)') ;
%         dyM = crossOp(ys(:), distReg,  ys(:)') ;
% %         D = sum( abs( stackMatrices(dxM, dyM)), 3 );
%         D = magsqr( stackMatrices(dxM, dyM), 3 ) + eps;
%     end

    % current algorithm: sum over phases
    if isstruct(OSP)
        R = OSP.R;
    else
        R = OSP;
    end
    OS = sum(R,3);
    [nOri, nSp] = size(OS);
    pval = .01;
    
    switch tunedMethod
        case 'energy'
    
            C = OS(:) * OS(:)';

            diagInds = (1:nOri*nSp).^2;
            C(diagInds) = 0;

            C = C/max(C(:));
            D = distMatrixWrapX(nOri, nSp);
            w = log( sum( C(:) ./ D(:) ) );

        case 'smoothness'
            ...
                    
            mu = mean(OS(:));

            s = window_size; n = s*2-1;
            window_uniform = ones(1,n);
            window_gaussian = gaussian( [-s:s], 0, s/2);

            window = window_gaussian;    
            window_2d = window(:) * window(:)';
            window_2d = window_2d / sum(window_2d(:));

            OS_smoothed = conv2_circ(OS, window_2d, 'same', [1]); % CHECK that: orientation is dimension 1


            diffs = abs(OS_smoothed - OS);
            w = sum( (diffs(:)/mu).^(pow) );

            if (nargin > 2) && (show)
                figure(30);
                X = [OS(:); OS_smoothed(:); diffs(:)];
                cax = [min(X), max(X)];
                subplot(2,2,1); imagesc(OS); caxis(cax);
                subplot(2,2,2); imagesc(OS_smoothed); caxis(cax); colorbar;
                subplot(2,2,3); imagesc(diffs); caxis(cax); 
                subplot(2,2,4); imagesc(diffs/mu); caxis(cax); colorbar
                3;
            end
        
        case 'pval',
            st = OSP.stats;
%             ori_sel_str = st.orientationSelectivityStrength;
            ori_sel_pval = st.orientationSelectivePval;
%             spf_sel_p = st.spatialfrequencySelectivePval;
%             alpha = .01;
            ori_rep_p = st.orientationReproduciblePval;
            ori_rep_str = st.orientationReproducibleStrength;
            spf_rep_p = st.spatFreqReproduciblePval;
            spf_rep_str = st.spatFreqReproducibleStrength;
%             w = -log(ori_sel_p);

%             base_ori = -log10(ori_rep_p);
%             w = -log10(ori_sel_pval);

            w = ori_rep_str + spf_rep_str;
            %             if (ori_rep_p < alpha) && (spf_rep_p < alpha) 
%                 w = ori_rep_str + spf_rep_str;
%             else
%                 w = 0;
%             end
%             w = ori_sel_str;
%             w = (ori_rep_p < pval) + (spf_rep_p < pval);

    end


end