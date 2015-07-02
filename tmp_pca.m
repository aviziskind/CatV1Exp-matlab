% 
% 
% figure(1);
% 
% 
% 
% spk = dbGetSpikes(4470, [], [], 1);
% tet = spk(:,3:6);
% v_orig = normV(tet,2);
% 
% dbSpk 

            [nT, nChannels, nSpk] = size(spk_wvfms_raw);
            spk_wvfm_cat = reshape(spk_wvfms_raw, [nT*nChannels, nSpk])';
            meanWaveform = mean(spk_wvfm_cat,1);
            spk_wvfm_cov = cov(spk_wvfm_cat); % covariance matrix
            [wvfm_eigvc, wvfm_ev] = eig(spk_wvfm_cov);
            wvfm_ev = diag(wvfm_ev);
            
            [tmp, idx_top_eigs] = sort(wvfm_ev, 'descend');
            
%             figure(343);
%             plot(
%                         
            nEigs = 4;
            pca_comps = wvfm_eigvc(:, idx_top_eigs(1:nEigs) );
            spk_wvfm_cat_msub = bsxfun(@minus, spk_wvfm_cat, meanWaveform);
            3;
            
            wvfm_eigvc' * e1
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            