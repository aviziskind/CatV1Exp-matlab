function [R_pred, nlinFunc, cc] = predictOSPfromMID(Gid, MID, R, figid1, OSP)
    
%     P_s_v
%     P_s_v_spk    
    returnJustSqr = 0;

    nlinType = 'powerlaw';
    nlinType = 'hist';

    [frameStimIds, uOri, uSpf, uPh] = getStimulusFrameSequence(Gid);            
    [nOri, nSpf, nPh] = size(R);
    
    [uStimIds, idx_firstOccurence] = unique(frameStimIds, 'first');    
    nStim = length(uStimIds);
    stim_firstFrameIds = idx_firstOccurence; %reshape(idx_firstOccurence);    
        
    
    getFrame = getFrameRetrieverFunction(Gid);    
    getFrame('load', Gid, 'scaling', 'aroundZero');    
    frameDims = getFrame('size');
    dsamp = frameDims(1) / length(MID);
        
%     R_pred = zeros(nOri, length(uSpf), length(uPh));
    s_v = zeros(1,nStim);
    s_v_spk = zeros(1,nStim);
    MID = double(MID);
    R = R/max(R(:));

    allFrames = zeros([frameDims, nStim], 'single');

    % gather all frames;
    for stim_i = 1:nStim
        frm_id = stim_firstFrameIds(stim_i);         
        frm_i = double( getFrame(frm_id) );
        allFrames(:,:,stim_i)  = frm_i;
            
    end    
    
    getFrame('close');
    if dsamp > 1
        allFrames = downSample(allFrames, dsamp);
    end
    3;
    
%         % compute probabilities & non-linearities
    for stim_i = 1:nStim
        frm_i = allFrames(:,:,stim_i);
        s_v(stim_i)     = sum( MID(:) .* frm_i(:) );
        s_v_spk(stim_i) = sum( MID(:) .* frm_i(:) ) * R(stim_i);        
    end
    
    % POWER-LAW fit
%     plaw = @(b, x) b(1).*rectified(x).^(abs(b(2)));

%     b0 = [1 2];    
%     plaw_coef = nlinfit(s_v, s_v_spk, plaw, b0);
%     plaw_coef(2) = abs(plaw_coef(2));
%     % 	if plaw_coef(2) < 1
%     % 		plaw_coef(2) = 2;
%     % 	end		    
%     R_pred_plaw = reshape( plaw(plaw_coef, s_v ), size(R));    
    

    if returnJustSqr
        R_pred = reshape( s_v.^2, size(R));    
        return;        
    end

    % Histogram fit. (binning method... (?))
    L = lims([s_v(:); s_v_spk(:)]);
    nbins = 20;
    binE = linspace(L(1), L(2), nbins);
    binC = binEdge2cent(binE);

    binVal_s_v = histcnt(s_v, binE);
    binVal_s_v_spk = histcnt(s_v, binE, R);
%     binVal_s_v_spk = histcnt(s_v_spk, binE);
    binSum = (binVal_s_v + binVal_s_v_spk);

    P_s_v     = binVal_s_v    /sum(binVal_s_v);
    P_s_v_spk = binVal_s_v_spk/sum(binVal_s_v_spk);
    P_spk     = binSum ./ sum(binSum);
      
    P_spike_s_v = P_s_v_spk./(P_s_v+eps);
    assert(~any(isnan(P_spike_s_v)));
    P_spike_s_v_orig = P_spike_s_v/sum(P_spike_s_v);
    P_spike_s_v = gaussSmooth(P_spike_s_v_orig, 1);
    
    binC_ext = extend_edges(binC);
    P_spike_s_v_ext = extend_edges(P_spike_s_v);
    
    R_pred_bin = interp1(binC_ext, P_spike_s_v_ext, s_v, 'spline');    
    
    R_pred = reshape(R_pred_bin, size(R));
    
    nlinFunc = [binC_ext(:)'; P_spike_s_v_ext(:)'];
    cc = corr(R(:), R_pred(:));
    3;
    


    show = exist('figid1', 'var') && ~isempty(figid1);
    
    if show
        if nargin < 4
            figid1 = 90;
        end
        sd = siteDataFor(Gid);
        [xs, ys] = getStimulusXY(sd.stimulusInfo, dsamp);
        figure(figid1+1); imagesc(xs, ys, MID); axis square;      
        title('Most Informative Dimension')
        
%         figure(figid1+2); clf;
%         plot(s_v, s_v_spk, 'o'); hold on;
%         xlim(lims(s_v, .05));
%         fplot(@(x) plaw(plaw_coef, x), xlim, 'r');        
%         title(sprintf('y = %.2f x^{%.2f}', plaw_coef));
%         xlabel('s . v'); ylabel('( s . v ) x P(spike)')

        
        %{
            plot(binC, P_s_v, 'o-'); xlabel('s.v'); ylabel('P(s.v)');  xlim(lims(s_v, .05));
            plot(binC, P_s_v_spk, 'o-'); xlabel('s.v'); ylabel('P(s.v|spike)');  xlim(lims(s_v, .05));
        
        %}
        
        figure(figid1+3); clf;
        binC_itp = linspace(binC_ext(1), binC_ext(end), 200);
        P_spike_s_v_itp = interp1(binC_ext, P_spike_s_v_ext, binC_itp, 'spline');            
        
        plot(binC, P_spike_s_v_orig, 'ro'); hold on;
        plot(binC, P_spike_s_v, '.'); hold on;
        hh = plot(binC_itp, P_spike_s_v_itp, 'b-', 'markersize', 1);                
        xlabel('s.v'); ylabel('P(spike|s.v)');        
        xlim(lims(binC_ext, .05));        
        legend('original', 'smoothed', 'location', 'NW')
        title('Sampled input-output function')
        
        
        haveOSP = nargin > 4;
        figure(figid1+4);
        subplot(1,2,1);
        if haveOSP
            imageOSP(OSP, 'mean:ph', 'OSP', 'noLabels');
        else
            imagesc(mean(R,3));
        end
        title('R (observed)');

        subplot(1,2,2);
        if haveOSP
            OSP.R = R_pred;
            imageOSP(OSP, 'mean:ph', 'OSP', 'noLabels');
        else
            imagesc(mean(R_pred,3));
        end
        title('R (predicted)');
        
        drawnow;
        refresh;
        
        if figid1 > 85
            3;
        end
        
%         figure(82); clf;
%         plot(P_s_v, P_s_v_spk, 'o-')
% 
%         figure(83);
% %         plot(binC, P_s_v_spk .* P_spk ./ P_s_v, 'x-'); hold on;    
%         plot(binC, P_s_v_spk ./ (P_s_v.*P_spk), 'go-');                
    end


end


function y_ext = extend_edges(y)
    dy_1 = y(2)-y(1);
    dy_end = y(end)-y(end-1);
    
    y_ext = [y(1)-dy_1; y(:); y(end)+dy_end];
end

%{

    % gather all frames;
    for stim_i = 1:nStim
        frm_id = stim_firstFrameIds(stim_i);         
        frm_i = double( getFrame(frm_id) );
        if dsamp > 1
            frm_i = downSample(frm_i, dsamp);
        end
%         allFrames(:,:,stim_i)  =frm;
    
        % compute probabilities & non-linearities
%         frm_i = allFrames(:,:,stim_i);
        s_v(stim_i) = sum( MID(:) .* frm_i(:) );
        s_v_spk(stim_i) = R(stim_i) .* sum(MID(:) .* frm_i(:));        
        
    end    
    getFrame('close');

%}