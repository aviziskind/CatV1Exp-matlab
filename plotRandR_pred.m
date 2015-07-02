function plotRandR_pred(Gid, cellId)

    persistent allCells
    
    if isempty(allCells)
        allCells = load('indivCells_GLFcuw8_movie_fg.mat');
%         allCells = S.allCells;
        
    end
    
    fn = getName('celldata', Gid, cellId);
%     idx = find([allCells.Gid] == Gid & [allCells.cellId] == cellId. 1);
    
    C = allCells.(fn);
    MID = C.OSP.MID;
    R = C.OSP.R;
    
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