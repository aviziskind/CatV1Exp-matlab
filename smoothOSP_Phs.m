function OSP = smoothOSP_Phs(OSP, action, n)
    % action must be "gauss", "fermi", or "alias"'

    %input: OSP (with R and/or R_full fields), or R, or R_full, or single
    %vector
    OSP_orig = OSP;
    if isvector(OSP_orig)        
        OSP = reshape(OSP, [1 1 length(OSP)]);
    end
    %%
    % parse input    
    if isstruct(OSP)
        R = OSP.R;
    else
        R = OSP;
    end
    if size(R,4) > 1
        doR_full = true;
        doR = false;
        R_full = R;
    else
        doR = true;   
        doR_full = isstruct(OSP) && isfield(OSP, 'R_full');
        if doR_full
            R_full = decompress( OSP.R_full );
        end
    end
    
    % get size of array(s)
    if doR
        [nOri, nSpf, nPh] = size(R);    
    end    
    if doR_full        
        [nOri, nSpf, nPh, nTrialsMax] = size(R_full);
    end
    
    % possible nPhases for drifting gratings: nPh = 40, 44, 60, 120
    % There are 2 ways of "compressing" the number of phases:
    % 1) using a smoothing kernel to smooth the phase tuning curves ("smooth") 
    % 2) averaging over a few phases at a time. ("average")
    
    if nargin < 2
        action = 'gauss';
    end
    action = lower(action);
    
    if any(strcmp(action, {'fermi'}));
        actionMode = 'indiv';
        if nargin >= 3
            n_smooth = n;
        else
            n_smooth = 2;
        end
        
    elseif any(strcmp(action, {'gauss', 'alias'}))
        actionMode = 'together';
        if nargin >= 3
            n_smooth = n;
        else
            n_smooth = nPh/2;
        end        
    else
        error('Unknown action type: must be "gauss", "fermi", or "alias"');
    end
        
    switch action
        case 'gauss', smoothFunc = @(r) gaussSmooth(r, n_smooth, 3, 'circ');
            nPh_cmp = nPh;
        case 'fermi',
            if length(n_smooth) == 1,  n_smooth(2) = n_smooth(1)/50;   end                
            smoothFunc = @(r) fermiSmooth(r, n_smooth(1), n_smooth(2), 'circ');
            nPh_cmp = nPh;
        case 'alias', smoothFunc = @(R) alias(R, n_smooth, 3);
            nPh_cmp = n;
    end
            
                        
    switch actionMode
        case 'indiv'
            if doR
                R_cmp = zeros(nOri, nSpf, nPh_cmp);
                for i = 1:nOri
                    for j = 1:nSpf                        
                        R_cmp(i,j,:) = smoothFunc(R(i,j,:));
                    end  
                end
            end
            if doR_full
                R_full_cmp = zeros(nOri, nSpf, nPh_cmp, nTrialsMax);

                for i = 1:nOri
                    for j = 1:nSpf                        
                        for k = 1:nTrialsMax
                            R_full_cmp(i,j,:,k) = smoothFunc(R_full(i,j,:,k));
                        end
                    end  
                end

            end
            
        case 'together'
            if doR
                R_cmp = smoothFunc(R);
            end
            if doR_full
%                 R_full_cmp = zeros(nOri, nSpf, nPh_cmp, nTrialsMax);
%                 for k = 1:nTrialsMax
%                     R_full_cmp(:,:,:,k) = smoothFunc(R_full(:,:,:,k));
%                 end
                R_full_cmp = smoothFunc(R_full);
            end
            
            
    end
        
    
    % return phase-averaged result
    if isstruct(OSP)
        OSP.R = R_cmp;
        if doR_full
            OSP.R_full = R_full_cmp;
        end
        ph = linspace(0,360, nPh_cmp+1); ph = ph(1:nPh_cmp);
        OSP.ph = ph;        
    else
        if doR
            OSP = R_cmp;
            if isvector(OSP_orig)
                OSP = squeeze(OSP)';
            end
            
        elseif doR_full
            OSP = R_full_cmp;
        end
    end

end


    
%{
        case 'alias'
            if nargin < 3
% %                 switch nPh
% %             %         case {40, 60, 120}, nPh_cmp = 10;
% %         %             case 40, nPh_cmp = 10;
% %         %             case 60, nPh_cmp = 12;
% %         %             case 120, nPh_cmp = 15;
% %         %             case 44, nPh_cmp = 11;
% % 
% %                     case 40, nPh_cmp = 8;
% %                     case 60, nPh_cmp = 6;
% %                     case 120, nPh_cmp = 8;
% %                     case 44, nPh_cmp = 4;
% %                     otherwise
% %                         return;
% %                 end    
%             else
%                 nPh_cmp = n;
%             end
%             nFramesToAverage = nPh/nPh_cmp;                            
                nFramesToAverage = 2;
            else
                nFramesToAverage = n;
                if (n ~= round(n))
                    error('n must be a whole number');
                end
            end
            nPh_cmp = nPh/nFramesToAverage;
                        
            % initialize
            if doR
                R_cmp = zeros(nOri, nSpf, nPh_cmp);
            end
            if doR_full
                R_full_cmp = zeros(nOri, nSpf, nPh_cmp, nTrialsMax);
            end

            % partially alias over phases
            for idx_av = 1:nPh_cmp
                idx_full = (idx_av-1)*nFramesToAverage + 1 : min(idx_av*nFramesToAverage, nPh);
                if doR
                    R_cmp(:,:,idx_av) = nanmean(R(:,:,idx_full), 3);
                end
                if doR_full
                    R_full_cmp(:,:,idx_av,:) = nanmean(R_full(:,:,idx_full,:), 3);
                end
            end

%}