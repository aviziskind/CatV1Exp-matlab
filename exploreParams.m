function exploreParams(gratingType_s, id)
    if nargin < 1
        [gratingType, gratingType_s] = curGratingType;  
    end
    if nargin < 2
        id = 1;
    end

    [pt_ids, pairTypes] = curPairTypes;
    nP = length(pairTypes);
    pairTypes_str = [pairTypes{:}];
    
    paramSearchFilename = [gratingType_s 'ParameterSearch' num2str(id) '_tuning_' pairTypes_str '.mat'];
    S = load(paramSearchFilename);
    [X, Y, Z, measures, data, data_std, NumPairs, pvals, pvalLabels] = deal(...
        S.X, S.Y, S.Z, S.measures_h, S.data, S.data_std, S.numPairs, S.pvals, S.pvalLabels_h );
    S.negLogPvals = cellfun(@(x) -log10(x), pvals, 'un', 0);
    nM = length(measures);
    xs = X.vals;
    ys = Y.vals;

    tmp = zeros(length(xs), length(ys));
    nColMp = 100;
    jet1 = [1 1 1; jet(nColMp)];
    size_vec = [length(xs), length(ys)];
    
    for pr_i = 1:nP
        figure(pr_i); clf;
        set(pr_i', 'name', [gratingType_s ' : ' pairTypes{pr_i}]);
        for ms_i = 1:nM    
            hIm_ax(pr_i,ms_i) = subplot(1,nM+1,ms_i); %#ok<*AGROW>
            hIm(pr_i,ms_i) = imagesc(xs, ys, tmp);    
%             xlabel(X.thName, 'interpreter', 'none'); ylabel(Y.thName, 'interpreter', 'none');
            h_tit(pr_i, ms_i) = title(measures{ms_i});
            axis xy;
            colormap(jet1);
            hCol(pr_i,ms_i) = colorbar;
        end
        hIm_ax(pr_i,nM+1) = subplot(1,nM+1,nM+1);
        hIm(pr_i,nM+1) = imagesc(xs, ys, tmp);    
        xlabel(X.thName, 'interpreter', 'none'); ylabel(Y.thName, 'interpreter', 'none');
        title('log_{10} NumPairs');
        axis xy;
        colormap(jet1);
        hCol(pr_i,ms_i+1) = colorbar;        
    end
        
    
    if ~isempty(Z)
        if iscell(Z.vals)        
            zvals = Z.vals{1};
        else
            zvals = Z.vals;
        end
    else
        zvals = [];
    end
        
    function updateImages(show, z_val, pvalLim_in, N_th, varargin)
             

        pvalNamesToUse = reshape(varargin, nP, nM);
                        
        if isempty(Z)
            z_id = 1;
        else
            [tmp, z_id] = min(abs(zvals - z_val));                    
        end     
        
        
        
        for pr_j = 1:nP
            
            for ms_j = 1:nM
               
                if any(strcmp(show, {'pvals', 'negLogPvals'})) || ~strcmp(pvalLim_in, 'none')   % get p-values if needed.
                    pval_id = find(strcmp(pvalLabels{pr_j, 1, ms_j}, pvalNamesToUse{pr_j, ms_j}), 1);                    
                    pvals_here = S.pvals{pr_j,1,ms_j}(:,:,z_id, pval_id)';
                end
                if any(strcmp(show, {'pvals', 'negLogPvals'})) 
                    dim4_id = pval_id;
                else
                    dim4_id = 1;
                end
                    
                data_p_m = S.(show){pr_j,1,ms_j}(:,:,z_id, dim4_id)';  % (if not pval, then pval_id == 1 doesn't do anything).

                idxToBlank = false( size_vec );
                
                if ~strcmp(pvalLim_in, 'none')  % white-out values above a certain p-value threshold
                    pvalLim = eval(pvalLim_in);                    
                    idxToBlank = idxToBlank | (pvals_here > pvalLim);
                end
                if N_th > 0                      % white-out values with N below a certain threshold
                    numPairs_here = NumPairs{pr_j,1,ms_j}(:,:,z_id)';
                    idxToBlank = idxToBlank | (numPairs_here < N_th);
                end           
                                
                data_p_m(idxToBlank) = nan;
                set(hIm(pr_j,ms_j), 'cdata', data_p_m);        
                
                tit_str = measures{ms_j};
                if any(strcmp(show, {'pvals', 'negLogPvals'}))
                    tit_str = [tit_str ' : ' pvalNamesToUse{pr_j, ms_j}];
                end
                set(h_tit(pr_j, ms_j), 'string', tit_str)
            
                set(hIm_ax(pr_j,ms_j), 'climmode', 'auto');
                nullVal = iff(any(strcmp(measures{ms_j}, {'cc', 'rho'})), 0, 90);                                        
                
                clims = get(hIm_ax(pr_j,ms_j), 'clim');
                if strcmp(show, 'data')
                    clims = (1 + 1/nColMp)*[-1, 1]*max(abs(clims-nullVal))+nullVal; % make symmetric around null-value; expand a little for colormap
                else
                    d = diff(clims)/(nColMp-1);
                    clims = clims + [-1, 1]*d; % expand a little for colormap
                end
                set(hIm_ax(pr_j,ms_j), 'clim', clims);                                    
            end

            NPairs = NumPairs{pr_j,1,1}(:,:,z_id)'; % use NPairs from the first plot.
%             NPairs(idxToBlank) = nan;
            set(hIm(pr_j,nM+1), 'cdata', log10( NPairs ) );            
            
        end
            
    end

    showOptions = {'data', 'data_std', 'numPairs', 'pvals', 'negLogPvals'};
    if ~isempty(Z)        
        zval0 = zvals(1);
        z_arg = { Z.name, zvals, zval0 };
    else
        z_arg = {'z', 1};
        zval0 = 1;
    end
    
    pvalLims = {'none', '0.05', '0.01', '0.005', '0.001'};
    [p_k, m_k] = meshgrid(1:nP, 1:nM);
    plotNames = arrayfun(@(p,m) [pairTypes{p} '_' measures{m}], p_k, m_k, 'un', 0)';
    pvalNmOptions = pvalLabels(:);
    pvalNmOptions0 = cellfun(@(C) C{1}, pvalNmOptions, 'un', 0);
    pvalNameArgs = cellfun(@(pn, pn_op) {pn, pn_op}, plotNames(:), pvalNmOptions(:), 'un', 0 );
    
    
%     pval
    
    3;
    N_th0 = 0;
    args = { {'show', showOptions}, ...
             z_arg, ...
             {'pval_th', pvalLims}, ...
             {'N_th', [0:30], N_th0}, pvalNameArgs{:}, ...
           };
    
    updateImages(showOptions{1}, zval0, pvalLims(1), N_th0 ,pvalNmOptions0{:});
    manipulate(@updateImages, args, 'FigId', 100);



end



    %{
    if doSurf
        figure(1); clf;
%         h3D_ax = subplot(1,nTot,1);    
        h3D = surf(xs, ys, tmp); 
        xlabel(X.thName, 'interpreter', 'none'); ylabel(Y.thName, 'interpreter', 'none');
    end    


        if doSurf
            set(h3D, 'zdata', Data{1,1,1}(:,:,z_id)');
        end

    %}
