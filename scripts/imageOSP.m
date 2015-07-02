function [h_ax, h_im] = imageOSP(OSP, plotOptions, axOrder, varargin)
    % syntax:
    %    imageOSP(OSP, plotOptions, axOrder)
    %  ** if you want to preserve all 3 dimensions (ori/spf/phase), plotOptions should be of the form:
    %       '[subplots]/sameaxes:[horizontal]/vertical/diagonal' 
    %    otherwise you can indicate how you want to collapse one dimension,
    %    with:
    %       'mean/max/pref:ori/sp/ph'
    %    eg. 
    %       'mean:ori'  collapses the ori axes by taking the mean across orienations. 
    %       'pref:sp'   collapses the sp axes by picking out the preferred spatial frequency.
    %   Note: We assume that the dimensions of the input Response profile  are 
    %       always originally ordered by (1) Ori, then (2) SpFreq, then (3) Phase.
    %  ** axOrder should be the three letters O,S,P, in the order you wish them to be plotted.
    %       eg.  'OSP' or 'SPO', or 'OS'
    %       the axes will be (1) vertical, (2) horizontal, (3) 3rd dimension (specified by plotmode)
    %        you can leave out the 3rd letter if that dimension is
    %        collapsed somehow.
    % options for additional arguments: 
    %   * 'noLabels' - doesn't label the axes
    %   * 'noTicks' - doesn't put tickmarks on the axes;
    %   * h_ax / h_fig - update a plot that exists already
    fontsize = 12;
    
    if (nargin < 2) || isempty(plotOptions)
        plotOptions = 'subplots:horizontal'; % alternative: 'sameaxes' (put all data on one axes),  'mean:Phase/Ori/SpatFreq'
    end

    if (nargin < 3) || isempty(axOrder);
        axOrder = 'OSP';
    end   

    fd = 0; % flashed/drifting = unknown
    
    % Parse OSP input
    if isstruct(OSP)
        if isfield(OSP, 'OSP'),
            R = OSP.OSP;
        else
            R = OSP.R;
        end
        [oris, sp_pix, phases] = deal(OSP.ori, OSP.sp, OSP.ph);                
        if isfield(OSP, 'type')
            fd = iff(strcmp(OSP.type, 'flashed'), 1, 2);
        end
        %ori / direction
        ori_ticklabels = [];
        if ~isfield(OSP, 'tf_Hz') || (OSP.tf_Hz == 0) || (OSP.tf_Hz == inf) 
            ori_ticks = [0:45:180];
            ori_label = 'Orientation';
        else
            ori_label = 'Direction';
            if length(oris) > 10
                ori_ticks = [0:90:359];
            else
                ori_ticks = 1:length(oris); % just a few oris - label them manually
                ori_ticklabels = num2cell(oris);
            end
        end
        ori_forPlot = iff(~isempty(ori_ticklabels), ori_ticks, oris );
        
        % spatial frequency
        nSpf = length(sp_pix);
        if any(sp_pix == 0)
            sp_pix(sp_pix == 0) = .1;
        end
        if isfield(OSP, 'degPerPix')
            degPerPixel = OSP.degPerPix;            
        else
            degPerPixel = 1;
        end
        sp_degPerCycle = (sp_pix * degPerPixel); % (pixels/cycle)* (deg/pixel) = degPerCycle;            
        spf_cycPerDeg = 1./sp_degPerCycle;        
        spf_ticklabels = arrayfun(@(n) num2str(n, '%.2f'), spf_cycPerDeg, 'un', 0);
        
        if nSpf > 2
            spf_ticks = [1,round(nSpf/2),nSpf];  % will then change labels using spf_ticklabels
        else
            spf_ticks = [1,nSpf];
        end
        spfreqs = 1:nSpf;
        spf_forPlot = iff(~isempty(spf_ticklabels), spf_ticks, spfreqs);        
    else
        R = OSP;
        ori_label = 'Orientation';
        [nOri, nSpf, nPh] = size(R);
        oris   = (0:nOri-1)*180/nOri;
        spfreqs = linspace(1, 10, nSpf);                
        phases = (0:nPh-1)*360/nPh;
        
        ori_ticks = [0:45:180];  % spfreqs    = [5 8 11 17 25 38 56 85 127 192];                
        
        ori_forPlot = oris;
        spf_forPlot = spfreqs;

        if nSpf > 2
            spf_ticks = [1,round(nSpf/2),nSpf];  % will then change labels using spf_ticklabels
        else
            spf_ticks = [1,nSpf];
        end        
        
        [ori_ticklabels, spf_ticklabels] = deal([]);
    end
        
    phTicks = [0:90:359];
    axvals = {oris, spfreqs, phases};
    axvals_forPlot = {ori_forPlot, spf_forPlot, phases};
    axticks = {ori_ticks, spf_ticks, phTicks};
    axticklabels = {ori_ticklabels, spf_ticklabels, []};
    
    sizeR = size(R);
    if length(sizeR) == 2, sizeR(3) = 1; end
    
    axlabels = {ori_label, 'Spatial Frequency', 'Phase'};    
    axnames = {'ori','sp','phase'};
    axsym1 = {'\theta', 'f', '\phi'};
    axsym2 = {'\circ', '', '\circ'};        
    
    noTicks = any(strcmpi('noTicks', varargin));
    noLabels = any(strcmpi('noLabels', varargin));
    arg_h = cellfun(@(arg) isscalar(arg) && isnumeric(arg) && ishandle(arg), varargin);
    existingPlot = any(arg_h);
         
%     if existingPlot
%         h_ax = arg_h(find(arg_h,1));
%         if strcmp(get(h_ax, 'type'), 'axes')
%     end
    
    
    if noLabels;
        axlabels = {'', '', ''};
    end
    if noTicks
        axticks = {[], [], []};
    end
    doColorbar = ~( noTicks || noLabels );
    
    % 1. Determine the *order* in which to plot the axes.
    origOrder = 'OSP';
    if ~isempty( setdiff(axOrder, 'OSP'))
        error('Only "O", "S" and "P" allowed in axisOrder parameter');
    end
    ospOrder = arrayfun(@(s) find( s == origOrder, 1), axOrder);
%     ospOrder = zeros(1,3);
%     for i = 1:length(axOrder) 
%         ospOrder(i) = find( axOrder(i) == origOrder, 1);    % eg. "OSP" -> 123
%     end                                                     %     "SPO" -> 231
    if length(ospOrder) < 3                                   %     "OS"  -> 120  [--> 123]
        missingDims = setdiff([1:3], ospOrder);
        ospOrder([end+1:end+length(missingDims)]) = missingDims; % put missing dims at the end;
    end
    
    % 1b. and then switch order of axes if necessary.        
    R = permute(R, ospOrder);           
    axvals = axvals(ospOrder);        
    axvals_forPlot = axvals_forPlot(ospOrder);
    axlabels = axlabels(ospOrder);
    axnames = axnames(ospOrder);
    axsym1 = axsym1(ospOrder);
    axsym2 = axsym2(ospOrder);
    axticks = axticks(ospOrder);
    axticklabels = axticklabels(ospOrder);
%     axticks_init = axticks_init(ospOrder);
    sizeR = sizeR(ospOrder);

    
    
    % 2. Checking *which* dimensions to use.
    % Check how many non-singleton dimensions there are. (There might be less than 3.)     
    nonSingletons = size(R) > 1;
    if length(nonSingletons) == 2    	
        nonSingletons = [nonSingletons, 0];
    end
    nDataDims = nnz(nonSingletons);
    collapseDim = [];    
    
    if nDataDims == 3   % use the user-defined method (or default) of showing the 3rd dimension
    
        % determine how the 3rd dimension is to be plotted.        
        [plot_mode, sub_mode] = strtok(plotOptions, ':');
        if any(strcmp(plot_mode, {'subplots', 'sameaxes'}))
            hv_mode = sub_mode(2:end);
        elseif any(strcmp(plot_mode, {'mean', 'max', 'pref'}))
            collapseDimName = sub_mode(2:end);
            collapseDim = find( strncmpi(collapseDimName,axnames,2) );
        else
            error('Invalid option for "plotOptions" parameter');
        end
        
    else        % don't need to collapse 1 dimension - 1 is already singleton
        collapseDim = find(~nonSingletons, 1);   
        collapseDimName = axnames{collapseDim};
        plot_mode = ['mean:' collapseDimName];        
    end
           


    % 3. Sometimes, we only have to draw one set of axes:
    if ~isempty(collapseDim) % only need to do one plot        
        
        collapseDim = find(strcmp(collapseDimName, axnames), 1);
        idx = 1:3;
        idx(collapseDim) = [];        
        
        
        [h_ax, h_im] = imagesc3(axvals_forPlot, R, collapseDim, plot_mode);
        xlabel(axlabels{idx(2)}, 'fontsize', fontsize );
        ylabel(axlabels{idx(1)}, 'fontsize', fontsize );
        
        set(gca, 'xtick', axticks{idx(2)});
        set(gca, 'ytick', axticks{idx(1)});        
        setTickLabels('x', axticklabels{idx(2)}, axticks{idx(2)});
        setTickLabels('y', axticklabels{idx(1)}, axticks{idx(1)});               
        
        if doColorbar
            colorbar;
        end
        axis tight xy;
        return;

    end
        
        
    % otherwise, have to do multiple axes: subplots
    
    switch plot_mode
        case 'subplots'
%             if strcmp(hv_mode , 'auto')
%                 if
            
            switch hv_mode 
                case {'horizontal', 'vertical'}
                    h_im = zeros(1,sizeR(3)+1);
                    h_ax = zeros(1,sizeR(3));
                    % show data for each phase
                    for i3 = 1:sizeR(3)
                        if strcmp(hv_mode, 'horizontal')
                            h_ax(i3) = subplot(1, sizeR(3)*2+1, i3*2-1:i3*2 ); 
                        else                 
                            m = 7;
                            h_ax(i3) = subplot(sizeR(3), m, (i3-1)*m+1 : i3*m-1 ); 
                        end
                        
                        h_im(i3) = imagesc(axvals{2}, axvals{1}, R(:,:,i3));
                        
                        % xyticks, xylabels
                        if strcmp(hv_mode, 'horizontal')
                            title( [axsym1{3} ' = ' num2str(axvals{3}(i3))  axsym2{3} ], 'fontsize', fontsize );
                            
                            if (i3 == 1) 
                                ylabel( axlabels{1}, 'fontsize', fontsize );
                                set(gca, 'xtick', [], 'ytick', axticks{1});
                                setTickLabels('y', axticklabels{1}, axticks{1});
                            elseif (i3 == round(sizeR(3)/2))
                                xlabel( axlabels{2}, 'fontsize', fontsize );
                                set(gca, 'xtick', axticks{2}, 'ytick', []);
                                setTickLabels('x', axticklabels{2}, axticks{2});                                
                            else
                                set(gca, 'xtick', [], 'ytick', []);
                            end
                                                                                
                        elseif strcmp(hv_mode, 'vertical')                            
                            if (i3 == 1)                                
                                set(gca, 'xtick', [], 'ytick', axticks{1});
                                setTickLabels('y', axticklabels{1}, axticks{1});
                            elseif (i3 == round(sizeR(3)/2) )
                                ylabel( axlabels{1}, 'fontsize', fontsize);
                                set(gca, 'xtick', [], 'ytick', []);                                
                            elseif (i3 == sizeR(3) )
                                xlabel( axlabels{2}, 'fontsize', fontsize );                                
                                set(gca, 'xtick', axticks{2}, 'ytick', []);
                                setTickLabels('x', axticklabels{2}, axticks{2});                                
                            else
                                set(gca, 'xtick', [], 'ytick', []);
                            end
                            
                            
                        end
                    end

                    if doColorbar
                        %show colorbar in last space
                        if strcmp(hv_mode, 'horizontal')
                            subplot(1, sizeR(3)*2+1, sizeR(3)*2+1); 
                        else
                            subplot(sizeR(3), m, m:m:m*sizeR(3) ); 
                        end
                        h_im(sizeR(3)+1) = imagesc([min(R(:)), max(R(:))]);
                        set([h_im(sizeR(3)+1), get(h_im(sizeR(3)+1), 'Parent')], 'Visible', 'off');
                        ch = colorbar;
                        [l, b, w, h] = elements(get(ch, 'Position'));
                        set(ch, 'Position', [l*.95, b, w*3, h]);
                        matchAxes('C', h_im);
                    end
                    
                    
                case 'diagonal' % for slides purposes only
        
                    if size(R, 2) > 10
                        q = 4;
%                         [R, inds] = compressAlongDim(R, 2, q); 
%                         axvals{2} = axvals{2}(inds);
%                         sizeR(2) = sizeR(2)/q;
                    end
                    maxDiags = 5;                    
%                     inds3 = floor(linspace(1,sizeR(3), maxDiags));
%                     N3 = sizeR(3);
                    N3 = maxDiags;
                    [tmp, imax3] = max( sum(sum(R, 1),2));
%                     [i1, i2, imax3] = elements(indmax);
%                     inds3 =  mod(imax3 + [round(-maxDiags/2):round(maxDiags/2)], sizeR(3));
%                     inds3 = inds3(1:maxDiags);
                    skp = max(1, floor(sizeR(3)/maxDiags));
                    inds3 = mod( imax3-1 + skp*[0:maxDiags-1], sizeR(3) )+1;
                    
%                     inds3 = 1:5;%sizeR(3); 
                    h_ax = zeros(1,N3); 
                    h_im = zeros(1,N3);
                    for i3 = 1:length(inds3)
                        h_ax(i3) = subplot(1, N3, i3); 
                        h_im(i3) = imagesc(axvals{2}, axvals{1}, R(:,:,inds3(i3)));
%                         set(h_ax(i3), 'fontsize', fontsize);
%                         axis tight xy;
%                         daspect( [xy_ratio, (length(axvals{2})*diff(ylim))/(length(oris)*diff(xlim)), 1]);
                        y_ratio = (sizeR(2))*diff(ylim)/(sizeR(1)*diff(xlim));
                        if fd == 2
                            y_ratio = y_ratio/4;
                        end
                        daspect( [1, y_ratio, 1]);
%                        
                        title( [axsym1{3} ' = ' num2str(axvals{3}(inds3(i3)))  axsym2{3} ], 'fontsize', fontsize );
                        if i3 == 1,
                            xlabel( axlabels{2}, 'fontsize', fontsize );
                            ylabel( axlabels{1}, 'fontsize', fontsize );
                            set(gca, 'ytick', axticks{1}, 'xtick', axticks{2});
                            setTickLabels('y', axticklabels{1}, axticks{1});
                            setTickLabels('x', axticklabels{2}, axticks{2});
                        else
                            set(gca, 'ytick', [], 'xtick', []);
                        end
                        
                    end 
                    if doColorbar
                        hcol = colorbar;
                        h_ax = [h_ax hcol];
                    end
                    set(gcf, 'children', h_ax); % reorder
                    
                    set(h_ax(1), 'Units', 'Pixels'); s = get(h_ax(1), 'Position'); 
                    dim_ratio = s(3)/s(4);
                    wid = 100;
                    dims = [wid, wid/dim_ratio];
                    
                    dx = 80;
                    dy = dx/3;
                    lb = [-10 -180];
                    for i3 = 1:N3
                        set(h_ax(i3), 'Units', 'Pixels', 'Position', [lb +  s(1:2) + [dx*(i3-1), dy*(i3-1)], dims]);                            
                    end                
                    matchAxes('C', h_im);                    
                    
                                        
                    
            end
    
        case 'sameaxes' % don't really use this anymore
            if strcmp(hv_mode, 'vertical')
                flatOSP = zeros(sizeR(1)*sizeR(3), sizeR(2));
                for i3 = 1:sizeR(3)
                    flatOSP( [1:sizeR(1)] + (i3-1)*sizeR(1), [1:sizeR(2)]) = R(:,:,i3);
                end

                imagesc(flatOSP);
                set(gca, 'ytick', ([1:i3]-0.5)*sizeR(1), 'yticklabel', num2str([1:sizeR(3)]')  );
                set(gca, 'xtick', []);

                drawHorizontalLine( (1:i3-1)*sizeR(1), 'LineWidth', 2);
            elseif strcmp(hv_mode, 'horizontal')

            end
       
%         case 'mean:Phase'
%             h_im = imagesc(axvals{2}, axvals{1}, mean(R, 3));
%             axis tight xy;
%             daspect( [xy_ratio, (length(axvals{2})*diff(ylim))/(length(axvals{1})*diff(xlim)), 1]);
%             set(gca, 'ytick', oris_disp, 'xtick', []);
%             xlabel('Spatial Frequency');                                
%             ylabel('Orientation');
% %                         
% %             , 'mean:Orientation', 'mean:SpatFreq'}
% %             switch 
            
            
            
        otherwise
            error(['Unknown plot_mode parameter "' plot_mode '"']);
    end
    
end


% function [Y, inds] = compressAlongDim(X, dim, n)
%     assert(mod(n,2) == 0);
%     sizeX = size(X);
%     sizeX(dim) = sizeX(dim)/n;
%     Y = zeros(sizeX);
%     sizeR = size(X, dim);
%     idxY = arrayfun(@(n) 1:n, size(X), 'Un', false);
%     idxX = idxY;
%         
%     for i = 1:sizeR/2    
%         idxY{dim} = i;        
%         idxX{dim} = i*2-1:i*2;        
%         Y(idxY{:}) = sum(X(idxX{:}), dim);
%     end
%     inds = 1:n:sizeR;
% end



function setTickLabels(ax, axlab, axlabInds, axticks)
    if ~isempty(axlab)
        set(gca, [ax 'tick'], axlabInds);
        set(gca, [ax 'ticklabel'], axlab(axlabInds));        
    else
%         set(gca, [ax 'tick'], axticks)
    end
end


% function imageOSP(R)
%     function xd = deltafunAtMax(x)
%         xd = zeros(size(x));
%         if any(x > 0)
%             idxmax = indmax(x);
%             xd(idxmax) = 1;
%         end
%     end
% 
%     [nOri, nSpf, nPh] = size(R);
%     reshpOSP = reshape(R, [nOri*nSpf, nPh])';
% %     for i = 1:nOri*nSpf
% %         reshpOSP(:,i) = deltafunAtMax(reshpOSP(:,i));
% %     end
%     imagesc(reshpOSP);
% z    colormap('gray');
%     
% end
% 

