function gb = findMostInformativeGabor(cellInfo, gaborGuess0)

    dbug = false;

    pLabels = {'A', 'mu_x', 'sig_x', 'mu_y', 'sig_y', 'k', '\phi', '\theta', 'c'};

    Gid = cellInfo.Gid;
%     cellId = cellInfo.cellId;
    Pspike_s = cellInfo.Pspike_s;
    getFrame = getFrameRetrieverFunction(Gid);
    
    getFrame('load', Gid, 'scaling', 'aroundZero');
    nStimulusFrames = getFrame('nTotalFrames');
	[nrows ncols] = elements(getFrame('size'));
    xs = 1:ncols;
    ys = 1:nrows;
    npix = nrows*ncols;
    
    %     [nrows, ncols, nframes, framesPerSecond, averIntensity] =
    %     getFieldsFromDatabaseTable(hnd, {'LNG_N_ROWS', 'LNG_N_COLUMNS',
    %     'LNG_N_FRAMES', 'DBL_FRAME_RATE_HZ', 'DBL_AVER_INTENSITY'}, ...
    %         'TBL_MOVIES', {'TXT_MOVIE_FILE_NAME', movieFileName} );

%     npixels = nrows * ncols;

    nframes = length(Pspike_s);
    assert(nframes == nStimulusFrames);
    [xs_grid, ys_grid] = meshgrid(xs, ys);

    function varargout = getDistributionProjections(v)
%        {PvX_sp_and_PvX}      = getDistributionProjections(v); 
%        [PvX_sp, PvX]         = getDistributionProjections(v);
%        [PvX_sp, PvX, Pxbins] = getDistributionProjections(v);

        if isvector(v)
            V = gabor(v, {xs_grid, ys_grid});
        elseif isamatrix(v)
            V = v;
        end

        xvalues = zeros(1, nframes);
        xvalues_sp = zeros(1, nnz(Pspike_s));
        spikeInd = 1;

        for frm_i = 1:nframes
            s = getFrame(frm_i);
            x = sum(s(:) .* V(:));
            xvalues(frm_i) = x;

            if (Pspike_s(frm_i) > 0)
                xvalues_sp(spikeInd) = x;% * Pspike_s(frm_i);    % ** double check
                spikeInd = spikeInd + 1;
            end
        end;

        % to get bins appropriate for both distributions
        [tmp, xbins] = hist([xvalues; xvalues_sp], 40);
        PvX = histc(xvalues, xbins);
        PvX_sp = histc(xvalues_sp, xbins);
        PvX = PvX / sum(PvX);
        PvX_sp = PvX_sp / sum(PvX_sp);

        if (nargout == 1)
            varargout =    { {PvX_sp, PvX} };
        elseif (nargout == 2)
            varargout =    {PvX_sp, PvX};
        elseif (nargout == 3)
            varargout =    {PvX_sp, PvX, xbins};
        end

    end


    function gradI = getInformationGradient(v)
        GABOR_SEARCH = 1;
        FULL_SEARCH = 2;
        if isvector(v)
            searchType = GABOR_SEARCH;
            V = gabor(v, {xs_grid, ys_grid});
            dGabor_dv = gabor(v, [xs_grid(:), ys_grid(:)], 'gradient');
            
        elseif isamatrix(v)
            searchType = FULL_SEARCH;
            V = v;
        end
                
        [PvX_sp, PvX, Pxbins] = getDistributionProjections(V);
        Pspike = mean(Pspike_s);
        Ps = 1/nframes;
        Ps_spike =  Pspike_s * (Ps / Pspike);   % ****
                
        s_x    = zeros(npix, length(PvX));
        s_x_sp = zeros(npix, length(PvX));
        
        for frm_i = 1:nframes
            s = getFrame(frm_i); 
            s = s(:);
            
            x = sum(s .* V(:));
            xbin = find( Pxbins <= x, 1, 'last' );  % xbin = binarySearch(Pxbins, x, -1, -1);
            if isempty(xbin), xbin = 1; end

            s_x(:,xbin) = s_x(:,xbin) + s * Ps;
            if (Pspike_s(frm_i) > 0)
                s_x_sp(:,xbin) = s_x_sp(:,xbin) + s * Ps_spike(frm_i);
            end
        end
        s_x    = s_x    / bsxfun(@rdivide, s_x, PvX);
        s_x_sp = s_x_sp / bsxfun(@rdivide, s_x, PvX_sp);
        
        s_x_sp__s_x = s_x_sp - s_x;
        
        if searchType == GABOR_SEARCH
            % (npix x nbins) x (npix x 8) -> (npix x nbins x 8) -> (1 x nbins x 8) -> (nbins x 8) -> (8 x nbins) 
            s_x_sp__s_x = squeeze( sum( bsxfun(@times, s_x_sp__s_x,  dGabor_dv), 1))'; 
        end
                
        ddx_PvXsp_PvX = derivative( Pxbins , PvX_sp ./ (PvX + eps) );
        integrand = bsxfun( @times,  s_x_sp__s_x,   PvX .* ddx_PvXsp_PvX );

        
        gradI = sum(integrand,2);
        gradI = gradI / norm(gradI);
        
    end


    gradI = @(v) getInformationGradient(v);
    I = @(v) mutualInfo(getDistributionProjections(v));

    R = Pspike_s / mean(Pspike_s);
    [PvX_sp, PvX] = getDistributionProjections(gaborGuess0);

            if dbug
                figure(5);
                bar(xbins, PvX, 'EdgeColor', 'b', 'FaceColor', 'none'); hold on
                bar(xbins, PvX_sp, 'EdgeColor', 'g', 'FaceColor', 'none'); hold off
                legend('P_v(x)', 'P_v(x|spike)');


                figure(6)
                bar(xbins, PvX_sp ./ PvX);
            end

    Nrep = 1;

    I_spike = mean( xlogy(R, R) );
    biasValue = 1/(mean(Pspike_s)*Nrep*2*log(2));

    %     xmax = gradientAscent(gaborv, I, gradI);

    %confirm that true gabor is at the peak of the information;
    testMax = false;
    if testMax
        Itrue = I(gaborGuess0); %%% use final value of gabor estimate
        delta = .5;
        for pi = 1:length(actualGaborParams);
            paramPlus = actualGaborParams;   paramPlus(pi) = paramPlus(pi) + paramPlus(pi)*delta;
            paramMinus = actualGaborParams;  paramMinus(pi) = paramMinus(pi) - paramMinus(pi)*delta;
            Iplus = I(paramPlus);
            Iminus = I(paramMinus);
            if (Itrue >= Iplus)
                disp([' == I_true > I_plus for ' pLabels{pi} ' by ' num2str(Itrue-Iplus)]);
            elseif (Itrue < Iplus)
                disp([' ** I_true < I_plus for ' pLabels{pi} ' by ' num2str(Iplus-Itrue)]);
            end

            if (Itrue >= Iminus)
                disp([' == I_true > I_minus for ' pLabels{pi} ' by ' num2str(Itrue-Iminus)]);
            end
            if (Itrue < Iminus)
                disp([' ** I_true < I_minus for ' pLabels{pi} ' by ' num2str(Iminus-Itrue)]);
            end
            fprintf('\n');
        end
    end
        
    gaborWithMaxInfo = gradientAscent(I, gradI, gaborGuess0, {xs, ys}, I_spike, biasValue);
   
    getFrame('close');
    gb = gaborWithMaxInfo;
end





%                 Sx = Sx + (s * x) * (1/nframes);
%                 if (Pspike_s(frm_i) > 0)
%                     Sx_sp = Sx_sp + (s * x) * Ps_spike(frm_i);
%                 end
