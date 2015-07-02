function [fit_func, fit_S] = getCorticalState(Gid, cellId, fitType, nTerms, redoFlag)

    persistent allCorticalStates saveCount

    if nargin < 3 || isempty(fitType)
        fitType = 'sin';
%         fitType = 'gauss';
    end
    if nargin < 4 || isempty(nTerms)
        nTerms = 8;
    end
    
    opt.fitType = fitType;
    opt.nTerms = nTerms;

%     opt.freq_upper_bound_Hz = 0.5/60; % half a cycle per minute (1 cycle every 2 minutes);
    opt.freq_upper_bound_Hz = 1/60; % one cycle per minute ;
    
    opt.windowSize_sec = 5;
    opt.nStepsPerWindow = 10; 
    
    opt.excludeBckgSpikes = 1;
    opt.excludeFromWindow_ms = 100;
    opt.fitToRectifiedCurves= true;
%     opt.rectifyAt = 0.05;
    opt.rectifyAt = 0.1;
    
%     opt.rectType = 'hard';
    opt.rectType = 'softplus';
    opt.softPlusScale = 10;
    
%     isDrifting = flashedOrDrifting(Gid) == 2;
%     if isDrifting
    opt.fracStimInFullCycle = 0.5;
%     end
    
    corticalStates_file = [CatV1Path 'MatLabDB_avi' filesep 'allCorticalStates.mat'];
        
    redo_all = 0;
    redo_current = 0 || exist('redoFlag', 'var') && ~isempty(redoFlag);
    saveCountSpacing = 25;
    
    if strcmp(Gid, 'save')
        %%
        save(corticalStates_file, 'allCorticalStates', '-v6');                
        saveCount = 0;
        return;
    end        
    
    if isempty(allCorticalStates)        
        if exist(corticalStates_file, 'file') && ~redo_all
            S_file = load(corticalStates_file);
            allCorticalStates = S_file.allCorticalStates;
        else
            allCorticalStates = struct;
        end        
        saveCount = 0;
    end
    
%     if nargin < 2
%         cellId = 100;
%     end
%     siteState = cellId == 100;
    
    funcFitName = sprintf('%s%d', opt.fitType, nTerms);
    group_fld_name = sprintf('CorticalState_Gid_%d_cellId_%d__%s', Gid, cellId, funcFitName);
    group_fld_name = strrep(group_fld_name, '-1', 'n1');    
    
    if (~isfield(allCorticalStates, group_fld_name) || redo_current || ~isequal(allCorticalStates.(group_fld_name).opt, opt) )
        fit_S = calcCorticalState(Gid, cellId, opt);
        
        allCorticalStates.(group_fld_name) = fit_S;            

        saveCount = saveCount + 1;
        if saveCount > saveCountSpacing
            save(corticalStates_file, 'allCorticalStates', '-v6');        
            saveCount = 0;
        end                

    end
    
    fit_S = allCorticalStates.(group_fld_name);
    fit_func = fit_S.func_fit;
    
    
end


function corticalState_S = calcCorticalState(Gid, cellId, opt)

    if (nargin < 2) || isempty(cellId)
        cellId = 100;
    end

    freq_upper_bound_Hz = opt.freq_upper_bound_Hz;

%     show = 0; %exist('show_flag', 'var') && ~isempty(show_flag);
    if isfield(opt, 'fracStimInFullCycle')
        fracStimuli = opt.fracStimInFullCycle;
        recur_t = getStimulusRecurrenceTime(Gid, fracStimuli);
        
        freq_upper_bound_stim_Hz = 1 / recur_t; 
        
        freq_upper_bound_Hz = min(freq_upper_bound_Hz, freq_upper_bound_stim_Hz);
        
    end
    
    freq_upper_bound_rad = (freq_upper_bound_Hz*2*pi); % for sine fit
    freq_lower_bound_sec = 1/freq_upper_bound_Hz;      % for gauss fit

    
    excludeFromWindow_sec = opt.excludeFromWindow_ms/1000;
        
    allSpk = getSpikes(Gid, cellId, 'sec');
    
    [~, movieStarts_sec, movieEnds_sec] = getSyncs(Gid, 'sec');
    allMovieLengths_sec = movieEnds_sec - movieStarts_sec;
    minMovieLength_sec = min(allMovieLengths_sec);
    
%     tf_bckgSpike = allSpk < startTicks(1) | allSpk > endTicks(end);
    
    %%
    if opt.excludeBckgSpikes
        tf_bckgSpike = true(size(allSpk));
        for i = 1:length(movieStarts_sec)
            tf_bckgSpike = tf_bckgSpike & ~ibetween(allSpk, movieStarts_sec(i)+excludeFromWindow_sec, movieEnds_sec(i));
        end
        
        allSpk = allSpk(~tf_bckgSpike);
    end
    %%
        
    windowSize_sec = opt.windowSize_sec;
    windowSize_sec = min(windowSize_sec, minMovieLength_sec/2);
    
    
    windowStep_sec = windowSize_sec/opt.nStepsPerWindow;
    
    sd = siteDataFor(Gid);
    t_end_sec = sd.dataFileInfo.duration_sec;
    
    slidingWindows_start_all = [0 : windowStep_sec : t_end_sec-windowSize_sec];
    slidingWindows_end_all   = slidingWindows_start_all + windowSize_sec;
    if opt.excludeBckgSpikes
        tf_useWindow = false(size(slidingWindows_start_all));
        for i = 1:length(movieStarts_sec)
            tf_useWindow = tf_useWindow | ibetween(slidingWindows_start_all, movieStarts_sec(i)+excludeFromWindow_sec, movieEnds_sec(i)) ...
                                        & ibetween(slidingWindows_end_all,   movieStarts_sec(i)+excludeFromWindow_sec, movieEnds_sec(i));
        end        
    else
        tf_useWindow = true(size(slidingWindows_start_all));
    end
    slidingWindows_start = slidingWindows_start_all(tf_useWindow);
    slidingWindows_end   = slidingWindows_end_all(  tf_useWindow);        
    
    %%
    show = 0;
    if show
       %%
       figure(93); clf; hold on;
       plot(slidingWindows_start_all(tf_useWindow), ones(1, nnz(tf_useWindow))*.9, 'o')
       plot(slidingWindows_end_all(tf_useWindow), ones(1, nnz(tf_useWindow))*1, 'ro')
       plot(slidingWindows_start_all(~tf_useWindow), ones(1, nnz(~tf_useWindow))*0, 'o')
       plot(slidingWindows_end_all(~tf_useWindow), ones(1, nnz(~tf_useWindow))*.1, 'ro')
       drawVerticalLine(movieStarts_sec)
       drawVerticalLine(movieEnds_sec)
        
    end
    
    binCents = (slidingWindows_start + slidingWindows_end)/2;
    nSpks_binned = elementsInRange(allSpk, [slidingWindows_start(:), slidingWindows_end(:)], 'count');
    spkRate_binned = nSpks_binned / windowSize_sec;
    
%     w = (windowSize_sec)/diff(binCents(1:2));
%     w = (windowSize_sec)/diff(binCents(1:2));
    w = freq_lower_bound_sec/5;
    spkRate_bin_sm = gaussSmooth_nu(binCents, spkRate_binned, w);
    
%     spkRate_bin_toFit = spkRate_binned;
    spkRate_bin_toFit = spkRate_bin_sm;
    
    scale_factor = 1/mean(spkRate_bin_toFit(:));
    spkRate_bin_toFit_scl = spkRate_bin_toFit(:) * scale_factor;
        
    nTerms = opt.nTerms;
    doBothFits = 0;
%     doGaussFit = 

    doSineFit = strcmp(opt.fitType, 'sin') || doBothFits;
    doGaussFit = strcmp(opt.fitType, 'gauss') || doBothFits;
    if doSineFit
        % put upper limit on the frequency of each cosine
        sine_fitFunc_name = sprintf('sin%d', nTerms);
        sine_upper_bnds = inf(1, nTerms*3); sine_upper_bnds([1:nTerms]*3-1) = freq_upper_bound_rad;  
        sine_fit_options = fitoptions(sine_fitFunc_name, 'upper', sine_upper_bnds);
        
        sine_fit = fit(binCents(:), spkRate_bin_toFit_scl, sine_fitFunc_name, sine_fit_options);
        sine_yfit = feval(sine_fit, binCents);
        
        orig_sine_fit = [];
        min_unrectfitValue = min(sine_yfit);
        fracRectified = nnz(sine_yfit < opt.rectifyAt) / numel(sine_yfit);
        
        
        doRectifiedSinusoids = any(sine_yfit < opt.rectifyAt) && opt.fitToRectifiedCurves;
        if doRectifiedSinusoids
            orig_sine_fit = sine_fit;
            
            sine_terms = arrayfun(@(i) sprintf('a%d*sin(b%d*x+c%d)', i, i, i), 1:nTerms, 'un', 0);
            switch opt.rectType 
                case 'hard'
                    rectSine_str = sprintf('rectified( %s - %.3f)+%.3f', cellstr2csslist(sine_terms, ' + '), opt.rectifyAt, opt.rectifyAt);
                case 'softplus';
                    rectSine_str = sprintf('softplus( %s, %f, %f)', cellstr2csslist(sine_terms, ' + '), opt.rectifyAt, opt.softPlusScale );
            end

            
            fit_coeffs_C = arrayfun(@(i) {['a' num2str(i)], ['b' num2str(i)], ['c' num2str(i)]}, 1:nTerms, 'un', 0);
            fit_coeffs = [fit_coeffs_C{:}];
           
            rect_sin_fittype = fittype(rectSine_str, ...
                'dependent', {'y'}, 'independent', {'x'}, 'coefficients', fit_coeffs );

%             startPoint = coeffvalues(sine_fit);
            sineStartPoint = sinnstart(binCents(:), spkRate_bin_toFit_scl, nTerms);
            sine_fit = fit(binCents(:), spkRate_bin_toFit_scl, rect_sin_fittype, 'startpoint', sineStartPoint, 'upper', sine_upper_bnds, 'lower', sine_fit_options.Lower);
        end
        3;
    end
    
%     sine_yfit = feval(sine_fit, binCents);        
   
    %%
%     doGaussFit = doGaussFit || (doSineFit && doRectifiedSinusoids);
    if doGaussFit
        % put lower limit on the width of each gaussian
        %%
        gauss_fitFunc_name = sprintf('gauss%d', nTerms);
        gauss_lower_bnds = -inf(1, nTerms*3); gauss_lower_bnds([1:nTerms]*3) = freq_lower_bound_sec/4;  
        gauss_fit_options = fitoptions(gauss_fitFunc_name, 'lower', gauss_lower_bnds);
        gauss_fit = fit(binCents(:), spkRate_bin_toFit_scl, gauss_fitFunc_name, gauss_fit_options);
        
        gauss_yfit = feval(gauss_fit, binCents);
        
        %%
        orig_gauss_fit = [];
        doRectifiedGaussians = any(gauss_yfit < 0) && opt.fitToRectifiedCurves;
        
        if doRectifiedGaussians
            orig_gauss_fit = gauss_fit;
            
            gauss_terms = arrayfun(@(i) sprintf('a%d*exp(-((x-b%d)/c%d)^2)', i, i, i), 1:nTerms, 'un', 0);
            
            switch opt.rectType 
                case 'hard'
                    rectGauss_str = sprintf('rectified( %s - %.3f)+%.3f', cellstr2csslist(gauss_terms, ' + '), opt.rectifyAt, opt.rectifyAt);
                case 'softplus';
                    rectGauss_str = sprintf('softplus( %s, %f, %f)', cellstr2csslist(gauss_terms, ' + '), opt.rectifyAt, opt.softPlusScale );
            end
            
            
            
            fit_coeffs_C = arrayfun(@(i) {['a' num2str(i)], ['b' num2str(i)], ['c' num2str(i)]}, 1:nTerms, 'un', 0);
            fit_coeffs = [fit_coeffs_C{:}];            
            
            rect_gauss_fittype = fittype(rectGauss_str, 'dependent', {'y'}, 'independent', {'x'}, 'coefficients', fit_coeffs );

%             startPoint = coeffvalues(sine_fit);
            gaussStartPoint = gaussnstart(binCents(:), spkRate_bin_toFit_scl, nTerms);
            gauss_fit = fit(binCents(:), spkRate_bin_toFit_scl, rect_gauss_fittype, 'startpoint', gaussStartPoint, 'lower', gauss_fit_options.Lower);
        end
    end
    

    switch opt.fitType
        case 'sin',  func_fit = sine_fit;  
            func_fit_unrect = orig_sine_fit;
            
        case 'gauss',func_fit = gauss_fit;  
            func_fit_unrect = orig_gauss_fit;
    end

    show = 1;

%     if (doSineFit && doRectifiedSinusoids) || (doGaussFit && doRectifiedGaussians)
%         beep;
%         fprintf('Gid = %d. cellId = %d. NEGATIVE STATE\n', Gid, cellId);
%         show = 1;
%     end
%     S = cell2struct(num2cell(coeffvalues(state_func))', coeffnames(state_func), 1);

    corticalState_S.opt = opt;
    corticalState_S.binCents = binCents;
    corticalState_S.nSpks_binned = nSpks_binned;
    corticalState_S.spkRate_binned = spkRate_binned;
    corticalState_S.spkRate_binned_scl = spkRate_bin_toFit_scl;
    corticalState_S.func_fit = func_fit;
    corticalState_S.func_fit_unrect = func_fit_unrect;
    
    corticalState_S.min_unrectfitValue = min_unrectfitValue;
    corticalState_S.fracRectified = fracRectified;
    
    
    if show
        %%
%         figure(Gid*1000+cellId + 1);
        figure(1001);
        clf; hold on; box on;
%         col = iff(cellId == 100, 'b', color_s(cellId));
        col = iff(cellId == 100, 'b', 'b');
        
        plot(binCents, spkRate_binned, 'b');
        plot(binCents, spkRate_bin_sm, 'k.', 'linewidth', 2);

%         f_rect = @(a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4,a5,b5,c5,a6,b6,c6,a7,b7,c7,a8,b8,c8,x) ...
%             rectified(a1.*sin(b1.*x+c1) + a2.*sin(b2.*x+c2) + a3.*sin(b3.*x+c3) + a4.*sin(b4.*x+c4) + a5.*sin(b5.*x+c5) + a6.*sin(b6.*x+c6) + a7.*sin(b7.*x+c7) + a8.*sin(b8.*x+c8) -.05)+.05;
%         f_norect = @(a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4,a5,b5,c5,a6,b6,c6,a7,b7,c7,a8,b8,c8,x) ...
%             (a1.*sin(b1.*x+c1) + a2.*sin(b2.*x+c2) + a3.*sin(b3.*x+c3) + a4.*sin(b4.*x+c4) + a5.*sin(b5.*x+c5) + a6.*sin(b6.*x+c6) + a7.*sin(b7.*x+c7) + a8.*sin(b8.*x+c8));
%         
%         [a1A,b1A,c1A,a2A,b2A,c2A,a3A,b3A,c3A,a4A,b4A,c4A,a5A,b5A,c5A,a6A,b6A,c6A,a7A,b7A,c7A,a8A,b8A,c8A] = dealV(coeffvalues(orig_sine_fit));
%         [a1B,b1B,c1B,a2B,b2B,c2B,a3B,b3B,c3B,a4B,b4B,c4B,a5B,b5B,c5B,a6B,b6B,c6B,a7B,b7B,c7B,a8B,b8B,c8B] = dealV(coeffvalues(sine_fit));
%         
%         f_noRect_fitNoRect = @(x) f_norect(a1A,b1A,c1A,a2A,b2A,c2A,a3A,b3A,c3A,a4A,b4A,c4A,a5A,b5A,c5A,a6A,b6A,c6A,a7A,b7A,c7A,a8A,b8A,c8A, x);
%         f_rect_fitNoRect   = @(x) f_rect(  a1A,b1A,c1A,a2A,b2A,c2A,a3A,b3A,c3A,a4A,b4A,c4A,a5A,b5A,c5A,a6A,b6A,c6A,a7A,b7A,c7A,a8A,b8A,c8A, x);
%         f_noRect_fitRect   = @(x) f_norect(a1B,b1B,c1B,a2B,b2B,c2B,a3B,b3B,c3B,a4B,b4B,c4B,a5B,b5B,c5B,a6B,b6B,c6B,a7B,b7B,c7B,a8B,b8B,c8B, x);
%         f_rect_fitRect     = @(x) f_rect(  a1B,b1B,c1B,a2B,b2B,c2B,a3B,b3B,c3B,a4B,b4B,c4B,a5B,b5B,c5B,a6B,b6B,c6B,a7B,b7B,c7B,a8B,b8B,c8B, x);
%         
        
%         rect_sine_yfit = feval(rect_sine_fit, binCents);
        hold on;
        legs = {};
%         plot(binCents, gauss_yfit, 'c', 'linewidth', 2)   
        if doSineFit
            sine_yfit = feval(sine_fit, binCents);        
            plot(binCents, sine_yfit / scale_factor, 'r.', 'linewidth', 2)   
            legs = [legs, {'Sine fit'}];
            if doRectifiedSinusoids
                
                nonrect_sine_yfit = feval(orig_sine_fit, binCents);        
                plot(binCents, nonrect_sine_yfit / scale_factor, 'm.', 'linewidth', 2)   
                legs = [legs, {'non-Rect-sine fit'}];
            end
        end
        
        if doGaussFit
            gauss_yfit = feval(gauss_fit, binCents);        
            plot(binCents, gauss_yfit / scale_factor, 'g', 'linewidth', 2)   
            legs = [legs, {'Gauss fit'}];
            
            
            if doRectifiedGaussians
                
                nonrect_gauss_yfit = feval(orig_gauss_fit, binCents);        
                plot(binCents, nonrect_gauss_yfit / scale_factor, 'r', 'linewidth', 2)   
                legs = [legs, {'non-Rect-Gauss fit'}];
            end
            
        end
                
         if exist('recur_t', 'var')
            period_lower_bound_sec_init = 1/opt.freq_upper_bound_Hz;
            drawVerticalLine(period_lower_bound_sec_init, 'color', 'k')
        end
        period_lower_bound_sec = 1/freq_upper_bound_Hz;        
        drawVerticalLine(period_lower_bound_sec, 'color', 'r')
        
        xlim([0 t_end_sec]);
        title(sprintf('Spikes for Group %d, cell %d', Gid, cellId));
        xlabel('Time (seconds)'); 
        ylabel('Firing rate (Hz)');
        legend({'Binned firing rate', 'Smoothed firing rate', legs{:}}, 'location', 'best')
        3;
        drawnow;
        pause(.1);
    end    
    
    
    
end

%             sd = siteDataFor('Gid', Gid, 1);
%             cellIds = sd.cellIds;
%             cell_linFuncs = cell(1,length(cellIds));
%             for i = 1:length(cellIds)
%                 cell_linFuncs{i} = getCellExcitabilityFunction(Gid, cellId, sine_fit);
%                 
%             end
%             
%             S = struct('sine_fit', sine_fit, 'cell_linFuncs', cell_linFuncs);




%{
allGids = getAllGids;
progressBar('init-', length(allGids), 60);
for i = 1:length(allGids)
    state_func = getCorticalState(allGids(i));
    progressBar(i);
end
getCorticalState('save');

%}




%---------------------------------------------------------
%  SINNSTART
%---------------------------------------------------------
function start = sinnstart(x,y,n)

% SINNSTART start point for SINN library function.
% START = SINNSTART(X,Y,N) computes a start point START based on X.
% START is a column vector, e.g. [p1, p2, ... pn+1,q1,...,qm];
 
% By running the y data through a Fast Fourier transform and then locating
% peak(s) in results, we can find the starting value fo the frequency
% variable 'b'. Because a phase-shifed sine function is separable and can
% be converted to a sum of sine and cosine functions, starting values for
% amplitude 'a' and phase shift 'c' can also be found.

% Get x and y values on a lattice, so we can use fft to find a frequency
[x,y] = getxygrid(x,y);

% Data size too small, cannot find starting values
if length(x) < 2
    start = rand(3*n,1);
    return;
end

% Loop for sum of sines functions
start = zeros(3*n,1);
oldpeaks = [];
lengthx = length(x);
freqs = zeros(n,1);
res = y;   % residuals from fit so far
for j=1:n
    % Apply fast fourier transform to the current residuals
    fy = fft(res);     % don't subtract mean, no constant term
    fy(oldpeaks) = 0;  % omit frequencies already used

    % Get starting value for frequency using fft peak
    [~,maxloc] = max(fy(1:floor(lengthx/2)));
    oldpeaks(end+1) = maxloc; %#ok<AGROW>
    w = 2*pi*(max(0.5,maxloc-1))/(x(end)-x(1));
    freqs(j) = w;
 
    % Compute Fourier terms using all frequencies we have so far
    X = zeros(lengthx,2*j);
    for k=1:j
        X(:,2*k-1) = sin(freqs(k)*x);
        X(:,2*k)   = cos(freqs(k)*x);
    end
    
    % Fit these terms to get the non-frequency starting values
    ab = X \ y(:);
    
    if j<n
        res = y - X*ab;    % remove these components to get next frequency
    end
end

% All frequencies found, now compute starting values from all frequencies
% and the corresponding coefficients
for k=1:n
    start(3*k-2) = sqrt(ab(2*k-1)^2 + ab(2*k)^2);
    start(3*k-1) = freqs(k);
    start(3*k)   = atan2(ab(2*k),ab(2*k-1));
end
end  % sinnstart


%---------------------------------------------------------
function start  = gaussnstart(x,y,n)
% GAUSSNSTART start point for GAUSSN library function.
% START = GAUSSNSTART(X,Y,N) computes a start point START based on X.
% START is a column vector, e.g. [p1, p2, ... pn+1,q1,...,qm];

% The main idea in this computation is to compute one peak at a
% time. Assuming all data is coming from one gaussian model, find out
% the estimated corresponding coefficients for this gaussian, use them as
% the coefficient for the first peak. Then, subtract this peak
% (evaluated at all x) from y data, repeat the procedure for the second
% peak. When we cannot continue for some reason (such as not enough
% significant data, we break, and assign random numbers for the rest of
% the starting point.

x = x(:); y = y(:);
if any(diff(x)<0) % sort x
    [x,idx] = sort(x);
    y = y(idx);
end
p = []; 
q = []; 
r = [];
while length(p) < n
    k = find(y == max(y),1,'last');
    a = y(k);
    b = x(k);
    id = (y>0)&(y<a);
    if ~any(id)
        break
    end
    c = mean(abs(x(id)-b)./sqrt(log(a./y(id))))/(2*n-length(p));
    % 0<y(id)<a , length(p) < 2n ==> c > 0.
    y = y - a*exp(-((x-b)/c).^2);
    p = [p b]; %#ok<AGROW>
    q = [q a]; %#ok<AGROW>
    r = [r c]; %#ok<AGROW>
end
if length(p) < n
    % Unable to find a full set of starting points.
    
    % The number of values found is
    numFound = length( p );
    % The number of points to append is
    numToAdd = n - numFound;
    
    % Choose centers by equally spacing over the domain of the data
    b = linspace( min( x ), max( x ), numToAdd+2 );
    p((numFound+1):n) = b(2:(end-1)); 
    
    % Choose the heights to be the maximum of the data
    q((numFound+1):n) = max( y );  
    
    % Choose to spreads to be the same as the spacing between the centers.
    r((numFound+1):n) = (b(2)-b(1)); 
end
start = zeros(3*n,1);
start(1:3:3*n) = q;
start(2:3:3*n) = p;
start(3:3:3*n) = r;
end  % gaussnstart



%---------------------------------------------------------
%  Utility Functions
%---------------------------------------------------------
function [x,y] = getxygrid(x,y)
% Determining the range of x values.
lengthx = length(x);

% Checking number of data points to be > 2
if lengthx < 2
    iTwoDataPointsRequiredError();
end

% Sorting data points to be in order of increasing x.
diffx = diff(x);
if any(diffx < 0)
    [x,idx] = sort(x);
    y = y(idx);
    diffx = diff(x);
end

% To avoid dividing by zero, we will get rid of repeated x entries.
tol = eps^0.7;
idx = [(diffx < tol); false];
idx2 = [false ; idx(1:end-1)];
x(idx) = (x(idx) + x(idx2)) / 2;
x(idx2) = [];
y(idx) = (y(idx) + y(idx2)) / 2;
y(idx2) = [];
lengthx = length(x);

% Data size too small, cannot find fit
if lengthx > 2
    % Checking to see whether the set of data points [x, y] are equally spaced

    % If idx has contains non-zero elements, then data points are scattered
    % Applying interpolation on the data points
    % if (sum(idx) > 0.0001)
    if all(abs(diffx-diffx(1)) < tol*max(diffx))
        % [newx, newy] = interpolate1(x, y);
        newx = linspace(min(x), max(x), numel(x));
        newy = interp1(x, y, newx);
        x = newx(:);
        y = newy(:);
    end
end
end  % getxygrid


%{
1787,1

1753,1
1753,2
1753,3

1857,2

1929,1
1929,2
1929,4


%}


%             
%         

