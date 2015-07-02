function x_max = gradientAscent(maxF, gradF, x0, Xs, maxOfF, tolMax)
    % maxF is the maximizing function (the function which we want to maximize)
    % gradF is the gradient of maxF.
    % x0 is the initial value of maxF, hopefully close to the max
    % Xs should be a cell containing the dimensions over which maxF is to
    %       be searched (eg {x}, or {x, y})
    % the function returns x_max, where maxF(x_max) is the (hopefully 
    %    global) maximum of maxF.

    global actualGaborParams;
    Np = length(x0);
    pXLabels = {'A', 'mx', 'sx', 'my', 'sy', 'k', 'p', 't', 'c'};
%     pLabels = {'A', '\mu_x', '\sigma_x', '\mu_y', '\sigma_y', 'k', '\phi', '\theta', 'c'};

    gradient_fig = 23; 
    params_fig = 24;
    prevF = maxF(x0);
    maxIterations = 4;
    numPointsPerIteration = 5;
    [xs, ys] = elements(Xs);  % assume 2 dimensions: x & y 
    [xs_grid, ys_grid] = meshgrid(xs, ys);
    
    showProgress = true;
    
    function xbnd = xbounded(x, lower, upper)
        xbnd = x;
        if xbnd < lower
            xbnd = lower;
        elseif xbnd > upper
            xbnd = upper;
        end
    end

    
    function c = findBestC(x, dx, searchRange, iteration_id)

        cs = linspace(searchRange(1), searchRange(2), numPointsPerIteration);
        ysTrials = zeros(size(cs));
        for ti = 1:length(cs)
            progressBar;
            ind = find( csResults(:,1) == cs(ti) );
            if ~isempty(ind)        % already computed for this value of c (no need to do again)
                ysTrials(ti) = csResults(ind,2);
            else
                ysTrials(ti) = maxF( x + cs(ti) * dx );
                csResults = [csResults; cs(ti), ysTrials(ti)]; %#ok<AGROW>
            end
            ysTrials(ti) = maxF( x + cs(ti) * dx );
        end
        [c_max, c_ind] = max(ysTrials);
        
                if showProgress
                    figure(22); 
                    plot(cs, ysTrials, ['o-' color(iteration_id)]);
                end
        
        if (iteration_id == maxIterations)
            c = c_max;
        else
            nextRangeInds = [ xbounded(c_ind-1, 1, length(cs)),  xbounded(c_ind+1, 1, length(cs)) ];
            c = findBestC(x, dx, cs(nextRangeInds), iteration_id+1);
        end     
        
    end

    hrz_fig = 3;
    count = 1;
    show = 1;
    initialSearchRange = [0, 10];
    x = x0;
    
    if showProgress
        figure(params_fig); clf; hold on;
        xlim([0, Np+1]);
        drawHorizontalLine(1);
        set(gca, 'XTick', 1:Np, 'XTickLabel', pXLabels(1:Np))
    end

    tol = 0.05;
    thisF = prevF*2;
    
    while abs(thisF-prevF)/prevF > tol
        dx = gradF(x);
        progressBar('init-', maxIterations * numPointsPerIteration, 20);
        csResults = zeros(0,2);
                if showProgress
                    figure(22); clf; hold on;
                    xlim(initialSearchRange);
                    drawHorizontalLine(maxOfF);
                end

%         c = 1;
        c = findBestC(x, dx, initialSearchRange, 1);
        progressBar('done');

                if showProgress
                    figure(gradient_fig);
                    subplot(1,hrz_fig, show);
                    graySquareImage(xs, ys, gabor(x, {xs_grid, ys_grid}) , num2str(count));

                    figure(params_fig);
                    if exist('actualGaborParams', 'var') && ~isempty(actualGaborParams)
                        plot(1:Np, (x+c*dx+eps) ./ (actualGaborParams + eps), ['o' color(count)]); 
                    else
                        plot(1:Np, (x+c*dx+eps) ./ (x0 + eps), ['o' color(count)]); 
                    end                    
                    show = cycle(show, [1:hrz_fig]);
                end
        
        x = x + c*dx;

        prevF = thisF;
        thisF = maxF(x);
        
        figure(1);
        count = count + 1;
        drawnow;
        
    end

    x_max = x;
end