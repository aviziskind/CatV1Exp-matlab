function y = gabor(gaborParameters, X, gradientFlag) %#ok<INUSD>
    % gaborParameters is a vector of 9 parameters:  (A, mu_x, sig_x, mu_y, sig_y, k, phi, theta, const)
    % X can be:
        % - a single co-ordinate [x;y]
        % - a list of co-oordinates [x1, x2, ... xn; y1, y2, ... yn]
        %        (in columns or rows)
        % - a cell containing a matrix of x and y co-ordinates {X, Y}
    % providing a third parameter ('gradientFlag') makes the function return the
    %    gradient of the gabor function instead of just the gabor function.
    constAllowed = true;
    
    global Xc;
    if isempty(Xc)
        Xc = [0;0];
    end;
    if constAllowed && (length(gaborParameters) == 9)
        const = gaborParameters(9);
    else
        gaborParameters = gaborParameters(1:8);
        const = 0;        
    end
	[A, mu_x, sig_x, mu_y, sig_y, k, phi, theta] = elements(gaborParameters);
    
    gaborCenter = [mu_x; mu_y];
    expsqr = @(x, mu, sigma)  exp( -((x-mu).^2)./(2*sigma^2));
    shiftAndRotateToGaborFrame = @(X)  rotationMatrix(-theta) * (X - repmat(gaborCenter+Xc, 1,size(X,2)));
    
    % GABOR FUNCTION        
    function G = gaborFunction(XY)
        XYp = shiftAndRotateToGaborFrame(XY);
        [Xp, Yp] = rowsOf(XYp);
        G = A * expsqr(Xp, 0, sig_x) .* expsqr(Yp, 0, sig_y) .* cos( k * Xp + phi ) + const;    
    end

    % GABOR GRADIENT FUNCTION
    function gradG = gaborGradient(XY)
        XYp = shiftAndRotateToGaborFrame(XY);
        [X,  Y]  = rowsOf(XY);
        [Xp, Yp] = rowsOf(XYp);
        
        exp_Xp = expsqr(Xp, 0, sig_x);%  exp( -(Xp.^2)./(2*sig_x^2));
        exp_Yp = expsqr(Yp, 0, sig_y);%  exp( -(Yp.^2)./(2*sig_y^2));
        cos_Xp = cos( k*Xp + phi );
        sin_Xp = sin( k*Xp + phi );
        ex_ey_cx = exp_Xp .* exp_Yp .* cos_Xp;
        ex_ey_sx = exp_Xp .* exp_Yp .* sin_Xp;
        
        dG_dxp = A * (ex_ey_cx .* (-Xp/(sig_x^2)) - ex_ey_sx * k);
        dG_dyp = A *  ex_ey_cx .* (-Yp/(sig_y^2)) ;
        
        dxp_dmu_x = -cos(theta);
        dyp_dmu_x =  sin(theta);        
        dxp_dmu_y = -sin(theta);
        dyp_dmu_y = -cos(theta);
        dxp_dtheta = -(X-mu_x)*sin(theta) + (Y-mu_y)*cos(theta);
        dyp_dtheta = -(X-mu_x)*cos(theta) - (Y-mu_y)*sin(theta);
        
        dG_dA     = ex_ey_cx;
        dG_dmu_x  = dG_dxp .* dxp_dmu_x + dG_dyp .* dyp_dmu_x;
        dG_dmu_y  = dG_dxp .* dxp_dmu_y + dG_dyp .* dyp_dmu_y;
        dG_dsig_x = A * ex_ey_cx .* (Xp.^2 / sig_x^3);
        dG_dsig_y = A * ex_ey_cx .* (Yp.^2 / sig_y^3);
        dG_dk     = A * ex_ey_sx .* (-Xp);
        dG_dphi   = A * ex_ey_sx * (-1);
        dG_dtheta = dG_dxp .* dxp_dtheta + dG_dyp .* dyp_dtheta;

        if constAllowed
            dG_dconst = ones(1, size(X,2));
        else
            dG_dconst = zeros(1, size(X,2), 0);
        end
        
%         if size(dG_dA, 1) == 1
%             stackdim = 1;
%         elseif size(dG_dA, 2) == 1
%             stackdim = 2;
%         else
%             stackdim = 3;
%         end
        stackdim = 3;
        gradG = cat(stackdim, ... 
            dG_dA, dG_dmu_x, dG_dsig_x, dG_dmu_y, dG_dsig_y, dG_dk, dG_dphi, dG_dtheta, dG_dconst);
        
    end 
        

    if nargin == 2  % regular gabor function
        functionToUse = @gaborFunction; 
    elseif (nargin == 3) % calculate the 8/9-element gradient
        functionToUse = @gaborGradient;
    end
        
    
    if iscell(X)
        % X can be a cell type, with two matrices of co-ordinates;
        xs = X{1}; ys = X{2};
        % assert(all(size(xs) == size(ys)));
        [m,n] = size(xs);
        X_list = [xs(:) ys(:)]';
        y_list = functionToUse(X_list);
        
        if nargin == 2
            reshapeDims = [m,n];
        elseif nargin == 3
            reshapeDims = [m,n, length(gaborParameters)];
        end
        y = reshape(y_list, reshapeDims);
    else
        % X can be a list of co-ordinates: with each column is a set of co-ordiantes;
        % if X is in rows, then X is temporarily flipped for calculations.
        [m,n] = size(X);
        if (m == 2)
            y = functionToUse(X);
        elseif (n == 2)
            y = functionToUse(X')';
        else
            error('X is not the right shape: dimensions must be mx2 or 2xn');
        end
    end
        
end



