function [xlims, ylims] = getEffectiveSpaceOfGabor(gaborParameters, numStdDevs)

    if nargin < 2
        numStdDevs = 3;
    end

%     global Xc;
	[A, mu_x, mu_y, sig_x, sig_y, theta, k, phi] = elements(gaborParameters);

    rotateAndShift = @(X, theta, Xcenter)  rotationMatrix(theta) * X - repmat(Xcenter, 1,size(X,2));

    X_center = [mu_x; mu_y];

    outerLimitsInGaborFrame = [ [ numStdDevs * sig_x; 0], ...
        [-numStdDevs * sig_x; 0], ...
        [0;  numStdDevs * sig_y], ...
        [0; -numStdDevs * sig_y] ];

    outerLimitsRealFrame = rotateAndShift(outerLimitsInGaborFrame, theta, -(X_center) );

    xlims = [  min( outerLimitsRealFrame(1,:) ),  max( outerLimitsRealFrame(1,:) )];
    ylims = [  min( outerLimitsRealFrame(2,:) ),  max( outerLimitsRealFrame(2,:) )];

    if (nargout == 1)
        xlims = [xlims ylims];
    end
    
%     fmesh(gaborFunction(p), lims(1:2), lims(3:4), 'calc', 'group');
%     hold on
%     stem3(X_center(1), X_center(2), 1);
%     stem3(outerLimitsRealFrame(1,1:2), outerLimitsRealFrame(2,1:2), .5*ones(2,1), 'b');
%     stem3(outerLimitsRealFrame(1,3:4), outerLimitsRealFrame(2,3:4), .5*ones(2,1), 'g');
%     hold off;

end
