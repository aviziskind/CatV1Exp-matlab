clear all
dim = 2;

N = 500;
x = 3*randn(N,1);
y = randn(N,1);
z = 10*randn(N,1);

nStd = 2;
% theta_orig = deg2rad(45);
% phi_orig   = deg2rad(30);
theta_orig = rand*2*pi;
phi_orig   = rand*pi;

M = [rand(dim,1)*5]-5

% rotationMatrix3D
figure(1); clf; 
if dim == 2
    X = [x, y]';
    X = [rotationMatrix(theta_orig)*X]';
    x = X(:,1); y = X(:,2);    
    
    plot(x,y, 'b.');

elseif dim == 3
    X = [x, y, z]';
    R1 = rotationMatrix3D(theta_orig, 'z');
    R2 = rotationMatrix3D(phi_orig,   'x');
    X = [R2 * R1 * X]';
    x = X(:,1); y = X(:,2); z = X(:,3);

    plot3(x,y, z, 'b.');
    zlabel('z');

end
xlabel('x'); ylabel('y'); 
axis equal;
hold on;
M = mean(X,1)';
C = cov(X,1);
[V,D] = eig(C);

switchEigs = false;
if switchEigs
   V = fliplr(V);
   D = fliplr(flipud(D));    
end

3;

if dim == 2
    
%     r = det(C);
    a = sqrt( D(1,1) );
    b = sqrt( D(2,2) );
    plot([0 V(1,1)*a], [0 V(2,1)*a], 'ro-');
    plot([0 V(2,1)*b], [0 V(2,2)*b], 'go-');

%     theta = atan2(V(2,1), V(1,1));
%     R = rotationMatrix(theta);    
            
    th = linspace(0, 2*pi, 50);

    v_elps = V * [a* cos(th); b*sin(th)] * nStd;
%     v_elps(1,:) = v_elps(1,:) * a*nStd;
%     v_elps(2,:) = v_elps(2,:) * b*nStd;
%     v_elps = R * v_elps;    
    v_elps = bsxfun(@plus, M, v_elps);
    
    el_x = v_elps(1,:);
    el_y = v_elps(2,:);

    plot(el_x, el_y, 'g.-');
    
    
elseif dim == 3
%     [X2, Y2, Z2] = sphere(30);
    [X,Y,Z] = sphere(30);

    a = sqrt( D(1,1) );
    b = sqrt( D(2,2) );
    c = sqrt( D(3,3) );
    plot3([0 V(1,1)*a], [0 V(2,1)*a], [0 V(3,1)*a], 'ro-', 'linewidth', 3);
    plot3([0 V(1,2)*b], [0 V(2,2)*b], [0 V(3,2)*b], 'go-', 'linewidth', 3);
    plot3([0 V(1,3)*c], [0 V(2,3)*c], [0 V(3,3)*c], 'ko-', 'linewidth', 3);

    [theta, phi] = cart2sph(V(1,1), V(2,1), V(3,1));
%     rad2deg(theta)
%     R1 = rotationMatrix3D(theta, 'z');
%     R2 = rotationMatrix3D(-phi, 'y');
    
    v_elps = [a*X(:)';b*Y(:)';c*Z(:)']*nStd;
    v_elps = V * v_elps;    
    v_elps = bsxfun(@plus, M, v_elps);
    
    el_X = reshape(v_elps(1,:), size(X));
    el_Y = reshape(v_elps(2,:), size(X));
    el_Z = reshape(v_elps(3,:), size(X));

    h = surf(el_X, el_Y, el_Z);
    set(h, 'faceColor', 'r', 'faceAlpha', .15, 'edgeAlpha', .8, 'edgecolor', 'r');
%     Xc = X(:); Y = Y(:); Z = Z(:);

    axis equal
end
% ellipse(xc, xr, yc, yr);


