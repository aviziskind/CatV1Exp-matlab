% projection in 2D onto a line
do_2Dline = 0;
do_3Dline = 0;
do_3D_normalPlane = 0;
do_3D_normalPlane_view = 0;
do_3D_spannedPlane = 0;
do_realTest = 1;

N = 10;

if do_2Dline
    D = 2;
    figure(4); clf;
    
    X = randn(D,N);
    plot(X(1,:), X(2,:), '.')
    axis([-3, 3, -3, 3])

    theta = pi/3;
    v = [cos(theta); sin(theta)];
    A = v*v';
    
%     R = rotationMatrix(theta);
%     P = [1 0; 0 0];
%     A = R*P*R';
    Y = A*X;

    hold on;
    plot(Y(1,:), Y(2,:), 'g.')

    [X1, Y1] = linesFromAtoB(X(1,:),X(2,:), Y(1,:),Y(2,:));
    plot(X1, Y1, 'ro-')

end


if do_3Dline
    D = 3;
    
    figure(4); clf;
    X = randn(D,N);
    
    plot3(X(1,:), X(2,:), X(3,:), 'b.');   hold on;  
    xlabel('x'); ylabel('y'); zlabel('z');
    axis([-3, 3, -3, 3, -3, 3])
    
    v = [1 .3 3]'; 
    plot3([0, v(1)], [0, v(2)], [0, v(3)], 'k*-')
    v = v/norm(v);
    
    A = v*v';    
    Y = A*X;

    plot3(Y(1,:), Y(2,:), Y(3,:), 'g.');

    [X1, Y1, Z1] = linesFromAtoB(X(1,:),X(2,:), X(3,:), Y(1,:),Y(2,:), Y(3,:));
    plot3(X1, Y1, Z1, 'ro-')

end

% compAinDirOfB = @(A, B) ( dot(B,A)/dot(B,B) ) * B;
        

if do_3D_normalPlane
    figure(5); clf;
        
    X = randn(D,N);
    plot3(X(1,:), X(2,:), X(3,:), 'b.');    
    xlabel('x'); ylabel('y'); zlabel('z');
    axis([-3, 3, -3, 3, -3, 3])

    v = [1 .3 3]';
    hold on;
    zeroaxes;
    plot3([0, v(1)], [0, v(2)], [0, v(3)], 'k*-')
        

    P = eye(3); P(1,1) = 0;    
    Q = GramSchmidt(v);
    A = P*Q';    
    Y = A*X;
    
%     A = -cross2mtx(v)*cross2mtx(v)/norm(v)^2;
%     Y = A*X;
    
    plot3(Y(1,:), Y(2,:), Y(3,:), 'g.');

    [X1, Y1, Z1] = linesFromAtoB(X(1,:),X(2,:),X(3,:),  Y(1,:),Y(2,:),Y(3,:));
    plot3(X1, Y1, Z1, 'ro-')

    f = @(x,y) x*(-v(1)/v(3)) +y*(-v(2)/v(3));
    h = fmesh(f, [-3:.2:3], [-3:.2: 3]);
    set(h, 'edgecolor', [.7 1 .7]);
    

end



if do_3D_normalPlane_view
    figure(5); clf; hold on;
    figure(6); clf; hold on;
    
    X = randn(D,N);
    figure(5);
    plot3(X(1,:), X(2,:), X(3,:), 'b.');    
    xlabel('x'); ylabel('y'); zlabel('z');
    axis([-3, 3, -3, 3, -3, 3])
    zeroaxes;
    h3_1 = plot3(0, 0, 0, 'ro-');
    h3_2 = plot3(0, 0, 0, 'go-');
    h3_3 = plot3(0, 0, 0, 'ko-');
    
    figure(6); 
    h2d = plot(0,0, '.');
    
    while true
        yn = input('Rotate to the desired view, then press enter.');
        if strcmpi(yn, 'q')
            return;
        end        
        [az, el] = view;
        v_r = [deg2rad(az)-pi/2, deg2rad(el), 1];
        [v_x, v_y, v_z] = sph2cart(v_r(1), v_r(2), v_r(3));
        set(h3_1, 'xdata', [0 v_x], 'ydata', [0 v_y], 'zdata', [0 v_z]);            
        
        Q = GramSchmidt([v_x; v_y; v_z]);
        set(h3_2, 'xdata', [0 Q(1,2)], 'ydata', [0 Q(2,2)], 'zdata', [0 Q(3,2)]);            
        set(h3_3, 'xdata', [0 Q(1,3)], 'ydata', [0 Q(2,3)], 'zdata', [0 Q(3,3)]);            
                
        P = eye(3); P(1,1) = 0;    
        A = P*Q';    
        Y = A*X;
        
        set(h2d, 'xdata', Y(2,:), 'ydata', Y(3,:));
    end
    
end



if do_3D_spannedPlane
        
%     X = randn(D,N);
    v1 = [1 .5 .2]';
    Q = GramSchmidt(v1);
    v2 = Q(:,2);
    v3 = Q(:,3);

    P = eye(3); P(3,3) = 0;    
    A = Q*P*Q';    
    Y = A*X;
    
    A2 = P*Q';
    X2 = A2*X;
    
    figure(5); clf;
    plot3(X(1,:), X(2,:), X(3,:), 'b.');  hold on;
    xlabel('x'); ylabel('y'); zlabel('z');
    axis([-3, 3, -3, 3, -3, 3])
    
    zeroaxes;
    plot3([0, v1(1)], [0, v1(2)], [0, v1(3)], 'ro-')
    plot3([0, v2(1)], [0, v2(2)], [0, v2(3)], 'go-')
            
    plot3(Y(1,:), Y(2,:), Y(3,:), 'g.');

    [X1, Y1, Z1] = linesFromAtoB(X(1,:),X(2,:),X(3,:),  Y(1,:),Y(2,:),Y(3,:));
    plot3(X1, Y1, Z1, 'ro-')

    3;
    f = @(x,y) x*(-v3(1)/v3(3)) +y*(-v3(2)/v3(3));
    h = fmesh(f, [-3:.2:3], [-3:.2: 3]);
    set(h, 'edgecolor', [.7 1 .7]);
    
    figure(6); clf;
    plot(X2(1,:), X2(2,:), '.'); axis equal square; 
    axis([-2 2, -2 2]);

end


if do_realTest
    
    %xy plane -  circle
    %yz plane - square
    N = 100;
    th = linspace(0, 10*2*pi, N);
    x = cos(th); y = sin(th);
    z = linspace(-.5, .5, N);
    X = [x(:), y(:), z(:)]';
    Q = GramSchmidt(randn(3));
    v1 = Q(:,1); v2 = Q(:,2); v3 = Q(:,3);
    Y = Q*X;
    
    figure(5); clf;
    plot3(Y(1,:), Y(2,:), Y(3,:), 'b.-');  hold on; axis equal
    xlabel('x'); ylabel('y'); zlabel('z');
    plot3([0, v1(1)], [0, v1(2)], [0, v1(3)], 'ro-')
    plot3([0, v2(1)], [0, v2(2)], [0, v2(3)], 'go-')
    plot3([0, v3(1)], [0, v3(2)], [0, v3(3)], 'ko-')
%     axis([-2, 2, -2, 2, -2, 2])


    figure(6); clf;
    P = eye(3); P(3,3) = 0;
    A = P*Q';   % projects to plane spanned by v1, v2, in e1, e2 basis.
%     A2 = v1*v1' + v2*v2';   % this projects to the plane spanned by v1 
        % and v2, and keeps in the same coords.
%     
    3;
    Z = A*Y;   % = P*Q'*Q*X
    
    plot(Z(1,:), Z(2,:), '.');

    
    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    