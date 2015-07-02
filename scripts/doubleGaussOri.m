function [y, y1, y2] = doubleGaussOri(A1, A2, sigma, B, x, x_cent1)
%     y = A .* exp( -(x - dir_pref_deg).^2 ./(2*(sigma).^2)) + B;
    y = zeros(size(x));    
    
    if nargin < 6
        x_cent1 = 0;
    end
    
    x_cent2 = mod(x_cent1+180, 360);
    
    dist_0   = circDist(x, x_cent1, 360);    
    dist_180 = circDist(x, x_cent2, 360);

    y1 = A1 .* exp( -( dist_0   ).^2 ./(2*(sigma).^2)) + B;
    y2 = A2 .* exp( -( dist_180 ).^2 ./(2*(sigma).^2)) + B;
    
    y = y1+y2;
    
%     pref_idx = abs(x) <= 90;        
%     x_opp_dist = 180 - abs( x(~pref_idx));
%     y(pref_idx)  = A  .* exp( -( x(pref_idx) ).^2 ./(2*(sigma).^2)) + B;
%     y(~pref_idx) = A2 .* exp( -( x_opp_dist  ).^2 ./(2*(sigma).^2)) + B;
            
    
end


%     x_dist = x;
%     x_gt_90 = find(x > 90);
%     x_dist(x_gt_90) = 180 - abs( x_dist(x_gt_90));    
%     x_lt_90 = find(x < -90);
%     x_dist(x_lt_90) = 180 - abs( x_dist(x_lt_90) );
