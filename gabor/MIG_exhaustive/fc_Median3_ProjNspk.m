function [x_m3, y_m3] = fc_Median3_ProjNspk(x, y);
x_m3 = [];
y_m3 = [];

xy = [x, y];
xy_acs = sortrows(xy);
y = edbMedianSmooth(xy_acs(:,2),3);
xy_med = [xy_acs(:,1), y(:)];

x_m3 = xy_med(:,1);
y_m3 = xy_med(:,2);
