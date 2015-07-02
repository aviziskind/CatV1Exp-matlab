function v = fc_FriedmanStuedzleMerit(x, y);

xy = [x, y];
xy_acs = sortrows(xy);
% --- (1) run 3-pont median to get rid of outliers ---
y = edbMedianSmooth(xy_acs(:,2),3);
xy_med = [xy_acs(:,1), y(:)];

% --- (2) locally linear fit with constant bandwidth --- 
local_vars = edbLocallyLinearFit(x, y, nbin);
edbMovingAver(x, wnd);
