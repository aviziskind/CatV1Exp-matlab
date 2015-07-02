%%
tb_X = 15;
tb_Y = 15;
screenW = 30;

r_in = 50;
r_out = 80;

doLeft = 1;
doRight = 0;

if doRight
     CL = -screenW/2;
    CR = screenW/2;
    
else
    
    CL = 0;
    CR = screenW;
   
end


yL_out = @(x) -sqrt( rectified( r_out.^2 - (x-CL).^2 ));
yL_in  = @(x) -sqrt( rectified( r_in.^2 - (x-CL).^2 ));
yR_out = @(x) -sqrt( rectified( r_out.^2 - (x-CR).^2 ));
yR_in  = @(x) -sqrt( rectified( r_in.^2 - (x-CR).^2 ));

y_hi = @(x) min(yL_out(x), yR_out(x));
y_lo = @(x) max(yL_in(x), yR_in(x));


xlims_L = CL + r_out*[-1, 1];
xlims_R = CR + r_out*[-1, 1];
if doRight
    xlims = lims([xlims_L, xlims_R]);
else
    xlims = xlims_L;
end

figure(86); clf; hold on;
% xlim1 = -(screenW/2) - r_out;
% xlim2 = (screenW/2) + r_out;
% xlims = [xlim1, xlim2];

h = fplot(yL_out, xlims, 'b');
fplot(yL_in, xlims, 'b:')
if doRight
    fplot(yR_out, xlims, 'g')
    fplot(yR_in, xlims, 'g:')
    fplot(y_hi, xlims, 'r')
    fplot(y_lo, xlims, 'r:')
end

xlim(lims(xlims, .05));
axis equal
%%

qLR = quad(@(x) y_hi(x)-y_lo(x), xlims(1), xlims(2));
qR = quad(@(x) yR_out(x) - yR_in(x), xlims(1), xlims(2))
qL = quad(@(x) yL_out(x) - yL_in(x), xlims(1), xlims(2))