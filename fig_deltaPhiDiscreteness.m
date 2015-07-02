
nPh = 4;


figure(232); clf; hold on;
r = .8;
r_txt = r*1.5;
r_im = r*[2.3, 2.5]; shft = [-.1 0];
ph_deg = linspace(0, 360, nPh+1); ph_deg = ph_deg(1:nPh);
phPos_deg = 90-ph_deg;
ph_rad = deg2rad(phPos_deg);
drawCircle(r);
mainAx = gca;
quiver(zeros(1,nPh), zeros(1,nPh), r*cos(ph_rad), r*sin(ph_rad));
mainPh = 90;
set(mainAx, 'xtick', [], 'ytick', []); box on;
title([num2str(nPh) ' phases']);
for i = 1:nPh
    
    boldflag = iff(ph_deg(i) == mainPh, '\bf', '');
    d = circDist(ph_deg(i), mainPh, 360);
    txt = {['\fontsize{13}' boldflag num2str(ph_deg(i)) '\rm\circ'], ['\fontsize{9}(' num2str(d) ') ']};
    
    text( r_txt*cos(ph_rad(i)), r_txt*sin(ph_rad(i)), txt, 'horizontalalignment', 'center');
end
axis equal;
L = 2.5;
axis([-L+shft(1), L, -L, L]);

ori = 90;
sp = 32;
hs = zeros(1,nPh);
for i = 1:nPh    
    frm = generateGratingFrame([32 32], ori, 32, 180+rad2deg(ph_rad(i)));
    axes(mainAx);
    [x,y] = dsxy2figxy(r_im(1)*cos(ph_rad(i))+shft(1), r_im(2)*sin(ph_rad(i))+shft(2));
    hs(i) = axes('Position', [x,y, .1, .1]);
    set(hs(i), 'units', 'pixels');
    p = get(hs(i), 'position'); [xc, yc] = deal(p(1), p(2));
    set(hs(i), 'position', [xc-16, yc-16, 32, 32]);
    imagesc(frm); colormap('gray'); set(hs(i), 'xtick', [], 'ytick', [])
    
%     text( r_im*cos(ph_rad(i)), r_txt*sin(ph_rad(i)), [num2str(ph_deg(i)) '\circ'], 'horizontalalignment', 'center');
end
for i = 1:nPh
    axes(hs(i))
end