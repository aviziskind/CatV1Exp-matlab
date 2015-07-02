%%
dens_lo_mm3 = 40000;
dens_hi_mm3 = 60000;

mm_per_micron = 1e-3;



ball1_r_microns = 30;
ball2_r_microns = 50;
ball3_r_microns = 100;

nInBall = @(r_um, dens_mm3)  (4/3)*pi*r_um^3 * (mm_per_micron.^3) * dens_mm3;

nInBall1_lo = nInBall(ball1_r_microns, dens_lo_mm3);
nInBall1_hi = nInBall(ball1_r_microns, dens_hi_mm3);

nInBall2_lo = nInBall(ball2_r_microns, dens_lo_mm3);
nInBall2_hi = nInBall(ball2_r_microns, dens_hi_mm3);

nInBall3_lo = nInBall(ball3_r_microns, dens_lo_mm3);
nInBall3_hi = nInBall(ball3_r_microns, dens_hi_mm3);

round([nInBall1_lo, nInBall1_hi])
round([nInBall2_lo, nInBall2_hi])
round([nInBall3_lo, nInBall3_hi])


