function [prj_sin, prj_cos] = fc_GaborProj_SC(gb_params, sin_params, x0s, y0s, p0s)

  prj_sin = [];
  prj_cos = [];

  Ny0 = length(x0s);
  Nx0 = length(y0s);
  Nph0 = length(p0s);
  
  p0s = pi*p0s/180;

  [X0s, Y0s, P0s] = meshgrid(x0s, y0s, p0s);
  
  % --- Gabor parameters ---
	x0_gb = gb_params(1);
	y0_gb = gb_params(2);
	Ori_deg_gb = gb_params(3);
  Sp_pix_gb = gb_params(4);
	Sx2 = gb_params(5);
	Sy2 = gb_params(6);
	ph0_deg_gb = gb_params(7);

  k_gb = 2*pi/Sp_pix_gb;
  ori_gb = pi*Ori_deg_gb/180;
  ph0 = pi*ph0_deg_gb/180;
  
  % --- Grating parameters ---
  x0_sin = sin_params(1);
  y0_sin = sin_params(2);
	Ori_deg_sin = sin_params(3);
	Sp_pix_sin = sin_params(4);
	ph0_deg_sin = sin_params(5);

	ori_sin = pi*Ori_deg_sin/180;
	W = 2*pi/Sp_pix_sin;
	ph0_sin = pi*ph0_deg_sin/180;
  ph0_cos = ph0_sin - pi/2;
  
  Wy = W*sin(ori_sin);
  Wx = W*cos(ori_sin);
  
  POs_sin = zeros(Ny0, Nx0, Nph0);
  POs_cos = zeros(Ny0, Nx0, Nph0);
  
  POs_sin = ph0_sin + Wx*(X0s - x0_sin) - Wy*(Y0s - y0_sin);
  POs_cos = ph0_cos + Wx*(X0s - x0_sin) - Wy*(Y0s - y0_sin);

	Cs = cos(ori_sin - ori_gb);
	Sn = sin(ori_sin - ori_gb);

	mu = k_gb * W * Cs * Sx2;
	Ep = exp( -0.5*(k_gb^2*Sx2 + W^2*(Cs^2*Sx2 + Sn^2*Sy2) + 2*mu) );
	Em = exp( -0.5*(k_gb^2*Sx2 + W^2*(Cs^2*Sx2 + Sn^2*Sy2) - 2*mu) );
  
	prj_sin = cos(P0s + POs_sin)*Ep + cos(P0s - POs_sin)*Em;
	prj_cos = cos(P0s + POs_cos)*Ep + cos(P0s - POs_cos)*Em;

