function [X, Y] = getHyperplane(s, xlims, ylims)

    x0 = s(1); y0 = s(2);
    M_s = y0/x0;
%     C_s = y0 + 1/M_s*x0;

    M_h = -1/M_s;
    C_h = y0 - M_h*x0;
    
    xL = xlims(1);
    xR = xlims(2);
    
    yB = ylims(1);
    yT = ylims(2);
    
    y_xL = M_h*xL + C_h;
    y_xR = M_h*xR + C_h;

    x_yB = (yB-C_h)/M_h;
    x_yT = (yT-C_h)/M_h;

    pts = [xL, y_xL;
           xR, y_xR;
           x_yB, yB;
           x_yT, yT]';
       %%
    pt_ok = false(1,4);
    for i = 1:4
        pt_ok(i) = ibetween(pts(1,i), xlims) & ibetween(pts(2,i), ylims);
    end
    assert(nnz(pt_ok) == 2);
    X = pts(1, pt_ok);
    Y = pts(2, pt_ok);
    

end