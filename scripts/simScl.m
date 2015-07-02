function D = simScl(x1, x2, gam)
    % gam = 0 : insensitive to any scaling
    % gam = 1 : sensitive to differential scaling only.

    normx1 = norm(x1);
    normx2 = norm(x2);
    if (normx1 == 0) || (normx2 == 0)
        C = 0;
    else
        C = dot(x1, x2)/(normx1 * normx2);
    end

    D = C - gam/2 * [normx1 /normx2 + normx2/normx1]; 

end