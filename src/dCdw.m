function [dC, C, Q] = dCdw(param, weight, point3d)
    x = param(1);
    y = param(2);
    z = param(3);
    a = param(4);
    b = param(5);
    g = param(6);
    L = param(7);
    
    sa = sin(a); ca = cos(a);
    sb = sin(b); cb = cos(b);
    sg = sin(g); cg = cos(g);
    
    vi = [ca*cb;
             sa*cb;
                 -sb];
    vj = [ca*sb*sg-sa*cg;
            sa*sb*sg+ca*cg;
                              cb*sg];
    vk = [ca*sb*cg+sa*sg;
               sa*sb*cg-ca*sg;
                                cb*cg];
    
    A = [x;y;z];
    Qs = cell(8, 1);
    Qs{1} = A;
    Qs{2} = A + L*vi;
    Qs{3} = A + L*vi + L*vk;
    Qs{4} = A + L*vk;
    Qs{5} = A + L*vj;
    Qs{6} = A + L*vi + L*vj;
    Qs{7} = A + L*vi + L*vj + L*vk;
    Qs{8} = A + L*vj + L*vk;
    
    dA = [eye(3,3), zeros(3,4)];
    dvi = [0, 0, 0, -sa*cb,  -ca*sb, 0;
               0, 0, 0,   ca*cb,  -sa*sb, 0; 
               0, 0, 0,           0,       -cb, 0];
    dviL = [L*dvi, vi];
    
    dvj =  [0, 0, 0, -sa*sb*sg-ca*cg, ca*cb*sg, ca*sb*cg+sa*sg;
                 0, 0, 0,   ca*sb*sg-sa*cg, sa*cb*sg,  sa*sb*cg-ca*sg;
                 0, 0, 0,                            0,     -sb*sg,                   cb*cg];
    dvjL = [L*dvj, vj];
    
    dvk = [0, 0, 0, -sa*sb*cg+ca*sg, ca*cb*cg, -ca*sb*sg+sa*cg;
                0, 0, 0,   ca*sb*cg+ca*sg, sa*cb*cg,  -sa*sb*sg-ca*cg;
                0, 0, 0,                              0,     -sb*cg,                   -cb*sg];
    dvkL = [L*dvk, vk];
    
    dQs = cell(8, 1);
    dQs{1} = dA;
    dQs{2} = dA + dviL;
    dQs{3} = dA + dviL + dvkL;
    dQs{4} = dA + dvkL;
    dQs{5} = dA + dvjL;
    dQs{6} = dA + dviL+ dvjL;
    dQs{7} = dA + dviL+ dvjL + dvkL;
    dQs{8} = dA + dvjL + dvkL;
    
    C = 0; dC = zeros(1, 7); Q = zeros(3, 8);
    for k = 1:1:8
        Qk = Qs{k};
        Q(:, k) = Qk;
        Gk = point3d(:, k);
        dQk = dQs{k};
        rho_k = weight(k);
        C = C + 0.5 * rho_k * (Qk - Gk)' * (Qk - Gk);
        dC = dC + rho_k * (Qk - Gk)' * dQk;
    end
end