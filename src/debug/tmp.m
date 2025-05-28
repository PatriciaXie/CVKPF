% MATLAB ¼ÆËãÐý×ª¾ØÕó
syms a b c
Mx = [1, 0, 0;
            0, cos(a), -sin(a);
            0, sin(a), cos(a)];
My = [cos(b), 0, sin(b);
            0, 1, 0;
            -sin(b), 0, cos(b)];
Mz = [cos(c), -sin(c), 0;
            sin(c), cos(c), 0;
            0, 0, 1];
        
M = Mz * My * Mx;