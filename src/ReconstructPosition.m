function point3d = ReconstructPosition(CamProj1, CamProj2, p1, p2)
%ReconstructPosition 双目重建
%   input
%       p1，相机1上的点坐标
%       p2，相机2上的点坐标

	lineActual = Point2Matrix(p1) * CamProj1;
    lineDesire = Point2Matrix(p2) * CamProj2;
    
    A = [lineActual(:, 1:3); % 射线C1-P1
         lineDesire(:, 1:3)]; % 射线C2-P2
    b = [-lineActual(:, 4);
         -lineDesire(:, 4)];
     
    point3d = (A\b);
end

function mat = Point2Matrix(pt2d)
    mat = [  0       -1        pt2d(2)
             1        0       -pt2d(1)
            -pt2d(2)  pt2d(1)  0];
end

