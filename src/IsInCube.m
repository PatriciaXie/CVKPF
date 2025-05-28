function inFlag = IsInCube(point3d, cube)
%IsInCube 判断三维点是否在盒子内
%   input
%       cube，盒子
%       point3d，三维坐标
%       epsilon，误差

% global cube
    inFlag = 1;
    epsilon = norm(cube.d3Vertices(:, 1) - cube.d3Vertices(:, 2)) - 5; % 靠近边缘的5mm不要
    for p=1:6
        distance = (point3d-cube.d3PlaneCenters(:,p))' * cube.d3Norms(:,p);
        if ( distance > epsilon )
            inFlag = 0;
            return;
        end
    end
%     figure(2), clf
%     hold on
%     [Xb, Yb, Zb] = cube2line(cube);
%     plot3(Xb, Yb, Zb, 'k');
%     plot3(point3d(1), point3d(2), point3d(3), 'ro');
%     p = 1; plot3(cube.d3PlaneCenters(1,p), cube.d3PlaneCenters(2,p), cube.d3PlaneCenters(3,p),'bo');
%     plot3(cube.d3PlaneCenters(1,p)+200*[0, cube.d3Norms(1,p)], cube.d3PlaneCenters(2,p) + 200*[0, cube.d3Norms(2,p)], cube.d3PlaneCenters(3,p) + 200*[0, cube.d3Norms(3,p)],'g-');
%     axis equal
%     hold off
end
