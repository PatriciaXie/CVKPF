function inFlag = IsInCube(point3d, cube)
%IsInCube �ж���ά���Ƿ��ں�����
%   input
%       cube������
%       point3d����ά����
%       epsilon�����

% global cube
    inFlag = 1;
    epsilon = norm(cube.d3Vertices(:, 1) - cube.d3Vertices(:, 2)) - 5; % ������Ե��5mm��Ҫ
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
