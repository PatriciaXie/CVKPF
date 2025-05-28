function showScenario(stereoModel, cube)

% 1 相机模型
R1 = angle2dcm(0, 0, 0);
T1 = stereoModel.cams(1).center;
C1 = CentralCamera('focal', stereoModel.cams(1).intrinsic(1,1), 'pixel', 1e-6, ...
    'resolution', [2048 2048], 'centre', [stereoModel.cams(1).intrinsic(1,3) stereoModel.cams(1).intrinsic(2,1)], ...
    'distortion', [0, 0, 0, 0, 0],...
    'pose', SE3(R1, T1), 'name', 'C1', 'color', [0, 0, 0]);

tmp = inv(stereoModel.cams(2).intrinsic) * stereoModel.cams(2).projection;
R2 = tmp(:, 1:3);
T2 = stereoModel.cams(2).center;
C2 = CentralCamera('focal', stereoModel.cams(2).intrinsic(1,1), 'pixel', 1e-6, ...
    'resolution', [2048 2048], 'centre', [stereoModel.cams(2).intrinsic(1,3) stereoModel.cams(2).intrinsic(2,1)], ...
    'distortion', [0, 0, 0, 0, 0],...
    'pose', SE3(R2, T2), 'name', 'C2', 'color', [0, 0, 0]);

[Xb, Yb, Zb] = cube2line(cube);
% 2 相机摆放示意图
figure(1);
C1.plot_camera('scale', 100, 'color', 'r', 'label');
C2.plot_camera('scale', 100, 'color', 'b', 'label');
view(0, 0);
hold on
plot3(Xb, Yb, Zb, 'k');
axis([-2000, 500, -1000, 1000, -500, 2000]);
axis equal 
grid on
hold off
end