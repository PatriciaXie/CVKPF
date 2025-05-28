clear; clc;
addpath(genpath('src\'));
cubePath = 'data/realdata001-1/cube.mat';
% % datasetPath = 'data\realdata001-1'; 
% % imgHandle = ImageReader(2, datasetPath, 'bmp');
% % img1 = imgHandle.Read(1, 1); 
% % img2 = imgHandle.Read(2, 1); 

%% 1 初值确定
% 用处理过的图像
point2d_1 = [628, 557; 1549, 529; 1570, 629; 837, 642; 624, 1457; 1550, 1485; 1571, 1378; 834, 1364];
point2d_2 = [753, 534; 1628, 542; 1436, 646; 745, 649; 760, 1444; 1641, 1419; 1448, 1336; 754, 1354];
cameraPath = 'data\realdata001-1\stereoModel.mat'; 
load(cameraPath);
CamProj1 = stereoModel.cams(1).projection;
CamProj2 = stereoModel.cams(2).projection;
point3d = zeros(3, 8);
for i = 1:1:8
    point1 = [point2d_1(i, 1); point2d_1(i, 2)]; 
    point2 = [point2d_2(i, 1); point2d_2(i, 2)]; 
    point3d(:, i) = ReconstructPosition(CamProj1, CamProj2, point1, point2);
end

x0 = point3d(1,1); y0 = point3d(2,1); z0 = point3d(3,1);
AB = point3d(:, 2) - point3d(:, 1);
AE = point3d(:, 5) - point3d(:, 1);
L0 = norm(AB);
ab = AB/L0;
ae = AE/norm(AE);
ad = cross(ab, ae);
ae = cross(ad, ab);
i_y = ab(2); i_z = ab(3); j_z = ae(3);
beta0 = asin(-i_z);
alpha0 = acos(i_y/cos(beta0));
gamma0 = asin(j_z/cos(beta0));
param = [x0, y0, z0, alpha0, beta0, gamma0, L0];

%% 2 梯度下降
error = 10; eta = 1e-6; k = 1;
% weight = [0.1, 1, 0.1, 0.1, 1, 0.1, 1, 1];
% weight = [1, 1, 1, 1, 1, 1, 1, 1];
weight = [0, 1, 1, 1, 1, 1, 1, 1];
while error > 1e-5
    [dC, C, Q] = dCdw(param, weight, point3d);
    param = param - eta * dC;
    error = eta * norm(dC);
%     fprintf('第%d次迭代，误差为%f，代价函数为%.2f\n', k, error, C);
    k = k + 1;
end
cube = ComputeCube(Q);
save(cubePath, 'cube');
%% 3 立体验证
cube1.d3Vertices = point3d;
% 1 相机模型
R1 = angle2dcm(0, 0, 0);
T1 = stereoModel.cams(1).center;
C1 = CentralCamera('focal', stereoModel.cams(1).intrinsic(1,1), 'pixel', 1e-6, ...
    'resolution', [2000 2000], 'centre', [stereoModel.cams(1).intrinsic(1,3) stereoModel.cams(1).intrinsic(2,1)], ...
    'distortion', [0, 0, 0, 0, 0],...
    'pose', SE3(R1, T1), 'name', 'C1', 'color', [0, 0, 0]);

tmp = inv(stereoModel.cams(2).intrinsic) * stereoModel.cams(2).projection;
R2 = tmp(:, 1:3)';
T2 = -R2 * stereoModel.cams(2).center';
C2 = CentralCamera('focal', stereoModel.cams(2).intrinsic(1,1), 'pixel', 1e-6, ...
    'resolution', [2000 2000], 'centre', [stereoModel.cams(2).intrinsic(1,3) stereoModel.cams(2).intrinsic(2,1)], ...
    'distortion', [0, 0, 0, 0, 0],...
    'pose', SE3(R2, T2), 'name', 'C2', 'color', [0, 0, 0]);

[Xb1, Yb1, Zb1] = cube2line(cube1);
[Xb, Yb, Zb] = cube2line(cube);
% 2 相机摆放示意图
figure(3), clf;
C1.plot_camera('scale', 50, 'color', 'r', 'label');
C2.plot_camera('scale', 50, 'color', 'b', 'label');
view(0, 0);
hold on
plot3(Xb1, Yb1, Zb1, 'b-');
plot3(Xb, Yb, Zb, 'g-');
axis equal
axis([-300, 300, -500, 500, -50, 1100]);
legend('Initial box', 'Fitting box')
grid on
hold off

%% 4 重投影验证
% 相机1
figure(1), clf, imshow('data\realdata001-1\cam1\fr00001.bmp');
hold on
for k = 1:1:8
    point = CamProj1 * e2h(Q(:, k)); point = point(1:2)/point(3); % 最优值重投影，绿色
    point0 = CamProj1 * e2h(point3d(:, k)); point0 = point0(1:2)/point0(3); % 初值重投影，黄色
    point1 = point2d_1(k, :);
    plot(point(1), point(2), 'g*'); 
    plot(point0(1), point0(2), 'b+'); 
%     if k == 3 || k == 7
%         plot(point1(1), point1(2), 'ro');
%     else
%         plot(point1(1), point1(2), 'bo');
%     end
end
 hold off
 
% 相机1
figure(2), clf, imshow('data\realdata001-1\cam2\fr00001.bmp');
hold on
for k = 1:1:8
    point = CamProj2 * e2h(Q(:, k)); point = point(1:2)/point(3); % 重投影
    point0 = CamProj2 * e2h(point3d(:, k)); point0 = point0(1:2)/point0(3); % 初值重投影，黄色
    point1 = point2d_2(k, :);
    hold on, plot(point(1), point(2), 'g*'); 
    plot(point0(1), point0(2), 'b+'); 
%     if k == 1 || k == 4
%         plot(point1(1), point1(2), 'ro');
%     else
%         plot(point1(1), point1(2), 'bo');
%     end
end
hold off
