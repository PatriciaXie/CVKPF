% clear; clc; close all;
% load('result\realdata001-2\result_sm1.mat');
[Xb, Yb, Zb] = cube2line(cube);

% 2 相机摆放示意图
figure(3), clf;
hold on
plot3(Xb, Zb, Yb, 'k-');
for k = 1:1:length(tk.trackers)
    X = tk.trackers(k).states(1, :);
    Y = tk.trackers(k).states(2, :);
    Z = tk.trackers(k).states(3, :);
    plot3(X, Z, Y, 'b-');
end

% 之前的方法
% load('resultCore.mat');
% cube2.d3Vertices = box';
% [Xb2, Yb2, Zb2] = cube2line(cube2);
% plot3(Xb2, Zb2, Yb2, 'r-');
% for k = 1:1:height(tracklet_f)
%     X = tracklet_f.tXYZ{k}(:, 2); 
%     Y = tracklet_f.tXYZ{k}(:, 3);
%     Z = tracklet_f.tXYZ{k}(:, 4);
%     plot3(X, Z, Y, 'r-');
% end

xlabel('X'); ylabel('Z'); zlabel('Y');
 set(gca,'YDir','reverse');
axis equal
view(23, 29);
% axis([-300, 300, -500, 500, -50, 1100]);
grid on
hold off

%% 重投影
CamProj1 = stereoModel.cams(1).projection;
t = 1;
figure(1), clf, imshow(imgHandle.Read(1, t));
hold on
point3d = tk.trackers(1).states(1:3, t);
point = CamProj1 * e2h(point3d); point = point(1:2)/point(3);
plot(point(1), point(2), 'r+'); 
 hold off
 
% 相机1
CamProj2 = stereoModel.cams(2).projection;
figure(2), clf, imshow(imgHandle.Read(2, t));
hold on
point3d = tk.trackers(1).states(1:3, t);
point = CamProj2 * e2h(point3d); point = point(1:2)/point(3);
plot(point(1), point(2), 'r+'); 
hold off
