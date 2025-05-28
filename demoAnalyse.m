% clear; clc; close all;
% load('result\realdata001-1\result_sm2.mat');

FPS = 200;  k = 12;
perFPS = 1000/FPS;
time = (tk.trackers(k).start : tk.trackers(k).end) * perFPS;
X = tk.trackers(k).states(1, :);
Y = tk.trackers(k).states(2, :);
Z = tk.trackers(k).states(3, :);

% time = (1:1:100) * perFPS; X = []; Y = []; Z = [];
% for k = 1:1:(length(tk.trackers)-1)
%     X = [X, tk.trackers(k).states(1, 1:end-1)];
%     Y = [Y, tk.trackers(k).states(2, 1:end-1)];
%     Z = [Z, tk.trackers(k).states(3, 1:end-1)];
% end
% k = k + 1;
% X = [X, tk.trackers(k).states(1, 1:end)];
% Y = [Y, tk.trackers(k).states(2, 1:end)];
% Z = [Z, tk.trackers(k).states(3, 1:end)];

%% 1 1只果蝇的位置变化
fig = figure(1);
ax = axes(fig);
hold on
plot(time, X, 'LineWidth', 1.5);
plot(time, Y, 'LineWidth', 1.5);
plot(time, Z, 'LineWidth', 1.5);
ylabel('Position (mm)');
xlabel('Time (ms)');
legend('$x$', '$y$', '$z$');
hold off
prettyPlot(fig, ax, 10, 5, 2);

%% 2 1只果蝇的速度变化
Vx = diff(X); Vy = diff(Y); Vz = diff(Z); 
fig = figure(2);
ax = axes(fig);
hold on
plot(time(2:end), Vx, 'LineWidth', 2);
plot(time(2:end), Vy, 'LineWidth', 2);
plot(time(2:end), Vz, 'LineWidth', 2);
ylabel('Velocity (mm/Frame)');
xlabel('Time (ms)');
legend('$v_x$', '$v_y$', '$v_z$');
hold off
prettyPlot(fig, ax, 10, 5, 2);
