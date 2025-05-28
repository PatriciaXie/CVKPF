clear; clc; close all; warning off
addpath(genpath('src\'));
nCamera = 2; FPS = 500; dynamicModel = 1;
nParticle = 500; sigma = 3;
showInfo = true;

fprintf('Using model %d to track fly/flies...\n', dynamicModel);

datasetPath = 'data/realdata001-1'; 
cameraPath = sprintf('data/realdata001-1/stereoModel.mat'); 
cubePath = sprintf('data/realdata001-1/cube.mat'); 
modelPath = 'flyModel.mat';
imgHandle = ImageReader(nCamera, datasetPath, 'png');

global stereoModel cube morphModel
[x,y,z] = sphere(50);
sphereModel = 1.5*[x(:)';y(:)';z(:)'];
flyModel = load(modelPath);
flyModel.Xf = 1000*flyModel.Xf;
flyModel.Yf = 1000*flyModel.Yf;
flyModel.Zf = 1000*flyModel.Zf;
morphModel{1} = sphereModel;
morphModel{2} = flyModel;

load(cameraPath);
load(cubePath);

% % showScenario(stereoModel, cube);
% % 1 相机模型
% R1 = angle2dcm(0, 0, 0);
% T1 = stereoModel.cams(1).center;
% C1 = CentralCamera('focal', stereoModel.cams(1).intrinsic(1,1), 'pixel', 1e-6, ...
%     'resolution', [2048 2048], 'centre', [stereoModel.cams(1).intrinsic(1,3) stereoModel.cams(1).intrinsic(2,1)], ...
%     'distortion', [0, 0, 0, 0, 0],...
%     'pose', SE3(R1, T1), 'name', 'C1', 'color', [0, 0, 0]);
% 
% tmp = inv(stereoModel.cams(2).intrinsic) * stereoModel.cams(2).projection;
% R2 = tmp(:, 1:3);
% T2 = stereoModel.cams(2).center;
% C2 = CentralCamera('focal', stereoModel.cams(2).intrinsic(1,1), 'pixel', 1e-6, ...
%     'resolution', [2048 2048], 'centre', [stereoModel.cams(2).intrinsic(1,3) stereoModel.cams(2).intrinsic(2,1)], ...
%     'distortion', [0, 0, 0, 0, 0],...
%     'pose', SE3(R2, T2), 'name', 'C2', 'color', [0, 0, 0]);
% 
% [Xb, Yb, Zb] = cube2line(cube);
% % 2 相机摆放示意图
% figure(1);
% C1.plot_camera('scale', 100, 'color', 'r', 'label');
% C2.plot_camera('scale', 100, 'color', 'b', 'label');
% view(0, 0);
% hold on
% plot3(Xb, Yb, Zb, 'k');
% axis([-2000, 500, -1000, 1000, -500, 2000]);
% axis equal 
% grid on
% hold off

%% 初始化背景模型
tic
nFrame = 50;
for idxCamera = 1:1:nCamera
    imgs = imgHandle.Read(idxCamera, 1:1:nFrame); 
    bgms(idxCamera, 1) = GaussianBackgroundModel(imgs);
end
timeAll = toc;
fprintf('成功构建背景模型(%.1fs)\n', timeAll);

%% 跟踪
% nParticles, sigma, dynamicModel, dt, useMultiThread
tk = FlyTracker(nParticle, sigma, dynamicModel, 1/FPS, false); % 1静止模型、2匀速模型、3当前统计模型；% 是否用多线程，多线程还有bug，先别用
idxFrame = 1;
for idxCamera = 1:1:nCamera
    img = imgHandle.Read(idxCamera, idxFrame);
	bgm = bgms(idxCamera);
    msm(idxCamera) = Measurement(img, bgm); % showMeasure(img, msm(idxCamera).blobs);
    bgms(idxCamera) = bgm.Update(img);
end
fprintf('Frame = %d，检测到%d,%d个blob(%.1fs)\n', idxFrame, length(msm(1).blobs), length(msm(2).blobs), toc - timeAll);
timeAll = toc;
%%
for idxFrame = 2:1:imgHandle.endFrame
    msm_ = msm; % 把t量测赋值给t-1量测
    for idxCamera = 1:1:nCamera
        img = imgHandle.Read(idxCamera, idxFrame);
        bgm = bgms(idxCamera);
        msm(idxCamera) = Measurement(img, bgm); % showMeasure(img, msm(idxCamera).blobs); 
        bgms(idxCamera) = bgm.Update(img);
    end
    fprintf('Frame = %d，检测到%d,%d个blob(%.1fs)，', idxFrame, length(msm(1).blobs), length(msm(2).blobs), toc - timeAll);
    timeAll = toc;
    if idxFrame == 2 % 初始化目标
        newTrackers = tk.Detect(idxFrame, msm_, msm); % 新检测的tracker在tk.newTracker中
        tk.trackers = newTrackers;
        if showInfo
            fprintf('成功初始化%d个tracker(%.1fs)\n', length(newTrackers), toc - timeAll);
        end
        timeAll = toc;
    else
        % 跟踪目标
        trackers = tk.Track(idxFrame, msm); 
        if ~isempty(trackers)
            idxInActive = find([trackers.active]==0);
            tk.oldTrackers = [tk.oldTrackers, trackers(idxInActive)];
            trackers(idxInActive) = [];
            tk.trackers = trackers; tk.nTrackers = length(tk.trackers);
        end
        if showInfo
            fprintf('跟踪%d个(%.1fs)，', tk.nTrackers, toc - timeAll);
        end
        timeAll = toc;
        % 暂停
        newTrackers = tk.Detect(idxFrame, msm_, msm);
        tk.trackers = [tk.trackers, newTrackers];
        if showInfo
            fprintf('追加%d个(%.1fs)，', length(newTrackers), toc - timeAll);
            fprintf('共计活跃%d个，失效%d个\n', length(tk.trackers), length(tk.oldTrackers));
        end
        timeAll = toc;
    end
end
tk.trackers = [tk.oldTrackers, tk.trackers];
mkdir(sprintf('result/realdata001-2'));
% save(sprintf('result/realdata001-1/result_cm%d.mat', dynamicModel));
save(sprintf('result/realdata001-2/result_sm%d.mat', dynamicModel), 'tk', 'cube', 'stereoModel', 'imgHandle');