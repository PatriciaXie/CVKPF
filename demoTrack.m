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
% % 1 ���ģ��
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
% % 2 ����ڷ�ʾ��ͼ
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

%% ��ʼ������ģ��
tic
nFrame = 50;
for idxCamera = 1:1:nCamera
    imgs = imgHandle.Read(idxCamera, 1:1:nFrame); 
    bgms(idxCamera, 1) = GaussianBackgroundModel(imgs);
end
timeAll = toc;
fprintf('�ɹ���������ģ��(%.1fs)\n', timeAll);

%% ����
% nParticles, sigma, dynamicModel, dt, useMultiThread
tk = FlyTracker(nParticle, sigma, dynamicModel, 1/FPS, false); % 1��ֹģ�͡�2����ģ�͡�3��ǰͳ��ģ�ͣ�% �Ƿ��ö��̣߳����̻߳���bug���ȱ���
idxFrame = 1;
for idxCamera = 1:1:nCamera
    img = imgHandle.Read(idxCamera, idxFrame);
	bgm = bgms(idxCamera);
    msm(idxCamera) = Measurement(img, bgm); % showMeasure(img, msm(idxCamera).blobs);
    bgms(idxCamera) = bgm.Update(img);
end
fprintf('Frame = %d����⵽%d,%d��blob(%.1fs)\n', idxFrame, length(msm(1).blobs), length(msm(2).blobs), toc - timeAll);
timeAll = toc;
%%
for idxFrame = 2:1:imgHandle.endFrame
    msm_ = msm; % ��t���⸳ֵ��t-1����
    for idxCamera = 1:1:nCamera
        img = imgHandle.Read(idxCamera, idxFrame);
        bgm = bgms(idxCamera);
        msm(idxCamera) = Measurement(img, bgm); % showMeasure(img, msm(idxCamera).blobs); 
        bgms(idxCamera) = bgm.Update(img);
    end
    fprintf('Frame = %d����⵽%d,%d��blob(%.1fs)��', idxFrame, length(msm(1).blobs), length(msm(2).blobs), toc - timeAll);
    timeAll = toc;
    if idxFrame == 2 % ��ʼ��Ŀ��
        newTrackers = tk.Detect(idxFrame, msm_, msm); % �¼���tracker��tk.newTracker��
        tk.trackers = newTrackers;
        if showInfo
            fprintf('�ɹ���ʼ��%d��tracker(%.1fs)\n', length(newTrackers), toc - timeAll);
        end
        timeAll = toc;
    else
        % ����Ŀ��
        trackers = tk.Track(idxFrame, msm); 
        if ~isempty(trackers)
            idxInActive = find([trackers.active]==0);
            tk.oldTrackers = [tk.oldTrackers, trackers(idxInActive)];
            trackers(idxInActive) = [];
            tk.trackers = trackers; tk.nTrackers = length(tk.trackers);
        end
        if showInfo
            fprintf('����%d��(%.1fs)��', tk.nTrackers, toc - timeAll);
        end
        timeAll = toc;
        % ��ͣ
        newTrackers = tk.Detect(idxFrame, msm_, msm);
        tk.trackers = [tk.trackers, newTrackers];
        if showInfo
            fprintf('׷��%d��(%.1fs)��', length(newTrackers), toc - timeAll);
            fprintf('���ƻ�Ծ%d����ʧЧ%d��\n', length(tk.trackers), length(tk.oldTrackers));
        end
        timeAll = toc;
    end
end
tk.trackers = [tk.oldTrackers, tk.trackers];
mkdir(sprintf('result/realdata001-2'));
% save(sprintf('result/realdata001-1/result_cm%d.mat', dynamicModel));
save(sprintf('result/realdata001-2/result_sm%d.mat', dynamicModel), 'tk', 'cube', 'stereoModel', 'imgHandle');