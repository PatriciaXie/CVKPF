classdef FlyTracker
    %Tracker 轨迹片段类
    %   含有轨迹的定义、检测、跟踪
    
    properties
        nTrackers
        
        oldTrackers
        trackers
        nParticles
        epipolarError
        mergeThreshold
        useMultiThread
        dt
        dynamicModel
        sigma
    end
    
    methods
        function tk = FlyTracker(nParticles, sigma, dynamicModel, dt, useMultiThread)
            %Tracker 构造此类的实例
            %   定义Tracker参数
            tk.nTrackers = 0;
            tk.nParticles = nParticles;
            tk.epipolarError = [5, 3];
            tk.mergeThreshold = 1.3;
            tk.useMultiThread = useMultiThread;
            tk.oldTrackers = [];
            tk.dt = dt; 
            tk.dynamicModel = dynamicModel;
            tk.sigma = sigma;
        end
        
        function newTrackers = Detect(tk, t, msm_, msm)
            %Detect 检测未被解析的blobs，重建新目标
            %   time: 时间（帧号）
            %   msm_: t-1时刻的量测
            %   msm  : t时刻的量测
            if tk.nTrackers > 0
                timeEnd = cat(2, tk.trackers.end);
                idxTracked = find(timeEnd == t);
            else
                idxTracked = [];
            end
            nCamera = length(msm);
            temporalPairs = cell(nCamera, 1);
            for idxCamera = 1:1:nCamera
                temporalPairs{idxCamera} = TemporalAssociation(msm_, msm, idxCamera);
                if isempty(temporalPairs{idxCamera})
                    continue;
                end
                % 对时间关联帧进行筛选，去掉已经在跟踪阶段关联上的对
                if ~isempty(idxTracked)
                    temporPairsTracked = zeros(size(temporalPairs{idxCamera}, 1), 1);
                    ids_ = []; ids = [];
                    for k = 1:1:length(idxTracked)
                        idxes_ = tk.trackers(idxTracked(k)).info.idxBlob{idxCamera}(end-1);
                        idxes = tk.trackers(idxTracked(k)).info.idxBlob{idxCamera}(end);
                        if ~isempty(idxes_)
                            ids_ = [ids_; find(idxes_ == temporalPairs{idxCamera}(:, 1))];
                        end
                        if ~isempty(idxes)
                            ids = [ids; find(idxes == temporalPairs{idxCamera}(:, 2))];
                        end
                    end
                    id = intersect(ids_, ids);
                    temporalPairs{idxCamera}(id, :) = [];
                end
            end
            [targets, nTargets] = SpatialAssociation(msm_, msm, temporalPairs, tk.epipolarError, tk.mergeThreshold);
            for i = 1:1:nTargets
                target = targets(i);
                newTracker = CreateTracker(t, target, tk.nParticles, tk.dt, tk.sigma, tk.dynamicModel);
                tmpTrackers(i) = newTracker;
            end
            if nTargets > 0
                newTrackers = tmpTrackers;
            else
                newTrackers = [];
            end
        end
       
        function trackers = Track(tk, time, msm)
            global stereoModel cube morphModel
            sM = stereoModel; cb = cube; mM = morphModel;
            %Track 跟踪器匹配新的量测
            %   time: 时间（帧号）
            %   msm  : t时刻的量测
            if isempty(tk.trackers)
                trackers = [];
                return;
            end
            idxActive = find([tk.trackers.active]>0);
            if tk.useMultiThread
                % 将idxActive分成nWorker个Task
                nWorker = 4;
                if length(idxActive) < nWorker
                    nWorker = length(idxActive);
                end
                taskIdx = cell(nWorker, 1);
                for i = 1:1:nWorker
                    taskIdx{i} = idxActive(i:nWorker:end);
                end
                c = parcluster;
                job = createJob(c);
                
                for i = 1:1:nWorker
%                     createTask(job, @singlyTrack, 1, {tk.trackers(taskIdx{i}(1)), tk.dynamicModel, time, msm, stereoModel, cube, morphModel, tk.sigma});
                    createTask(job, @jobTrack, 1, {msm, tk.trackers(taskIdx{i}), tk.dynamicModel, time, sM, cb, mM, tk.sigma});
                end
                submit(job);
                wait(job);
                out = fetchOutputs(job);
                for i = 1:1:nWorker
                    for j = 1:1:length(taskIdx{i})
                        idx = taskIdx{i}(j);
                        newTracker = out{i}(j);
                        trackers(idx) = newTracker;
                    end
                end
                delete(job);
            else
                for idxTracker = idxActive
                    tracker_ = tk.trackers(idxTracker);
                    tracker = singlyTrack(tracker_, msm, tk.dynamicModel, time, stereoModel, cube, morphModel, tk.sigma);
                    trackers(idxTracker) = tracker;
                end
            end

        end
        
    end
end


function temporalPairs = TemporalAssociation(msm_, msm, idxCamera)
    temporalPairs = [];
    imageH = msm(idxCamera).imageSize(2); imageW = msm(idxCamera).imageSize(1);
    radius = 50;
    [hGrid, wGrid] = meshgrid(-radius:radius,-radius:radius);
    circleMask = find(hGrid.^2 + wGrid.^2 <= radius^2);
    circleIdxH = hGrid(circleMask);
    circleIdxW = wGrid(circleMask);
    
    if isempty(msm_(idxCamera).blobs) || isempty(msm(idxCamera).blobs)
        return
    end
    blobIdxRange_ = cat(2, msm_(idxCamera).blobs.idxRange); % msm_所有的blobs_像素范围，用于根据blob快速找到radius内的blob_
    for idxBlob = 1:1:length(msm(idxCamera).blobs)
        blob = msm(idxCamera).blobs(idxBlob);
        blobIdx1 = msm(idxCamera).blobs(idxBlob).idx1;
        blobImg1 = msm(idxCamera).image(blobIdx1);
        ellipse = msm(idxCamera).blobs(idxBlob).ellipse;
        
        circleIdx2_h = circleIdxH + round(blob.center(2));
        circleIdx2_w = circleIdxW + round(blob.center(1));
        h1 = min(circleIdx2_h); h2 = max(circleIdx2_h);
        w1 = min(circleIdx2_w); w2 = max(circleIdx2_w);
        if h1 >= 1 && h2 <= imageH && w1 >= 1 && w2 <= imageW
            circleIdx1 = sub2ind(msm(idxCamera).imageSize, circleIdx2_h, circleIdx2_w);
            possibleIdxBlob_ = find(blobIdxRange_(1,:) >= min(circleIdx1) & blobIdxRange_(2,:) <= max(circleIdx1));
            associateIdx = []; similarity = []; % 存储下一帧与该blob关联的blob_编号，即关联的相似度
            for idxBlob_ = possibleIdxBlob_
                blob_ = msm_(idxCamera).blobs(idxBlob_);
                overlapIdx = blob_.idx1(ismember(blob_.idx1, circleIdx1));
                if (~isempty(overlapIdx))
%                     fprintf('帧间关联：Camera(%d).blob(%d) @ (%.1f, %.1f)与blob(%d) @ (%.1f, %.1f)有重合区域\n', idxCamera, idxBlob, blob.center(1), blob.center(2), idxBlob_, blob_.center(1), blob_.center(2));
                    associateIdx = [associateIdx, idxBlob_];
                    blobIdx1_ = msm_(idxCamera).blobs(idxBlob_).idx1;
                    blobImg1_ = msm_(idxCamera).image(blobIdx1_);
                    nccValue = GetNcc(blobImg1, blobImg1_);
                    ellipse_ = msm_(idxCamera).blobs(idxBlob_).ellipse;
                    try
                        ellipseDistance = GetEllipseDistance(ellipse, ellipse_); % 计算两个椭圆的相似度
                    catch
                        a = 1;
                    end
                    similarity = [similarity, nccValue * (1/ellipseDistance)]; % 计算blob的外观相似度
                end
            end
            if ~isempty(associateIdx)
                [~, maxIdx] = max(similarity);
                temporalPairs = [temporalPairs; associateIdx(maxIdx), idxBlob];
            end
        else
%             fprintf('帧间关联：Camera(%d).blob(%d) @ (%.1f, %.1f)靠近边缘，不追溯其上一时刻的关联检测\n', idxCamera, idxBlob, blob.center(1), blob.center(2));
        end
    end

end

function [targets, nTargets] = SpatialAssociation(msm_, msm, temporalPairs, epipolarError, mergeThreshold)
%spatialAssociation 同一时刻测量的相机关联，极限约束
%   输入
%       msm_，t-1时刻的量测
%       msm，t时刻的量测
%       temporalPairs，各个相机的帧间匹配
%       epipolarError，极线误差
%       mergeThreshold，跟椭圆的pi*a*b/像素数有关的，不懂

%	点P(X,Y)到直线Ax+By+C=0的距离为 |AX+BY+C|  除以 根号下(A^2+B^2)

    global stereoModel cube
    
    nTargets = 0; targets = [];
    idxCamera1 = 1; idxCamera2 = 2; % 只考虑相机1、2之间第t帧的极线关系
    for pairIdx = 1:1:size(temporalPairs{idxCamera1}, 1)
        idx1_ = temporalPairs{idxCamera1}(pairIdx, 1);
        idx1 = temporalPairs{idxCamera1}(pairIdx, 2);
        point1 = msm(idxCamera1).blobs(idx1).center;
        desireCoefs = stereoModel.fundamentals{1} * [reshape(point1, 2, 1); 1];

        associateIdx = []; associateIdx_ = [];
        for i=1:size(temporalPairs{idxCamera2}, 1)
            idx2_ = temporalPairs{idxCamera2}(i, 1);
            idx2 = temporalPairs{idxCamera2}(i, 2);
            point2 = msm(idxCamera2).blobs(idx2).center;
            distance = ( [reshape(point2,1,2) 1] * desireCoefs ) / sqrt( desireCoefs(1)^2 + desireCoefs(2)^2 );
            if ( abs(distance) <= epipolarError(2) ) 
                associateIdx = [associateIdx; idx2]; % t时刻与相机1的point1关联的相机2点的量测编号
                associateIdx_ = [associateIdx_; idx2_]; % t-1时刻与相机1的point1关联的相机2点的t-1关联量测编号
            end
        end
        
        % 相机1上t
        blob1 = msm(idxCamera1).blobs(idx1);
%         gamma1 = blob1.ellipse.a / blob1.ellipse.b;
%         factor1 = blob1.ellipse.factor;
%         if factor1 > mergeThreshold
%             gamma1 = 1 / factor1;
%         end
        % 相机1上t-1
        blob1_ = msm_(idxCamera1).blobs(idx1_);
%         gamma1_ = blob1_.ellipse.a / blob1_.ellipse.b;
%         factor1_ = blob1_.ellipse.factor;
%         if factor1_ > mergeThreshold
%             gamma1_ = 1 / factor1_;
%         end
        
        for i = 1:1:length(associateIdx)
            % 相机2上t
            idx2 = associateIdx(i);
            blob2 = msm(idxCamera2).blobs(idx2);
%             gamma2 = blob2.ellipse.a / blob2.ellipse.b;
%             factor2 = blob2.ellipse.factor;
%             if factor2 > mergeThreshold
%                 gamma2 = 1 / factor2;
%             end
            % 相机2上t-1
            idx2_ = associateIdx_(i);
            blob2_ = msm_(idxCamera2).blobs(idx2_);
%             gamma2_ = blob2_.ellipse.a / blob2_.ellipse.b;
%             factor2_ = blob2_.ellipse.factor;
%             if factor2_ > mergeThreshold
%                 gamma2_ = 1 / factor2_;
%             end
            
            % 计算t时刻的三维坐标
            CamProj1 = stereoModel.cams(1).projection;
            CamProj2 = stereoModel.cams(2).projection;
            point1 = blob1.center; point2 = blob2.center;
            point3d = ReconstructPosition(CamProj1, CamProj2, point1, point2);
            inFlag = IsInCube(point3d, cube);
%             figure(1),hold on, plot3(point3d(1), point3d(2), point3d(3), 'yo'); hold off
            if inFlag              
                % 重建方向，gamma越大的越好
%                 figure(1),hold on, plot3(point3d(1), point3d(2), point3d(3), 'ro'); hold off
                ellipse1 = blob1.ellipse; ellipse2 = blob2.ellipse;
                orientation = ReconstructOrientation(CamProj1, CamProj2, ellipse1, ellipse2);
            else
%                 fprintf('空间关联2：t时刻的三维点在盒子外，不计算t-1时刻三维坐标\n');
                continue; % 如果t时刻计算的三维坐标在盒子外，则t-1、t这个关联对也无效了。否则，去计算t-1时刻的三维坐标。
            end
            % 验证相机1上t-1和相机2上t-1是不是满足极线约束
            point1_ = blob1_.center; point2_ = blob2_.center;
            desireCoefs = stereoModel.fundamentals{1} * [reshape(point1_, 2, 1); 1];
            distance = abs( ( [reshape(point2_,1,2) 1] * desireCoefs ) / sqrt( desireCoefs(1)^2 + desireCoefs(2)^2 ));
            if ( distance <= epipolarError(1) )
                % 计算三维坐标
                point3d_ = ReconstructPosition(CamProj1, CamProj2, point1_, point2_);
                % 验证三维坐标，不能在壁上或盒子外面
                inFlag_ = IsInCube(point3d_, cube);
                if inFlag_
                    % 重建方向，gamma越大的越好
                    ellipse1_ = blob1_.ellipse; ellipse2_ = blob2_.ellipse;
                    orientation_ = ReconstructOrientation(CamProj1, CamProj2, ellipse1_, ellipse2_);
                    nTargets = nTargets + 1;
%                     fprintf('空间关联4：确认目标 %d\n', nTargets);
                    targets(nTargets).point3d_ = point3d_;
                    targets(nTargets).point3d = point3d;
                    targets(nTargets).orientation_ = orientation_;
                    targets(nTargets).orientation = orientation;
                    targets(nTargets).appearance.ellipse{1} = blob1.ellipse;
                    targets(nTargets).appearance.ellipse{2} = blob2.ellipse;
                    targets(nTargets).appearance.intensity{1} = msm(idxCamera1).image(blob1.idx1);
                    targets(nTargets).appearance.intensity{2} = msm(idxCamera2).image(blob2.idx1);
                    targets(nTargets).info.idxBlob{1} = [idx1_, idx1];
                    targets(nTargets).info.idxBlob{2} = [idx2_, idx2];
                else
%                     fprintf('空间关联3：t时刻三维点在盒子内，但t-1时刻的三维点在盒子外，该关联对无效\n');
                    continue;
                end
            else
%                 fprintf('空间关联1：t-1时刻，Camera(1).blob(%d) @(%.1f, %.1f)和Camera(2).blob(%d) @(%.1f, %.1f)不满足极线约束\n',...
%                     idx1_, blob1_.center(1), blob1_.center(2), idx2_, blob2_.center(1), blob2_.center(2));
            end
        end
    end
end

function tracker = CreateTracker(t, target, nParticles, dt, sigma, dynamicModel)
    tracker.start = t-1; tracker.end = t; tracker.active = 1; tracker.missing = 0;
    [theta_, phi_, ~] = cart2sph(target.orientation_(1), target.orientation_(3), -target.orientation_(2));
    state_ = [target.point3d_; theta_; phi_; 0; 0; 0];
    
    [theta, phi, ~] = cart2sph(target.orientation(1), target.orientation(3), -target.orientation(2));
    state = [target.point3d; theta; phi; target.point3d_];
    tracker.states = [state_, state];
    
    newState = state;
    switch dynamicModel
        case 1
            % begin{无模型}
            tracker.particles = repmat(newState, 1, nParticles); % 无当前统计模型
            tracker.prediction = newState(1:3);
            % end{无模型}
            
        case 2
            % begin{FLE模型}
            matTransition = [
                2 0 0 0 0 -1  0  0;
                0 2 0 0 0  0 -1  0;
                0 0 2 0 0  0  0 -1;
                0 0 0 1 0  0  0  0;
                0 0 0 0 1  0  0  0;
                1 0 0 0 0  0  0  0;
                0 1 0 0 0  0  0  0;
                0 0 1 0 0  0  0  0]; %FLE 状态转移模型
            tmpState = matTransition * newState;
            tracker.particles = repmat(tmpState, 1, nParticles); % 无当前统计模型
            tracker.prediction = tmpState(1:3);
            % end{FLE模型}
            
        case 3
            % begin{当前统计模型}
            % 建立当前统计模型并更新
            alpha = 0.01; aMax = 5; R = 1e-9; T = dt;
            tracker.csm.state{1} = [tracker.states(1, 1); 0; 0];
            tracker.csm.state{2} = [tracker.states(2, 1); 0; 0];
            tracker.csm.state{3} = [tracker.states(3, 1); 0; 0];
            tracker.csm.P{1} = 1e-6 * eye(3);
            tracker.csm.P{2} = 1e-6 * eye(3);
            tracker.csm.P{3} = 1e-6 * eye(3);
            tracker.csm.filter{1} = FilterCSMKCF(alpha, aMax, R, T); % x通道当前统计模型
            tracker.csm.filter{2} = FilterCSMKCF(alpha, aMax, R, T); % y通道当前统计模型
            tracker.csm.filter{3} = FilterCSMKCF(alpha, aMax, R, T); % z通道当前统计模型
            tracker.csm.aAvg = [0, 0, 0];
            for i = 1:3
                filterObj = tracker.csm.filter{i};
                X = tracker.csm.state{i};
                P = tracker.csm.P{i};
                aAvg = tracker.csm.aAvg(i);
                Z = tracker.states(i, 2);
                [X, P, ~, ~, aAvg] = filterObj.Filter(X, P, aAvg, Z); % Xm、Pm没有用到观测
                tracker.csm.aAvg(i) = aAvg;
                tracker.csm.state{i} = X;
                tracker.csm.P{i} = P;
            end
            
            % 进行预测
            csmState = [0; 0; 0; newState(4:8)];
            for i = 1:3
                filterObj = tracker.csm.filter{i};
                X = tracker.csm.state{i};
                P = tracker.csm.P{i};
                aAvg = tracker.csm.aAvg(i);
                Z = newState(i, end);
                [~, ~, Xm, ~, ~] = filterObj.Filter(X, P, aAvg, Z); % Xm、Pm没有用到观测
                csmState(i) = Xm(1);
            end
            tracker.prediction = csmState(1:3);
            tracker.particles = repmat(csmState, 1, nParticles); 
            % end{当前统计模型}
        otherwise
                fprintf('没有运动模型%d\n', dynamicModel);
    end
    
    tracker.particles(1:3, :) = tracker.particles(1:3, :) + sigma*randn(3, nParticles); % gaussian
    tracker.weights = repmat(1/nParticles, 1, nParticles); % 初始时，每个粒子的权重相同
    tracker.appearance = target.appearance;
    tracker.info = target.info;

end

function newTracker = singlyTrack(tracker, msm, dynamicModel, t, stereoModel, cube, morphModel, sigma)
%     global morphModel 
    nParticles = length(tracker.weights);
%     sigma = 0.003; pNoise = sigma*randn(3, nParticles); % gaussian
%     sigma = 0.001; pNoise = sigma*randn(3, nParticles);
    % 测试一下500个粒子的高斯分布
%     trajectory = tracker.states(1:3, max(1, end-10):end); % 最近10帧
    tmpParticles = zeros(size(tracker.particles));
    tmpWeights = zeros(size(tracker.weights));
    nMiss = 0;
    
% %  debug begin    
%     imgSize = msm(1).imageSize;
%     img{1} = zeros(imgSize);
%     img{2} = zeros(imgSize);
%     center{1} = [0; 0]; center{2} = [0; 0];
%     for idxCamera = 1:1:2
%         for n = 1:1:nParticles
%             point3d = tracker.particles(1:3, n);
%             sphere = morphModel{1} + point3d;
%             pixels2 =stereoModel.cams(idxCamera).projection * [sphere; ones(1,size(sphere,2))];
%             pixels2(1:2, :) = pixels2(1:2, :) ./ repmat(pixels2(3, :), 2, 1);
%             pixels2 = floor(pixels2(1:2,:));
%             pixels2 = unique(pixels2','rows')';
%             pixels1 = sub2ind(imgSize, pixels2(2,:), pixels2(1, :));
%             img1 = showRePorjection(pixels1, msm(idxCamera).blobs(1), imgSize);
%             img1 = img{idxCamera} + img1/nParticles;
%             img1(img1>1) = 1;
%             img{idxCamera} = img1;
%             center{idxCamera} = center{idxCamera} + mean(pixels2, 2)/nParticles;
%         end
%         xs = round(center{idxCamera}(2)) - 100 : round(center{idxCamera}(2)) + 100;
%         ys = round(center{idxCamera}(1)) - 100 : round(center{idxCamera}(1)) + 100;
%         figure(idxCamera+4)
%         imshow(img{idxCamera}(xs, ys, :));
%         set(gcf,'unit','centimeters','position',[1,1,11,11]);
%         set(gca,'unit','centimeters','position',[1,1,9,9]);
%         hold on
% %         plot(101, 101, 'y*');
%         hold off
%     end
% % debug end   

    % 粒子更新
    for n = 1:1:nParticles
        labelParticles(n).associated = 0;
        labelParticles(n).pal = -1; % pal 是什么？
        labelParticles(n).pol = -1;
        tmpParticles(:,n) =  tracker.particles(:,n);% + [pNoise(:,n); zeros(5,1)];
        point3d = tmpParticles(1:3,n);
        
        % 生成的粒子要在盒子内
        inFlag = IsInCube(point3d, cube);
        if ~inFlag
            nMiss = nMiss + 1;
            tmpWeights(n) = 0;
            continue;
        end
        
        % 生成的粒子在各视图上的投影要与现有的blobs有重叠（可能与0、1、2、...个blob重叠）重叠的，0则miss
        okFlag = 1;
        state = [point3d; tracker.states(4:8, end)];
        
        for idxCamera = 1:1:length(msm)
            [idxes{idxCamera}, ratios{idxCamera}] = GetOverlap(stereoModel.cams(idxCamera).projection, stereoModel.camResolution, state, msm(idxCamera).blobs, 'SphereShape', morphModel);
            if isempty(idxes{idxCamera})
                nMiss = nMiss + 1;
                tmpWeights(n) = 0;
                okFlag = 0;
                break;
            end
        end
        
        if okFlag % 该粒子在各视图上的投影均与现有的blobs重叠
            labelParticles(n).associated = 1; 
            for  idxCamera = 1:1:length(msm)
                labelParticles(n).msm{idxCamera} = [];
                for i = 1:1:length(idxes{idxCamera})
                    idx = idxes{idxCamera}(i);
                    blob = msm(idxCamera).blobs(idx);
                    blobImage = msm(idxCamera).image(blob.idx1);
                    blobEllipse = blob.ellipse;
                    labelParticles(n).msm{idxCamera} = [labelParticles(n).msm{idxCamera}, idx];
                    % 像素相似度
                    pal = ratios{idxCamera}(i) * GetNcc(tracker.appearance.intensity{idxCamera}, blobImage); 
                    % 姿态相似度
                    pol = GetEllipseDistance(tracker.appearance.ellipse{idxCamera}, blobEllipse);
% %                     debug begin
%                     imgSize = msm(1).imageSize;
%                     point3d = tmpParticles(1:3, n);
%                     sphere = morphModel{1} + point3d;
%                     pixels2 =stereoModel.cams(idxCamera).projection * [sphere; ones(1,size(sphere,2))];
%                     pixels2(1:2, :) = pixels2(1:2, :) ./ repmat(pixels2(3, :), 2, 1);
%                     pixels2 = floor(pixels2(1:2,:));
%                     pixels2 = unique(pixels2','rows')';
%                     pixels1 = sub2ind(imgSize, pixels2(2,:), pixels2(1, :));
%                     img1 = showRePorjection(pixels1, blob, imgSize);
%                     figure(1), imshow(img1);
% %                     debug end
            
                    labelParticles(n).weights{idxCamera}(i, :) = [pal, pol];
                end
            end
        end
    end
    
    if nMiss > nParticles - 2
        tracker.active = 0;
        newTracker = tracker;
        return;
    end
    
    % 对有associate的粒子进行分析  
    idxAssociated = find([labelParticles.associated]>0);
    if isempty(idxAssociated)
        nMiss = nParticles;
        tracker.missing = tracker.missing + 1;
        missThreshold = 1;
        if tracker.missing >  missThreshold % 如果miss
            tracker.active = 0;
            newTracker = tracker;
            return;
        end
    else
        for idxCamera = 1:1:2
            idxBlobAssociated{idxCamera} = [];
            for i = idxAssociated
                idxBlobAssociated{idxCamera} = [idxBlobAssociated{idxCamera}, labelParticles(i).msm{idxCamera}];
            end
            idxBlobAssociated{idxCamera} = unique(idxBlobAssociated{idxCamera});
        end
    end
    
    % 将idxAssociated的粒子与可能的msm.blobs关联的pal, pol计算出来
    for i = 1:1:length(idxAssociated)
        aParticle = labelParticles(idxAssociated(i));
        for idxCamera = 1:1:2
            aParticle.msm{idxCamera} = unique(aParticle.msm{idxCamera});
            for r = 1:1:length(idxBlobAssociated{idxCamera})
                idx =  find(aParticle.msm{idxCamera} == idxBlobAssociated{idxCamera}(r), 1); % 粒子关联的量测
                if isempty(idx)
                    weights{idxCamera}.pal(i, r) = 0;
                    weights{idxCamera}.pol(i, r) = 100;
                else
                    weights{idxCamera}.pal(i, r) = aParticle.weights{idxCamera}(idx, 1);
                    weights{idxCamera}.pol(i, r) = aParticle.weights{idxCamera}(idx, 2);
                end
            end
        end
    end
    
    % 按照像素相似度对粒子对应的blobs进行排序
    for idxCamera = 1:1:2
        [~, idxOrder{idxCamera}] = sort(sum(weights{idxCamera}.pal), 'descend');
    end
    
    imgSize = msm(1).imageSize;
    img{1} = zeros(imgSize);
    img{2} = zeros(imgSize);
    center{1} = [0; 0];
    center{2} = [0; 0];
    max_pal = [0, 0];
    % 更新粒子权重
    for i = 1:1:length(idxAssociated)
        aParticle = labelParticles(idxAssociated(i));
        for idxCamera = 1:1:2
% %                 debug begin
%                 point3d = tmpParticles(1:3, idxAssociated(i));
%                 sphere = morphModel{1} + point3d;
%                 pixels2 =stereoModel.cams(idxCamera).projection * [sphere; ones(1,size(sphere,2))];
%                 pixels2(1:2, :) = pixels2(1:2, :) ./ repmat(pixels2(3, :), 2, 1);
%                 pixels2 = floor(pixels2(1:2,:));
%                 pixels2 = unique(pixels2','rows')';
%                 pixels1 = sub2ind(imgSize, pixels2(2,:), pixels2(1, :));
%                 img1 = showRePorjection(pixels1, msm(idxCamera).blobs(idxOrder{idxCamera}(1)), imgSize);
%                 img1 = img{idxCamera} + img1/length(idxAssociated);
%                 img{idxCamera} = img1;
%                 center{idxCamera} = center{idxCamera} + mean(pixels2, 2)/length(idxAssociated);
% %                 if max_pal(idxCamera) < aParticle.weights{idxCamera}(idx, 1)
% %                     max_pal(idxCamera) =  aParticle.weights{idxCamera}(idx, 1);
% %                     center{idxCamera} = mean(pixels2, 2);
% %                     tmpWeights(tmpWeights>0) = 0;
% %                     tmpWeights(idxAssociated(i)) = 1;
% %                 end
% %                 debug end
            
            % 最相似的blob为
            mostBlob = idxBlobAssociated{idxCamera}(idxOrder{idxCamera}(1) );
            idx = find(aParticle.msm{idxCamera} == mostBlob);
            if isempty(idx)
                idx = ismember(aParticle.msm{idxCamera}, idxBlobAssociated{idxCamera}(idxOrder{idxCamera}(2:end)));
                idx = find(idx);
                labelParticles(idxAssociated(i)).pal = aParticle.weights{idxCamera}(idx(1), 1);
                labelParticles(idxAssociated(i)).pol = aParticle.weights{idxCamera}(idx(1), 2);
            else
                labelParticles(idxAssociated(i)).pal = 2*aParticle.weights{idxCamera}(idx,1); % 2和0.5应该是对像素和姿态进行加权
                labelParticles(idxAssociated(i)).pol = 0.5*aParticle.weights{idxCamera}(idx,2);
            end
        end
    end
%     % debug begin
%     xs = round(center{1}(2)) - 50 : round(center{1}(2)) + 50;
%     ys = round(center{1}(1)) - 50 : round(center{1}(1)) + 50;
%     figure(1)
%     imshow(img{1}(xs, ys, :));
%     set(gcf,'unit','centimeters','position',[1,1,11,11]);
%     set(gca,'unit','centimeters','position',[1,1,9,9]);
%     hold on
%     plot(51, 51, 'y*');
%     hold off
%     
%     xs = round(center{2}(2)) - 50 : round(center{2}(2)) + 50;
%     ys = round(center{2}(1)) - 50 : round(center{2}(1)) + 50;
%     figure(2)
%     imshow(img{2}(xs, ys, :)); 
%     set(gcf,'unit','centimeters','position',[15,1,11,11]);
%     set(gca,'unit','centimeters','position',[1,1,9,9]);
%     hold on
%     plot(51, 51, 'y*');
%     hold off
%     % debug end
    
    % 不加权
%     for i = 1:1:length(idxAssociated)
%         tmpWeights(idxAssociated(i)) = 1;
%     end
    % 加权
	pals = cat(2, labelParticles.pal); pols = cat(2, labelParticles.pol);
    w1 = pals(idxAssociated);
    if max(pols(idxAssociated)) ~= min(pols(idxAssociated))
        w2 = -(pols(idxAssociated)-max(pols(idxAssociated))) / (max(pols(idxAssociated))-min(pols(idxAssociated)));
    else
        w2 = pols(idxAssociated) / max(pols(idxAssociated));
    end
    tmpWeights(idxAssociated) = w1 .* exp(3*w2);
    
    tracker.missing = 0;
    tmpWeights = tmpWeights / sum(tmpWeights);
    newState = tmpParticles * tmpWeights'; % 综合的粒子可能与原粒子不同，因此要重新验证在不在盒子内，以及是否能与量测关联
    newPoint3d = newState(1:3);
    newState = [newPoint3d; tracker.states(4:8, end)]; % state只更新状态，不管姿态
    
    inFlag = IsInCube(newPoint3d, cube); % 验证求出的粒子是否在盒子内
    if ~inFlag
        tracker.active = 0;
        newTracker = tracker;
        return;
    end
    % 验证求出的粒子是否关联到量测，这次不用球，而用椭球模型
    clear idxes ratios
    for idxCamera = 1:1:2
        [idxes, ratios] = GetOverlap(stereoModel.cams(idxCamera).projection, stereoModel.camResolution, newState, msm(idxCamera).blobs, 'EllipsoidShape', morphModel);
        if isempty(idxes)
            tracker.active = 0;
            newTracker = tracker;
            return;
        end
        [~, idx] = max(ratios);
        idxOverlapBlob(idxCamera) = idxes(idx);
    end
    
    % idxOverlapBlob就是tracker和msm的关联
    CamProj1 = stereoModel.cams(1).projection;
    CamProj2 = stereoModel.cams(2).projection;
    blob1 = msm(1).blobs(idxOverlapBlob(1));
    blob2 = msm(2).blobs(idxOverlapBlob(2));
    point1 = blob1.center; point2 = blob2.center;
    point3d = ReconstructPosition(CamProj1, CamProj2, point1, point2);
%     newState(1:3) = newState(1:3); % 结果来自粒子滤波
% %     newState(1:3) = point3d; % 结果完全取决于检测
    
    % 重建姿态
    ellipse1 = blob1.ellipse; ellipse2 = blob2.ellipse;
    orientation = ReconstructOrientation(CamProj1, CamProj2, ellipse1, ellipse2);
    [theta, phi, ~] = cart2sph(orientation(1), orientation(3), -orientation(2));
    newState(4:8) = [theta; phi; tracker.states(1:3, end)];
    
%     % debug begin
%     for idxCamera = 1:1:2
%         ellipsoid = GetShape(morphModel{2}, newState);
%         pixels2 =stereoModel.cams(idxCamera).projection * [ellipsoid; ones(1,size(ellipsoid,2))];
%         pixels2(1:2, :) = pixels2(1:2, :) ./ repmat(pixels2(3, :), 2, 1);
%         pixels2 = floor(pixels2(1:2,:));
%         pixels2 = unique(pixels2','rows')';
%         pixels1 = sub2ind(imgSize, pixels2(2,:), pixels2(1, :));
%         img{idxCamera} = showRePorjection(pixels1, msm(idxCamera).blobs(idxOrder{idxCamera}(1)), imgSize);
%         xs = round(blob1.center(2)) - 50 : round(blob1.center(2)) + 50;
%         ys = round(blob1.center(1)) - 50 : round(blob1.center(1)) + 50;
%         figure(3), imshow(img{1}(xs, ys, :)); 
%         set(gcf,'unit','centimeters','position',[1,1,11,11]);
%         set(gca,'unit','centimeters','position',[1,1,9,9]);
%         xs = round(blob2.center(2)) - 50 : round(blob2.center(2)) + 50;
%         ys = round(blob2.center(1)) - 50 : round(blob2.center(1)) + 50;
%         figure(4), imshow(img{2}(xs, ys, :));
%         set(gcf,'unit','centimeters','position',[1,1,11,11]);
%         set(gca,'unit','centimeters','position',[1,1,9,9]);
%     end
%     % debug end
    
    tracker.end = t; tracker.missing = 0; tracker.active = 1;
    
    % 更新粒子
    switch dynamicModel
        case 1
            % begin{无模型}
            tracker.particles = repmat(newState, 1, nParticles); % 无当前统计模型
            tracker.prediction = [tracker.prediction, newState(1:3)];
            % end{无模型}
            tracker.states = [tracker.states, newState];
        case 2
            % begin{FLE模型}
            matTransition = [
                2 0 0 0 0 -1  0  0;
                0 2 0 0 0  0 -1  0;
                0 0 2 0 0  0  0 -1;
                0 0 0 1 0  0  0  0;
                0 0 0 0 1  0  0  0;
                1 0 0 0 0  0  0  0;
                0 1 0 0 0  0  0  0;
                0 0 1 0 0  0  0  0]; %FLE 状态转移模型
            tmpState = matTransition * newState;
            tracker.particles = repmat(tmpState, 1, nParticles); % 无当前统计模型
            tracker.prediction = [tracker.prediction, tmpState(1:3)];
            % end{FLE模型}
            tracker.states = [tracker.states, newState];
        case 3
            % begin{当前统计模型}
            % 更新
            if (tracker.end -tracker.start) > 10
                filterObj.R = 1e-3;
            end
            for i = 1:3
                filterObj = tracker.csm.filter{i};
                X = tracker.csm.state{i};
                P = tracker.csm.P{i};
                aAvg = tracker.csm.aAvg(i);
                Z = newState(i, end);
                [X, P, ~, ~, aAvg] = filterObj.Filter(X, P, aAvg, Z); % Xm、Pm没有用到观测
                tracker.csm.state{i} = X;
                tracker.csm.P{i} = P;
                tracker.csm.aAvg(i) = aAvg;
            end
            %预测
            csmState = [0; 0; 0; newState(4:8)];
            for i = 1:3
                filterObj = tracker.csm.filter{i};
                X = tracker.csm.state{i};
                P = tracker.csm.P{i};
                aAvg = tracker.csm.aAvg(i);
                Z = newState(i, end);
                [~, ~, Xm, ~, ~] = filterObj.Filter(X, P, aAvg, Z); % Xm、Pm没有用到观测
                csmState(i) = Xm(1);
            end
            tracker.prediction = [tracker.prediction, csmState(1:3)];
            tracker.particles = repmat(csmState, 1, nParticles); 
            % end{当前统计模型}
            newState(1) = tracker.csm.state{1}(1);
            newState(2) = tracker.csm.state{2}(1);
            newState(3) = tracker.csm.state{3}(1);
            tracker.states = [tracker.states, newState];
        otherwise
                fprintf('没有运动模型%d\n', dynamicModel);
    end
    tracker.particles(1:3, :) = tracker.particles(1:3, :) + sigma*randn(3, nParticles); % gaussian
    
% 	sigma = 0.008; tracker.particles(1:3, :) = tracker.particles(1:3, :) + sigma*(0.5-rand(3, nParticles)); % uniform
    % begin{测试粒子}
%     tmpx = tracker.particles(1, 1); tmpy = tracker.particles(2, 1); tmpz = tracker.particles(3, 1);
%     
    
%     fig = figure;
%     ax =axes(fig);
%     plot3(1000*tracker.particles(1, :), 1000*tracker.particles(2, :), 1000*tracker.particles(3, :), '.');
%     axis equal
%     grid on 
%     switch dynamicModel
%         case 1
%                 title('None Model');
%         case 2
%                 title('CV Model');
%         case 3
%                 title('CSM Model');
%     end
%     tmp = [-10, 10, -10, 10, -10, 10] + 1000 * [tmpx, tmpx, tmpy, tmpy, tmpz, tmpz];
%     axis(tmp);
%     xlabel('X'); ylabel('Y'); zlabel('Z'); 
% 	prettyPlot(fig, ax, 10, 10, 3);
    % end{测试粒子}
	tracker.weights = repmat(1/nParticles, 1, nParticles); % 初始时，每个粒子的权重相同
    
    % 更新图像模板
    appearance.ellipse{1} = blob1.ellipse;
    appearance.ellipse{2} = blob2.ellipse;
    appearance.intensity{1} = msm(1).image(blob1.idx1);
    appearance.intensity{2} = msm(2).image(blob2.idx1);
    tracker.appearance = appearance;
    tracker.info.idxBlob{1} = [tracker.info.idxBlob{1}, idxOverlapBlob(1)];
    tracker.info.idxBlob{2} = [tracker.info.idxBlob{2}, idxOverlapBlob(2)];
    
    newTracker = tracker;
end

function newTrackers = jobTrack(msm, trackers, dynamicModel, t, stereoModel, cube, morphModel, sigma)
	for i = 1:1:length(trackers)
        tracker = trackers(i);
        newTracker = singlyTrack(tracker, msm, dynamicModel, t, stereoModel, cube, morphModel, sigma);
        newTrackers(i) = newTracker;
	end
end