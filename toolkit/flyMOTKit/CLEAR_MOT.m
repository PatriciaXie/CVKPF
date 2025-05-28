function [recall, precision, IDSW, MOTA, MOTP] = CLEAR_MOT(gtMat, rsMat, threshold)
% compute CLEAR MOT
%
% metrics contains the following
% [1]   recall	- percentage of detected targets
% [2]   precision	- percentage of correctly detected targets
% [3]  IDSW	- number of id switches     (IDs)
% [4]  MOTA	- Multi-object tracking accuracy in [0,100]
% [5]  MOTP	- Multi-object tracking precision in [0,100] (3D) / [td,100] (2D)

% 把ID调整成连续的，从1开始的
[~, ~, ic] = unique(gtMat(:,2)); gtMat(:,2) = ic;
[~, ~, ic2] = unique(rsMat(:,2)); rsMat(:,2) = ic2;

%% 生成每一帧的匹配
frame_gt = unique(gtMat(:, 1)); % Gt涉及的帧号
% frame_res = unique(resMat(:, 1)); % Res涉及的帧号
% frame_min = min(min(frame_gt), min(frame_res)); % 两者的最大帧号
% frame_max = max(max(frame_gt), max(frame_res)); % 两者的最小帧号
pairsAll = cell(length(frame_gt), 1);

Gt = zeros(length(frame_gt), 1); % GT每一帧的目标数
Rt = zeros(length(frame_gt), 1); % Res每一帧的目标数
Mt = zeros(length(frame_gt), 1); % 每一帧的匹配数

frames = frame_gt;
for k = 1:1:length(frame_gt)
    iframe = frames(k);
    dataGt = gtMat(gtMat(:, 1) == iframe, :); Gt(k) = size(dataGt, 1);
    dataRes = rsMat(rsMat(:, 1) == iframe, :); Rt(k) = size(dataRes, 1);
    idsGt = dataGt(:, 2);
    idsRes = dataRes(:, 2);
    if iframe>frame_gt(1)
        lastPairs = pairsAll{k-1};
    else
        lastPairs = [];
    end
    pairs = [];
    
    % 进行轨迹演化
    if ~isempty(lastPairs)
        nPairs = size(lastPairs, 1);
        for i = 1:1:nPairs
            id1 = lastPairs(i, 1); id2 = lastPairs(i, 2);
            line1 = dataGt(dataGt(:, 2) == id1, :);
            line2 = dataRes(dataRes(:, 2) == id2, :);
            if ~isempty(line1) && ~isempty(line2)
                pos1 = line1(3:5); pos2 = line2(3:5);
                distance = norm(pos1 - pos2);
                if distance <= threshold
                    pairs = [pairs; id1, id2, distance];
                end
            end
        end
    end
    if ~isempty(pairs)
        idsLeftGt = setdiff(idsGt, pairs(:, 1));
        idsLeftRes = setdiff(idsRes, pairs(:, 2));
    else
        idsLeftGt = idsGt;
        idsLeftRes = idsRes;
    end
    
    % 不能演化的，生成cost矩阵，进行匈牙利匹配
    nLeftGt = length(idsLeftGt); nLeftRes = length(idsLeftRes);
    costs = zeros(nLeftGt, nLeftRes);
    for i = 1:1:nLeftGt
        id1 = idsLeftGt(i); 
        line1 = dataGt(dataGt(:, 2) == id1, :);
        pos1 = line1(3:5);
        for j = 1:1:nLeftRes
            id2 = idsLeftRes(j);
            line2 = dataRes(dataRes(:, 2) == id2, :);
            pos2 = line2(3:5);
            costs(i, j) = norm(pos1 - pos2);
        end
    end
    [Matching, Cost] = Hungarian(costs);
    for i = 1:1:nLeftGt
        id1 = idsLeftGt(i);
        tmp = find(Matching(i,:));
        if ~isempty(tmp)
            id2 = idsLeftRes(tmp);
            if costs(i, tmp) <= threshold
                pairs = [pairs; id1, id2, costs(i, tmp)];
            end
        end
    end
    Mt(k) = size(pairs, 1);
    pairsAll{k} = pairs;
end

%% 分析匹配
TPt = zeros(length(frame_gt), 1); 
% TP并不等于Mt，是匈牙利匹配的pair，但考虑了时间连续性后，不一定是TP，如Fig.9(c)
FPt = zeros(length(frame_gt), 1);
FNt = zeros(length(frame_gt), 1);
TPdt = zeros(length(frame_gt), 1); % 每一帧所有TP的总距离
IDSWt = zeros(length(frame_gt), 1);

%% 正确地计算TPt
% 以Res轨迹为参考，形成tracks
% idsRes = unique(rsMat(:, 2));
% tracksRes = cell(length(idsRes), 1);
% for i = 1:1:length(idsRes)
%     trackRes = zeros(length(frame_gt), 3); % id1, id2
%     id2 = idsRes(i);
%     for k = 1:1:length(pairsAll)
%         pairs = pairsAll{k};
%         if ~isempty(pairs)
%             tmp = pairs(pairs(:, 2) == id2, :);
%             if ~isempty(tmp)
%                 trackRes(k, :) = [frames(k), tmp(1:2)];
%             end
%         end
%     end
%     
%     trackRes(trackRes(:, 3)==0, :) = [];
%     tracksRes{i} = trackRes;
% end
% 
% % 分析track(:, 2)，若track(:, 2)保持不变，则是TP，否则后面的是FP
% for i = 1:1:length(idsRes)
%     trackRes = tracksRes{i};
%     if isempty(trackRes)
%         continue;
%     end
% 	id1_old = trackRes(1, 2);
%     for k = 1:1:size(trackRes, 1)
%         id1 = trackRes(k, 2);
%         frame = trackRes(k, 1);
%         if id1 == id1_old
%             TPt(frames == frame) = TPt(frames == frame) + 1;
%         else
%             FPt(frames == frame) =  FPt(frames == frame) + 1;
%         end
%     end    
% end

%% 正确地计算IDSWt
% 以Gt轨迹为参考，形成tracks
empTracks = [];
idsGt = unique(gtMat(:, 2));
kk = 0;
for i = 1:1:length(idsGt)
    trackGt = zeros(length(frame_gt), 3); % id1, id2
    id1 = idsGt(i);
    for k = 1:1:length(pairsAll)
        pairs = pairsAll{k};
        if ~isempty(pairs)
            tmp = pairs(pairs(:, 1) == id1, :);
            if ~isempty(tmp)
                trackGt(k, :) = [frames(k), tmp(1:2)];
            end
        end
    end
    trackGt(trackGt(:, 2)==0, :) = [];
    if isempty(trackGt)
        empTracks = [empTracks, i];
    else
        kk = kk + 1;
       tracksGt{kk} = trackGt;
    end
end
if ~isempty(empTracks)
    idsGt(empTracks) = [];
end
        
% 分析track(:, 3)有几次跳变，就是几次ID Switch
for i = 1:1:length(idsGt)
    trackGt = tracksGt{i};
    id2_old = trackGt(1, 3);
    for k = 2:1:size(trackGt, 1)
        id2 = trackGt(k, 3);
        if id2 ~= id2_old
            id2_old = id2;
            frame = trackGt(k, 1);
            IDSWt(frames == frame) =  IDSWt(frames == frame) + 1;
        end
    end    
end


%% 指标计算
for k = 1:1:length(pairsAll)
    TPt(k) = size(pairsAll{k}, 1);
    FPt(k) = Rt(k) - TPt(k);
    FNt(k) = Gt(k) - TPt(k);
    pairs = pairsAll{k};
    if ~isempty(pairs)
        TPdt(k) = sum(pairs(:, 3)); 
    else
        TPdt(k) = 0;
    end
end

TP = sum(TPt); FP = sum(FPt); FN = sum(FNt); G = sum(Gt);
IDSW=sum(IDSWt);
MOTP=sum(TPdt)/TP * 1000; % mm
MOTA=100 * (1-(FN+FP+IDSW)/G);
recall = 100 * TP/G;
precision = 100* TPt/(FP+TP);
end 
