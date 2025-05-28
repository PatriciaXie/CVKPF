function [idxes, ratios] = GetOverlap(camProj, imgSize, state, blobs, shape, model)
%GetOverlap 计算以point3d为球心的球在相机上的投影与blobs的重叠情况
%   Input
%       camProj 相机投影矩阵
%       state [x,y,z, phi, theta,  x_,y_,z_] 状态
%       blobs 待匹配量测
sphereModel = model{1};
flyModel = model{2};

idxes = []; ratios = [];
if isempty(blobs)
    return;
end
% 计算球在相机上的投影

if strcmp('EllipsoidShape', shape)
    ellipsoid = GetShape(flyModel, state);
    pixels2 =camProj * [ellipsoid; ones(1,size(ellipsoid,2))];
else
    point3d = state(1:3);
    sphere = sphereModel + point3d;
    pixels2 =camProj * [sphere; ones(1,size(sphere,2))];
end
pixels2(1:2, :) = pixels2(1:2, :) ./ repmat(pixels2(3, :), 2, 1);
pixels2 = floor(pixels2(1:2,:));
pixels2 = unique(pixels2','rows')';
% img = zeros(imgSize);
% pixels1 = sub2ind(imgSize, pixels2(2,:), pixels2(1, :));
% img(pixels1) = 1;
% imshow(img);
    
% 判断投影是否完全在相平面上
flag1 = isempty(find(pixels2(1,:)>imgSize(1), 1));
flag2 = isempty(find(pixels2(1,:)<1, 1));
flag3 = isempty(find(pixels2(1,:)>imgSize(2), 1));
flag4 = isempty(find(pixels2(1,:)<1, 1));

if flag1 && flag2 && flag3 && flag4
    pixels1 = sub2ind(imgSize, pixels2(2,:), pixels2(1, :));
%     showRePorjection(pixels1, blobs, imgSize); % 红是粒子，绿是量测

    regionLimit = cat(2, blobs.idxRange);
    tmpIdx = find((regionLimit(1,:)>=min(pixels1)) & (regionLimit(2,:)<=max(pixels1)));
    for i = 1:1:length(tmpIdx)
        idx = tmpIdx(i);
        overlapped = blobs(idx).idx1(ismember(blobs(idx).idx1, pixels1));
        if ~isempty(overlapped)
            idxes = [idxes, idx];
            ratios = [ratios, length(overlapped) / blobs(idx).area];
        end
    end
    
end

end

