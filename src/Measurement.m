classdef Measurement
    %Blobs �����ʾ
    %   �˴���ʾ��ϸ˵��
    
    properties
        blobs
        image
        imageSize
    end
    
    methods
        function msm = Measurement(image, bgm)
            %Measurement ��������ʵ��
            %   �˴���ʾ��ϸ˵��
            msm.image = image;
            msm.imageSize = size(image);
            bgMean = bgm.meanValue;
            bgStd = bgm.stdValue;
            % �������ݲ�������
%             mask = abs((image-bgMean) ./ bgStd) > 9;
            tmp = abs((image-bgMean) ./ bgStd); tmp = tmp/min(9, max(tmp(:)));
%             mask =  tmp > 0.95;
%             se = strel('disk',5);
            mask = imdilate(tmp > 0.95, strel('disk',5));
            mask(1:10, :) = 0;
            % ��ɸ
            minArea = 10;
            tmpBlobs = analyseMask(mask, minArea);
            
            % ��ɸ
            mask = zeros(size(image));
            stdError = 1; %1.7;
            for i = 1:1:length(tmpBlobs)
                idx1 = tmpBlobs(i).idx1;
                blobIntensity = image(idx1);
                blobMean = mean(blobIntensity);
                blobStd = std(blobIntensity);
                validIdx = blobIntensity <= blobMean - blobStd * stdError;
                idx1 = idx1(validIdx);
                if length(idx1) > 10 && length(idx1) < 200
                    mask(idx1) = 1;
                end
            end
            minArea = 10;
            tmpBlobs = analyseMask(logical(mask), minArea);
            
            % �����Բ
            weightMask = abs((image-bgMean) ./ bgStd);
            try
                tmpBlobs = analyseBlobs(tmpBlobs, weightMask);
            catch
                a = 1;
            end
            msm.blobs = tmpBlobs;
        end
        
    end
    
end


function blobs = analyseMask(mask, minArea)
    stats = regionprops(mask, 'Area', 'PixelList', 'Centroid');
    i = 1;
    for n=1:length(stats)
        if stats(n).Area>minArea
            tmpBlobs(i).center = stats(n).Centroid';
            tmpBlobs(i).isOccupied = false;
            tmpBlobs(i).area = stats(n).Area;
            tmpBlobs(i).idx2 = stats(n).PixelList'; %��һ�����У��ڶ�������
            tmpBlobs(i).idx1 = sub2ind(size(mask), tmpBlobs(i).idx2(2,:), tmpBlobs(i).idx2(1, :));
            tmpBlobs(i).idxRange = [max(tmpBlobs(i).idx1); min(tmpBlobs(i).idx1)];
            i = i+1;
        end
    end
    if i <= 1
        blobs = [];
    else
        blobs = tmpBlobs;
    end
end


function blobs = analyseBlobs(blobs, weightMask)
    emptyIdx = [];
    for i = 1:1:length(blobs)
        idx1 = blobs(i).idx1;
        weights = weightMask(idx1);
        N = size(blobs(i).idx2, 2);
        W = sum(weights); 
        if W == 0
            emptyIdx = [emptyIdx, i];
            continue
        end
        meanValue = sum(blobs(i).idx2 .* repmat(weights, 2, 1), 2) / W;
        tmp = blobs(i).idx2 - repmat(meanValue,1, N);
        covariance = tmp .* repmat(weights, 2, 1) * tmp' / W;
        covariance = covariance * 5; % ΪʲôҪ����5��
        [ve, va] = eig(covariance);
        [vaSorted, idxOrg] = sort(diag(va), 'ascend');
        veSorted = ve(:,idxOrg);
        radii = sqrt(vaSorted);
        ellipse.center = meanValue;
        ellipse.a = radii(2);
        ellipse.b = radii(1);
        ellipse.axis = veSorted;
        ellipse.theta = atan2(veSorted(2,2),veSorted(1,2)); % ���᷽������
        ellipse.covariance = covariance;
        ellipse.e = pinv(covariance); % inv, pinv
        ellipse.factor = pi * ellipse.a * ellipse.b / length(idx1);
        blobs(i).ellipse = ellipse;
    end
    blobs(emptyIdx) = [];
end
