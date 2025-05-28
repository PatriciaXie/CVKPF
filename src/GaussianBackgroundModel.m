classdef GaussianBackgroundModel
    properties
        meanValue
        stdValue
        imageSeries
        startFrame
        nFrame
    end
    
    methods
        function bgm = GaussianBackgroundModel(imageSeries)
            bgm.meanValue = mean(imageSeries, 3);
            tmpStd = std(imageSeries, 0, 3);
            tmpStd(tmpStd==0) = 1;
            bgm.stdValue = tmpStd;
            bgm.imageSeries = imageSeries;
            bgm.startFrame = 0;
            bgm.nFrame = size(imageSeries, 3);
        end
        
        function bgm = Update(bgm, newImage)
            idx = mod(bgm.startFrame, bgm.nFrame) + 1;
            tmpImageSeries = bgm.imageSeries;
            tmpImageSeries(:, :, idx) = newImage;
            bgm.meanValue = mean(tmpImageSeries, 3);
            tmpStd = std(tmpImageSeries, 0, 3);
            tmpStd(tmpStd==0) = 1;
            bgm.stdValue = tmpStd;
            bgm.imageSeries = tmpImageSeries;
            bgm.startFrame = bgm.startFrame + 1;
        end
    end
end