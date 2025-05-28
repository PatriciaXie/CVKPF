function nccValue = GetNcc(srcImg1, targetImg1)
% 计算一维图像的相关性
    lenSrc = length(srcImg1);
    lenTarget = length(targetImg1);
    meanSrc = mean(srcImg1);
    deltaLen = lenSrc - lenTarget;
    
    if deltaLen > 0
        meanTarget = mean(targetImg1);
        targetImg1 = [targetImg1, repmat(meanTarget, 1, deltaLen)];
    else
        if deltaLen < 0
            meanSrc = mean(srcImg1);
            srcImg1 = [srcImg1, repmat(meanSrc, 1, -deltaLen)];
        end
    end
    nccValue = ncc(srcImg1, targetImg1);
end

function m = ncc(w1, w2)
	denom = sqrt( sum(sum(w1.^2))*sum(sum(w2.^2)) );
	if denom < 1e-10
		m = 0;
	else
		m = sum(sum((w1.*w2))) / denom;
	end
end