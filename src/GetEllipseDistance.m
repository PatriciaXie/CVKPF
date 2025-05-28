function ellipseDistance = GetEllipseDistance(ellipse, ellipse_)
% 程昔恩的方法，先令两椭圆中心重合（都平移到原点），在ellipse上采样30个点，求得每个点处的法线与ellipse_的交点，计算30对点的距离之和
    nSample = 30+1;
    sampleAlpha = linspace(0, 2*pi, nSample);
    try 
        samplePoints = sqrtm(ellipse.covariance)*[cos(sampleAlpha); sin(sampleAlpha)];
    catch
        a = 1;
    end
    

    samplePoints_ = zeros(2, nSample);
    for i = 1:1:nSample-1
        p0 = samplePoints(:, i+1); p = samplePoints(:, i);
        gradientVector = (p0 - p)/norm(p0 - p);
        normalVector = [0, 1; -1, 0] * gradientVector;
        line.p0 = p0;
        line.gradientVector = gradientVector;
        line.d = normalVector;
%         p0_ = LineCrossEllipse1(ellipse_, line);
        p0_ = LineCrossEllipse2(ellipse_, line);
        samplePoints_(:, i+1) = p0_;
    end
    samplePoints_(:, 1) = samplePoints_(:, end);
    
%     showEllipseDistance(samplePoints, samplePoints_);
    
    ellipseDistance = 0;
    for i = 2:1:nSample
        p0 = samplePoints(:, i); p0_ = samplePoints_(:, i); 
        ellipseDistance = ellipseDistance + norm(p0 - p0_);
    end
    
end