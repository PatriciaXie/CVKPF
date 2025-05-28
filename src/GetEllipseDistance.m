function ellipseDistance = GetEllipseDistance(ellipse, ellipse_)
% �������ķ�������������Բ�����غϣ���ƽ�Ƶ�ԭ�㣩����ellipse�ϲ���30���㣬���ÿ���㴦�ķ�����ellipse_�Ľ��㣬����30�Ե�ľ���֮��
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