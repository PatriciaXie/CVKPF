function showEllipseDistance(samplePoints, samplePoints_)
    figure
    hold on
    plot(0,0,'k*');
    plot(samplePoints(1,:), samplePoints(2, :), 'k-o');
    plot(samplePoints_(1,:), samplePoints_(2, :), 'r-o');
    legend('≤ŒøºÕ÷‘≤t', 'Ωªµ„t-1');
    axis equal
    grid on
    hold off
end