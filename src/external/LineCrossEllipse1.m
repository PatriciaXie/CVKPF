function p0_ = LineCrossEllipse1(ellipse_, line)
    p0 = line.p0;
    gradientVector = line.gradientVector;
    syms x y
    a_ = ellipse_.a; b_ = ellipse_.b; theta_ = ellipse_.theta;
    equationEllipse = (x * cos(theta_) + y * sin(theta_))^2 / a_^2 + (-x * sin(theta_) + y * cos(theta_))^2 / b_^2 == 1;
    equationLine = [x-p0(1), y-p0(2)] * gradientVector == 0;
    sol = solve(equationEllipse, equationLine, x, y); % 直线和椭圆可能没交点吗？应该不可能
    xs = real(double(sol.x)); ys = real(double(sol.y));
    p1 = [xs(1); ys(1)]; p2 = [xs(2); ys(2)];
    d1 = norm(p1-p0); d2 = norm(p2-p0);
    if d1 < d2
        p0_ = p1;
    else
        p0_ = p2;
    end
end