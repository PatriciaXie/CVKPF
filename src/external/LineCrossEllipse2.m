function p0_ = LineCrossEllipse2(ellipse_, line)

    lp1 = line.p0(1); lp2 = line.p0(2); d1 = line.d(1); d2 = line.d(2); 
    a = ellipse_.e(1,1); b = ellipse_.e(2,2); c = ellipse_.e(1,2); ep1 = 0; ep2 = 0;
    lambda = [((c^2*d1^2*ep2^2 - 2*c^2*d1^2*ep2*lp2 + c^2*d1^2*lp2^2 - 2*c^2*d1*d2*ep1*ep2 + 2*c^2*d1*d2*ep1*lp2 + 2*c^2*d1*d2*ep2*lp1 - 2*c^2*d1*d2*lp1*lp2 + c^2*d2^2*ep1^2 - 2*c^2*d2^2*ep1*lp1 + c^2*d2^2*lp1^2 + 2*c*d1*d2 - a*b*d1^2*ep2^2 + 2*a*b*d1^2*ep2*lp2 - a*b*d1^2*lp2^2 + a*d1^2 + 2*a*b*d1*d2*ep1*ep2 - 2*a*b*d1*d2*ep1*lp2 - 2*a*b*d1*d2*ep2*lp1 + 2*a*b*d1*d2*lp1*lp2 - a*b*d2^2*ep1^2 + 2*a*b*d2^2*ep1*lp1 - a*b*d2^2*lp1^2 + b*d2^2)^(1/2) + a*d1*ep1 + b*d2*ep2 + c*d1*ep2 + c*d2*ep1 - a*d1*lp1 - b*d2*lp2 - c*d1*lp2 - c*d2*lp1)/(a*d1^2 + 2*c*d1*d2 + b*d2^2);
              -((c^2*d1^2*ep2^2 - 2*c^2*d1^2*ep2*lp2 + c^2*d1^2*lp2^2 - 2*c^2*d1*d2*ep1*ep2 + 2*c^2*d1*d2*ep1*lp2 + 2*c^2*d1*d2*ep2*lp1 - 2*c^2*d1*d2*lp1*lp2 + c^2*d2^2*ep1^2 - 2*c^2*d2^2*ep1*lp1 + c^2*d2^2*lp1^2 + 2*c*d1*d2 - a*b*d1^2*ep2^2 + 2*a*b*d1^2*ep2*lp2 - a*b*d1^2*lp2^2 + a*d1^2 + 2*a*b*d1*d2*ep1*ep2 - 2*a*b*d1*d2*ep1*lp2 - 2*a*b*d1*d2*ep2*lp1 + 2*a*b*d1*d2*lp1*lp2 - a*b*d2^2*ep1^2 + 2*a*b*d2^2*ep1*lp1 - a*b*d2^2*lp1^2 + b*d2^2)^(1/2) - a*d1*ep1 - b*d2*ep2 - c*d1*ep2 - c*d2*ep1 + a*d1*lp1 + b*d2*lp2 + c*d1*lp2 + c*d2*lp1)/(a*d1^2 + 2*c*d1*d2 + b*d2^2)];
   cross = [line.p0 + lambda(1)*line.d, line.p0 + lambda(2)*line.d];
   p0_ = real(cross(:, 1)');
end