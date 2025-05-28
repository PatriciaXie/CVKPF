function orientation = ReconstructOrientation(CamProj1, CamProj2, ellipse1, ellipse2)
%ReconstructOrientation ˫Ŀ�ؽ����峯��
%   Input
%       ellipse1�����1�ϵ���Բ
%       ellipse2�����2�ϵ���Բ
%   Output
%       orientation����������

    major(:,1) = ellipse1.center - ellipse1.axis(:,2) * ellipse1.a;
    major(:,2) = ellipse1.center + ellipse1.axis(:,2) * ellipse1.a;
    rayActual1 = GenerateRay(CamProj1,major(:,1));
    rayActual2 = GenerateRay(CamProj1,major(:,2));

    major(:,1) = ellipse2.center - ellipse2.axis(:,2) * ellipse2.a;
    major(:,2) = ellipse2.center + ellipse2.axis(:,2) * ellipse2.a;
    rayDesire1 = GenerateRay(CamProj2,major(:,1));
    rayDesire2 = GenerateRay(CamProj2,major(:,2));

    normActual = cross(rayActual1.d, rayActual2.d); 
    normDesire = cross(rayDesire1.d,rayDesire2.d);
    orientation = cross(normActual,normDesire); 
    
    orientation = unit(orientation);
    
    orientation = RectifyOrientation(orientation);
end

function ray3d = GenerateRay(CamProj, d2Point)
    Mi = inv(CamProj(1:3,1:3)); p4 = CamProj(:,4);
    ray3d.P0 = -Mi*p4;
    ray3d.d = unit(Mi*e2h(d2Point));
end


function d3Orientation = RectifyOrientation(d3Orientation)
    [~, phi, ~] = cart2sph(d3Orientation(1), d3Orientation(3), -d3Orientation(2)); % �ѿ�������ϵת��Ϊ��������ϵ
    if ( phi < 0 ) 
        d3Orientation = -d3Orientation; 
    end
end