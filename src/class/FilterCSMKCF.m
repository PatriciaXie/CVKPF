classdef FilterCSMKCF
    %FilterCSMKCF �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
        alpha % ����Ƶ��
        aMax % �����ٶ�
        T % �������
        Q1 % q11~q33���ɵľ���
        U % �������
        F % ״̬ת�ƾ���
        H % �۲����
        R % �۲���������
    end
    
    methods
        function this = FilterCSMKCF(alpha, aMax, R, T)
            %FilterCSMKCF ��������ʵ��
            %   �˴���ʾ��ϸ˵��
            this.alpha = alpha;
            this.aMax = aMax;
            
            q11 = 1/(2*alpha^5)*(1-exp(-2*alpha*T)+2*alpha*T+2*alpha^3*T^3/3-2*alpha^2*T^2-4*alpha*T*exp(-alpha*T));
            q12 = 1/(2*alpha^4)*(exp(-2*alpha*T)+1-2*exp(-alpha*T)+2*alpha*T*exp(-alpha*T)-2*alpha*T+alpha^2*T^2);
            q13 = 1/(2*alpha^3)*(1-exp(-2*alpha*T)-2*alpha*T*exp(-alpha*T));
            q22 = 1/(2*alpha^3)*(4*exp(-alpha*T)-3-exp(-2*alpha*T)+2*alpha*T);
            q23 = 1/(2*alpha^2)*(exp(-2*alpha*T)+1-2*exp(-alpha*T));
            q33 = 1/(2*alpha)*(1-exp(-2*alpha*T));
            this.Q1 = [q11, q12, q13; q12, q22, q23; q13, q23, q33];
            
            u1 = 1/alpha*(-T+alpha*T^2/2+(1-exp(-alpha*T))/alpha);
            u2 = T-(1-exp(-alpha*T))/alpha;
            u3 = 1-exp(-alpha*T);
            this.U = [u1; u2; u3];
            
            this.F =[1, T, (alpha*T-1+exp(-alpha*T))/alpha^2;
                          0, 1,  (1-exp(-alpha*T))/alpha;
                          0, 0, exp(-alpha*T)];
            this.H = [1, 0, 0];
            this.R = R;
        end
        
        function [X, P, Xm, Pm, aAvg] = Filter(this, X_, P_, aAvg_, z)
            %METHOD1 �˴���ʾ�йش˷�����ժҪ
            %   �˴���ʾ��ϸ˵��
            if aAvg_ > 0
                sigma2 = (4-pi)/pi * (this.aMax - aAvg_)^2;
            else
                sigma2 = (4-pi)/pi * (-this.aMax - aAvg_)^2;
            end
            Q = 2 * this.alpha * sigma2 * this.Q1; % CSMKCF
%             Q = diag([1,1,1]); % KCF
            Xm = this.F * X_ + this.U * aAvg_;
            zh = this.H * Xm;
            dz = z - zh;
            Pm = this.F * P_ * (this.F)' + Q;
            K = Pm * this.H' * inv(this.H*Pm*(this.H)' + this.R);
            X = Xm + K*dz;
            P = Pm - K * this.H * Pm;
            aAvg = X(3);
        end
    end
end

