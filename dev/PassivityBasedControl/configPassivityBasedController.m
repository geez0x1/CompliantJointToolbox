jointModel = WMBig10k_rigid_gearbox;%_coulomb;

jointModel.I_b = 1;

% jointModelSym = WMBig10k_rigid_gearbox_coulomb;
% jointModelSym.makeSym;

[A B C D ] = jointModel.getDynamicsMatrices;

p_is_Hz = abs(eig(A))/2/pi;
p_d_Hz = 1.2*[  10.1; 10.2; 10.3; 10.4];%M10*ones(size(p_is_Hz));


p_d_rad = -p_d_Hz*2*pi;


[K,prec,message] = place(A,B(:,1),p_d_rad);
% Q = ctrbf(A,B,C);
% Qinv = inv(Q);
% t1 = Qinv(end,:);

sys = ss(A,B(:,1),C,0);

Kp = K(1)
Kd = K(3)
KT = K(2)
KS = K(4)
KB = jointModel.k_b;
D = jointModel.d_b;
alpha = 1300;

% if (KB < alpha)
%     error('Joint stiffness insufficient!');
% end
% 
% if (KB+KT < 0)
%     error('KT too small!');
% end
% 
% if (Kp < alpha*(KT + KB)/(KB-alpha))
%     error('K_p too small!');
% end
% 
% if (Kd < (KS*KB - KT*D)^2/(4*KB*D*(KB+KT)));
%     error('Passivity condition violated. Increase Kd!');
% end

Ccl = [KB;-KB; 0; 0;].';
clsys = ss(A-B(:,1)*K,B(:,1),Ccl,0)
