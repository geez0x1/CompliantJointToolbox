function [ K_lqr, V, Aex, Bex, Cex] = mask_LQR_Full_Measured_State( jointObj, Q, R )
%MASK_LQR_DIRECT Summary of this function goes here
%   Detailed explanation goes here

%% Get state-space model
sys     = jointObj.getStateSpace();
sys     = ss(sys.A, sys.B, eye(size(sys.A)), 0);

% Shorthands
A       = sys.A;
B       = sys.B;
C       = sys.C;
D       = sys.D;

% Create system with 1 output
Ac   	= A;
Bc   	= B;
Cc    	= jointObj.k_b*C(2,:);
Dc   	= 0;
sysc   	= ss(Ac, Bc, Cc, Dc);


%% Design LQR controller
% Calculate LQR gain matrix K
[K_lqr, S, e] = lqr(sysc, Q, R);

Aex = Ac-Bc*K_lqr;

% Design premultiplication.
V = 1 / (Cc*inv(eye(size(Aex)) - Aex )*B);

Bex = Bc*V;
Cex = Cc;

end

