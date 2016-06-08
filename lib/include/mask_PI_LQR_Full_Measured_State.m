function [ K_lqr, V, Ar, Br, Cr] = mask_PI_LQR_Full_Measured_State( jointObj, Q, R )
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
% [K_lqr, S, e] = lqr(sysc, Q, R);

% Aex = Ac-Bc*K_lqr;
% % Design premultiplication.
% V = 1 / (Cc*inv(eye(size(Aex)) - Aex )*B);
% Bex = Bc*V;
% Cex = Cc;

% Augmented plant: here we augment by an integral action
r = size(Cc,1); % Number of controlled outputs
Ar = zeros(r, r); % Controller system matrix
Br = eye(r);      % Controller input matrix
Cr = eye(r);
Er = -1;

Ai = [    Ac,  zeros(size(Ac,1), r);
       Br*Cc,                    Ar];
Bi = [Bc; zeros(r, size(Bc,2))];
Ci = [Cc, zeros(size(Cc,1), r)];

sysi = ss(Ai,Bi,Ci,0);

[K_lqr, S, e] = lqr(sysi, Q, R);

% Cie = [                  C, zeros(r, r);
%        zeros(r, size(C,2)), eye(r)      ];


Aex = Ac-Bc*K_lqr(1:end-1);
% Design premultiplication.
V = 1 / (Cc*inv(eye(size(Aex)) - Aex )*Bc);
Bex = Bc*V;
Cex = Cc;

end

