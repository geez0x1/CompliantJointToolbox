function [ K_lqr, V, Ar, Br, Cr] = mask_PI_LQR_Full_Measured_State( jointObj, Q, R )
%MASK_PI_LQR_FULL_MEASURED_STATE Mask code for the PI-LQR with full
%measured state.
%
% [ K, V, A, B, C ] = mask_PI_LQR_Full_Measured_State(jointObj, Q, R )
%
% Inputs::
% jointObj: Object defining the model to use.
% Q:        Diagonal weighting matrix
% R:        Diagonal weighting matrix
%
% Outputs::
% K:        Gain matrix
% V:        Premultiplication
% A:        Closed loop system matrix
% B:        Closed loop input matrix
% C:        Closed loop output matrix
%
% Notes::
%
%
% Examples::
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also mask_LQR_Full_Measured_State.

% Copyright (C) 2016, by Joern Malzahn, Wesley Roozing
%
% This file is part of the Compliant Joint Toolbox (CJT).
%
% CJT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% CJT is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
% License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CJT. If not, see <http://www.gnu.org/licenses/>.
%
% For more information on the toolbox and contact to the authors visit
% <https://github.com/geez0x1/CompliantJointToolbox>

% Get state-space model
sys     = jointObj.getStateSpace();

% Shorthands
A       = sys.A;
B       = sys.B(:,1); % Use only current input
C       = sys.C;
D       = sys.D;

% Get state-space model
sys     = ss(A, B, eye(size(sys.A)), 0);

% Create system with 1 output
Ac      = A;
Bc      = B;
Cc      = jointObj.k_b * C(2,:);
Dc      = 0;
sysc    = ss(Ac, Bc, Cc, Dc);


%% Design LQR controller
% Calculate LQR gain matrix K
% [K_lqr, S, e] = lqr(sysc, Q, R);

% Aex = Ac-Bc*K_lqr;
% % Design premultiplication.
% V = 1 / (Cc*inv(eye(size(Aex)) - Aex )*B);
% Bex = Bc*V;
% Cex = Cc;

% Augmented plant: here we augment by an integral action
r   = size(Cc,1);       % Number of controlled outputs
Ar  = zeros(r, r);      % Controller system matrix
Br  = eye(r);           % Controller input matrix
Cr  = eye(r);
Er  = -1;

Ai  = [ Ac,     zeros(size(Ac,1), r); ...
        Br*Cc,  Ar                      ];
Bi  = [Bc; zeros(r, size(Bc,2))];
Ci  = [Cc, zeros(size(Cc,1), r)];

sysi = ss(Ai,Bi,Ci,0);

[K_lqr, S, e] = lqr(sysi, Q, R);

% Cie = [                  C, zeros(r, r);
%        zeros(r, size(C,2)), eye(r)      ];


Aex = Ac - Bc * K_lqr(1:end-1);

% Design premultiplication
V   = 1 / (Cc * inv(eye(size(Aex)) - Aex) * Bc);
Bex = Bc * V;
Cex = Cc;

end

