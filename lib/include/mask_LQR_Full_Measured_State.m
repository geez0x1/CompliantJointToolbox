function [ K_lqr, V, Aex, Bex, Cex] = mask_LQR_Full_Measured_State( jointObj, Q, R )
%MASK_LQR_DIRECT Mask code for the LQR controller block.
%
% [K, V, A, B, C ] = mask_LQR_Full_Measured_State(jointObj, Q, R )
%
% Inputs::
% jointObj: Object defining the model to use.
% Q:        Diagonal weighting matrix
% R:        Diagonal weighting matrix
%
% Outputs::
% K:    Gain matrix
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
% See also mask_PI_LQR_Full_Measured_State.

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

