%FULL_DYN Get linear dynamics matrices - default
%
% [A, B, C, D, I, R, K] = jointObj.full_dyn
%
% jointObj is the instance of the joint class object for which this
% function has been called.
%
% Outputs::
%   A:   System matrix
%   B:   Input matrix
%   C:   Output matrix
%   D:   Direct Feedthrough matrix
%   I:   Inertia matrix
%   R:   Damping matrix
%   K:   Stiffness matrix
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
% See also output_fixed, rigid_gearbox.

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


function [A, B, C, D, I, R, K] = full_dyn(obj)
    
    % The computations below assume a state vector definition according to:
    % x = [q_m, q_g, q_l, q_m_dot, q_g_dot, q_l_dot]', where 
    % q_m is the motor angle,
    % q_g is the gearbox output angle
    % q_l is the flange angle (output of the torsion bar)
    %
    % The '_dot' denotes the temporal derivative.

    % NOTE
    % The size of I, R, K needs to be identical, square, and the number of
    % rows needs to equal the number of derivatives in the state x.
    % This is important for the construction of B_nonlinear in
    % mask_JointMechanicalSubsystem.
    
    % Inertia matrix
    I = diag([obj.I_m, obj.I_g, obj.I_l]);

    % Damping matrix
    d_m     = obj.d_m;
    d_g     = obj.d_g;
    d_l     = obj.d_l;
    d_mg    = obj.d_mg;
    d_gl    = obj.d_gl; % shorthands %#ok<*PROP>
    R = [   d_m + d_mg,     -d_mg,                  0;
            -d_mg,          d_g + d_mg + d_gl,      -d_gl;
            0,              -d_gl,                  d_l + d_gl  ];

    % Stiffness matrix
    k_g = obj.k_g;
    k_b = obj.k_b; % shorthands %#ok<*PROP>
    K = [   k_g,            -k_g,           0;
            -k_g,           k_g + k_b,      -k_b;
            0,              -k_b,           k_b         ];

    % State-space matrices
    A = [   zeros(size(I)),     eye(size(I));
            -I\K,               -I\R            ];

    % Input
    % u = [tau_m, tau_e]
    B   = [ 0,          0;
            0,          0;
            0,          0;
            1/I(1,1),   0;
            0,          0;
            0,          1/I(3,3)	];
    
    % Output
    C = [   eye(size(A,2));                     % All states
            0, k_b, -k_b, 0, d_gl, -d_gl	];	% Torsion bar torque
        
    % Direct Feedthrough
    nIn = size(B,2);
    nOut = size(C,1);
    D = zeros(nOut,nIn);

end
