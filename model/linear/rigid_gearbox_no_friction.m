% RIGID_GEARBOX_NO_FRICTION Get linear dynamics matrices - rigid gearbox without
% friction
%
% [A, B, C, D, I, R, K] = jointObj.rigid_gearbox_no_friction
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
%  This function is identical to rigid_gearbox, but with the
%  difference that the bearing friction coefficients are set to zero.
%
% Examples::
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also full_dyn, rigid_gearbox.

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
function [A, B, C, D, I, R, K] = rigid_gearbox_no_friction(obj)
    
    % The computations below assume a state vector definition according to:
    % x = [q_g, q_l, q_g_dot, q_l_dot]', where 
    % q_g is the gearbox output angle
    % q_l is the flange angle (output of the torsion bar)
    %
    % The '_dot' denotes the temporal derivative.

    % Inertia matrix
    I = diag([obj.I_m + obj.I_g, obj.I_l]);

    % Damping matrix
    d_m     = 0 * obj.d_m;
    d_l     = 0 * obj.d_l;
    d_gl    = obj.d_gl; % shorthands %#ok<*PROP>
    R = [   d_m + d_g + d_gl,       -d_gl;
            -d_gl,                  d_l + d_gl      ];

    % Stiffness matrix
    k_b = obj.k_b; % shorthands %#ok<*PROP>
    K = [   k_b,        -k_b;
            -k_b,       k_b         ];

    % State-space matrices
    A = [   zeros(size(I)),     eye(size(I));
            -I\K,               -I\R            ];

    % Input
    % u = [tau_m, tau_e]
    k_t = obj.k_t;
    n   = obj.n;
    B   = [ 0,              0;
            0,              0;
            k_t*n/I(1,1),   0;
            0,              1/I(2,2)    ];
    
    % Output
    C = [   1,      0,      0, 0;           % motor position
            1,      0,      0, 0;           % gear position
            0,      1,      0, 0;           % link position
            0,      0,      1, 0;           % motor velocity
            0,      0,      1, 0;           % gear velocity
            0,      0,      0, 1;           % link velocity
            k_b,    -k_b,   d_gl, -d_gl	];	% Torsion bar torque
     
    % Direct Feedthrough
    nIn = size(B,2);
    nOut = size(C,1);
    D = zeros(nOut,nIn);
end
