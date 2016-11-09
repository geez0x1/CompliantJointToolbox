%FULL_DYN Get linear dynamics matrices - default
%
% [A, B, C, I, D, K] = jointObj.full_dyn
%
% jointObj is the instance of the joint class object for which this
% function has been called.
%
% Outputs::
%   A:   System matrix
%   B:   Input matrix
%   C:   Output matrix
%   I:   Inertia matrix
%   D:   Damping matrix
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
% See also full_dyn_no_friction, output_fixed.

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


function [A, B, C, I, D, K] = full_dyn(obj)
    
    % The computations below assume a state vector definition according to:
    % x = [q_m, q_g, q_l, q_m_dot, q_g_dot, q_l_dot]', where 
    % q_m is the motor angle,
    % q_g is the gearbox output angle
    % q_l is the flange angle (output of the torsion bar)
    %
    % The '_dot' denotes the temporal derivative.

    % Inertia matrix
    I = diag([obj.I_m, obj.I_g, obj.I_b]);

    % Damping matrix
    d_m     = obj.d_m;
    d_g     = obj.d_g;
    d_b     = obj.d_b;
    d_mg    = obj.d_mg;
    d_gb    = obj.d_gb; % shorthands %#ok<*PROP>
    D = [   d_m + d_mg,     -d_mg,                  0; ...
            -d_mg,          d_g + d_mg + d_gb,      -d_gb; ...
            0,              -d_gb,                  d_b + d_gb  ];

    % Stiffness matrix
    k_g = obj.k_g;
    k_b = obj.k_b; % shorthands %#ok<*PROP>
    K = [   k_g,            -k_g,           0; ...
            -k_g,           k_g + k_b,      -k_b; ...
            0,              -k_b,           k_b         ];

    % State-space matrices
    A = [   zeros(size(I)),     eye(size(I)); ...
            -I\K,               -I\D            ];

    % Input
    k_t = obj.k_t;
    n   = obj.n;
    B = [ 0, 0, 0, k_t*n/I(1,1),  0, 0; ...
            0, 0, 0, 0,             0, 1/I(3,3) ]';
    
    % Output
    C = [   eye(size(A,2)); ...         % All states
            0, k_b, -k_b, 0, 0, 0   ];  % Torsion bar torque

end
