%OUTPUT_FIXED_RIGID_GEARBOX Get linear dynamics matrices - output link
% fixed, gearbox rigid
%
% [A, B, C, D, I, R, K] = jointObj.output_fixed_rigid_gearbox
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
%  This function is identical to full_dyn, but with the difference that the
%  joint output is now considered to be fixed and the gearbox is considered
%  to be rigid. This leads to a reduced model structure by an order of two.
%  Instead of environmental torque acting on the load, the load motion
%  becomes a velocity input.
%
% Examples::
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also full_dyn, output_fixed.

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

function [A, B, C, D, I, R, K] = output_fixed_rigid_gearbox(obj)
    
    % The computations below assume a state vector definition according to:
    % x = [q_g, q_l, q_g_dot]', where
    % q_g is the gearbox output angle
    % q_l is the load angle
    %
    % The '_dot' denotes the temporal derivative.

    % Inertia matrix
    I = obj.I_m + obj.I_g;

    % Damping matrix
    d_m     = obj.d_m;
    d_g     = obj.d_g;
    %d_mg    = obj.d_mg;
    d_gl    = obj.d_gl; % shorthands %#ok<*PROP>
    R = d_m + d_g + d_gl;

    % Stiffness matrix
    k_b = obj.k_b;
    K   = k_b;

    % State-space matrices
    A = [   zeros(size(I)),     zeros(size(I,1),1),     eye(size(I));
            0,                  0,                      0;
            -I\K,               I\K,                    -I\R            ];
        
    % Input
    % u = [tau_m, q_l_dot]
    B   = [ 0,      0;
            0,      1;
            1/I,    d_gl/I  ];
    
    % Output
    C = [   1,      0,      0;          % motor position
            1,      0,      0;          % gear position
            0,      1,      0;          % link position
            0,      0,      1;          % motor velocity
            0,      0,      1;          % gear velocity
            0,      0,      0;          % link velocity
            k_b,    -k_b,   d_gl    ];  % Torsion bar torque

    % Direct Feedthrough
    nIn     = size(B,2);
    nOut    = size(C,1);
    D       = zeros(nOut,nIn);
    if obj.isSym
        D = sym(D);
    end
    D(6,2)  = 1;        % Link velocity
    D(7,2)  = -d_gl;    % Torque
    
end

