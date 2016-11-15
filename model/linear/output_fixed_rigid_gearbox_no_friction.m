% OUTPUT_FIXED_RIGID_GEARBOX_NO_FRICTION  Get linear dynamics matrices - output 
% link fixed, gearbox rigid, no friction
%
% [A, B, C, I, D, K] = jointObj.output_fixed_rigid_gearbox_no_friction
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
%  This function is identical to output_fixed_rigid_gearbox, but with the
%  difference that the bearing friction coefficients are set to zero.

%
% Examples::
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also full_dyn, output_fixed, output_fixed_rigid_gearbox.

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

function [A, B, C, I, D, K] = output_fixed_rigid_gearbox_no_friction(obj)
    
    % The computations below assume a state vector definition according to:
    % x = [q_g, q_g_dot,]', where 
    % q_g is the gearbox output angle
    %
    % The '_dot' denotes the temporal derivative.
    
    % Inertia matrix
    I = obj.I_m + obj.I_g;

    % Damping matrix
    D = 0 * (obj.d_m + obj.d_g + obj.d_gl);

    % Stiffness matrix
    k_b = obj.k_b;
    K   = k_b;

    % State-space matrices
    A = [   zeros(size(I)),     eye(size(I)); ...
            -I\K,               -I\D            ];
        
    % Input
    k_t = obj.k_t;
    n   = obj.n;
    B   = [ 0,              0, ...
            k_t*n/I(1,1),   0       ];
    
    % Output
    C = [   1,      0; ... % motor position
            1,      0; ... % gear position
            0,      0; ... % link position
            0,      1; ... % motor velocity
            0,      1; ... % gear velocity
            0,      0; ... % link velocity
            k_b,    0   ]; % Torsion bar torque

end

