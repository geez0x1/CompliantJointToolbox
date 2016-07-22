% RIGID_NO_FRICTION Get linear dynamics matrices - fully rigid joint without
% friction
%
% [A, B, C, I, D, K] = jointObj.rigid
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
%  This function is identical to rigid model, but with the
%  difference that all friction coefficients are set to zero.
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

function [A, B, C, I, D, K] = rigid_no_friction(obj)
    
    % The computations below assume a state vector definition according to:
    % x = [q_b, q_b_dot]', where 
    % q_b is the flange angle (output of the torsion bar)
    %
    % The '_dot' denotes the temporal derivative.

    % Inertia matrix
    I = obj.I_m + obj.I_g + obj.I_b;

    % Damping matrix
    D = 0;

    % There is no stiffness (rigid joint)
    K = 0;

    % State-space matrices
    A = [   zeros(size(I)),     eye(size(I)); ...
            -I\K,               -I\D            ];

    % Input
    k_t = obj.k_t;
    n   = obj.n;
    B	= [ 0, k_t*n/I(1,1); ...
            0, 1/I(1,1)         ]';
    
    % Output
    C = [1, 0;  ... % motor position
         1, 0;  ... % gear position
         1, 0;  ... % link position
         0, 1;  ... % motor velocity
         0, 1;  ... % gear velocity
         0, 1;];... % link velocity
            
end
