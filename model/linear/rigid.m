%RIGID Get linear dynamics matrices - fully rigid joint
%
% [A, B, C, D, I, R, K] = jointObj.rigid
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

function [A, B, C, D, I, R, K] = rigid(obj)
    
    % The computations below assume a state vector definition according to:
    % x = [q_l, q_l_dot]', where 
    % q_l is the flange angle (output of the torsion bar)
    %
    % The '_dot' denotes the temporal derivative.
    
    % Inertia matrix
    I = obj.I_m + obj.I_g + obj.I_l;

    % Damping matrix
    d_m = obj.d_m;
    d_g = obj.d_g;
    d_l = obj.d_l; % shorthands %#ok<*PROP>
    R = d_m + d_g + d_l;

    % There is no stiffness (rigid joint)
    K = 0;

    % State-space matrices
    A = [   zeros(size(I)),     eye(size(I));
            -I\K,               -I\R            ];

    % Input
    % u = [tau_m, tau_e]
    k_t = obj.k_t;
    n   = obj.n;
    B   = [ 0,              0;
            k_t*n/I(1,1),	1/I(1,1)         ];
    
    % Output
    C = [   1,  0;      % motor position
            1,  0;      % gear position
            1,  0;      % link position
            0,  1;      % motor velocity
            0,  1;      % gear velocity
            0,  1;      % link velocity
            0,  0   ];  % Torsion bar torque

    % Direct Feedthrough
    D = [0, 0];
end
