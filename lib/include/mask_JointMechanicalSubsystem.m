% MASK_JOINTMECHANICALSUBSYSTEM Mask code for JointMechanicalSubsystem
% model block
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also genericJoint.

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

% Get linear dynamics matrices
[A, B, C, D, I, R, K] = jointObj.getDynamicsMatrices();

% Build B matrix for nonlinear components
B_nonlinear = [zeros([size(B,1)-size(I,1),1]); 1./diag(I)];

% Get other properties
Ts      = jointObj.Ts;
Ts_elec = jointObj.Ts_elec;
k_b     = jointObj.k_b;
k_t     = jointObj.k_t;
n       = jointObj.n;
d_gl    = jointObj.d_gl;

% Handle MATLAB function block
handleMatlabFunctionBlock(gcb, jointObj);

% If input delay is enabled, set the flag to 1
input_delay_enable = 0;
if (input_delay > 0)
    input_delay_enable = 1;
end

% If output delay is enabled, set the flag to 1
output_delay_enable = 0;
if (output_delay > 0)
    output_delay_enable = 1;
end

% If noise is disabled, overwrite the given variances
if (input_noise_enabled == 0)
    var_u = 0;
end
if (output_noise_enabled == 0)
	var_y = 0;
end
