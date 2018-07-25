% MASK_KALMAN_OBSERVER Mask code for the Kalman state observer block
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

% Get linear dynamics matrices
[A, B, C, D, I, R, K] = jointObj.getDynamicsMatrices();

% Build vector of input indices
inputIdx = [];
if (inputIdx1_enable)
    inputIdx = [inputIdx(:) 1];
end
if (inputIdx2_enable)
    inputIdx = [inputIdx(:) 2];
end

% Get Kalman observer
% inputIdx may be modified by getKalman() to add the secondary input to
% enforce observability (see documentation for getKalman()).
[kest, L, Cc, inputIdx] = getKalman(jointObj, inputIdx, outputIdx, var_u, var_y);

% Discretize using Tustin transform
kest_tust               = c2d(kest, jointObj.Ts, 'tustin');

% Get demuxing indices
demux_N                 = size(kest.C,1);               % Kalman state-space block output size
demux_yIdx              = 1:length(outputIdx);        	% Estimated outputs indices
demux_stateIdx          = length(outputIdx)+1:demux_N;	% Estimated state indices

% Output
%disp('Kalman gain L:');
%disp(L);
%disp('Estimator kest:')
%kest
%disp('Original dynamic system sys:');
%jointObj.getStateSpace()
