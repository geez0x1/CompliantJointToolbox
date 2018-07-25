% MASK_LUENBERGER_OBSERVER Mask code for the Luenberger state observer block
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also mask_Kalman_observer.

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

% Get Luenberger observer
[sys_hat, L, Cc] = getObserver(jointObj, outputIdx, poles);

% Discretize using Tustin transform
sys_hat_tust = c2d(sys_hat, jointObj.Ts, 'tustin');

% Output
%disp('Observer gain L:');
%disp(L);
%disp('Estimator sys_hat:')
%sys_hat
%disp('Original dynamic system sys:');
%jointObj.getStateSpace()
%sys_hat.A - A
%sys_hat.B(:,1) - B