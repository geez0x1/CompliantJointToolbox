% MASK_LINEAR_TRANSFER_FUNCTION_DOB Mask code for the LTF DOB block
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

% Make sure local values are defined
Pc      = 0;
Q_td    = 0;
PQ_td   = 0;

% If filename is given, assume the model and filters can be found in there
if (loadData)
    if (isempty(modelFname) || ~exist(modelFname, 'file'))
        error('Error: Model/filters data file not found!');
    else
        % Load file
        clear jointObj; % Clear mask parameter
        load(modelFname);
    end
else
    % Get linear DOB based on model
    [Pc, Q_td, PQ_td] = getLinearDOB(jointObj, 2*pi*f_c, measIdx, doPlot);
end

% Check if the necessary variables exist
if (~exist('jointObj', 'var') || ~exist('Q_td', 'var') || ~exist('PQ_td', 'var'))
    error('Linear Transfer Function DOB error: Not all variables exist');
end

% Get linear dynamics matrices
[A, B, C, D, I, R, K] = jointObj.getDynamicsMatrices();

% Discretize using Tustin transform
Q_td_tust   = c2d(Q_td, jointObj.Ts, 'tustin');
PQ_td_tust  = c2d(PQ_td, jointObj.Ts, 'tustin');
