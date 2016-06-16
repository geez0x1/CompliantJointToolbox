%%GETSTATESPACE Template for obtaining a state space model of the linear
%% dynamics of a joint.
%%
%% sys = jointObj.getStateSpace
%%
%% jointObj is the instance of the joint class object for which this
%% function has been called.
%%
%% Outputs::
%%
%%   sys: State space model of the compliant joint
%%
%% Notes::
%%
%%
%% Examples::
%%
%% Author::
%%  Joern Malzahn
%%  Wesley Roozing
%%
%% See also genericJoint, full_dyn, getNonlinearDynamics.

%% Copyright (C) 2016, by Joern Malzahn, Wesley Roozing
%%
%% This file is part of the Compliant Joint Toolbox (CJT).
%%
%% CJT is free software: you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation, either version 3 of the License, or
%% (at your option) any later version.
%%
%% CJT is distributed in the hope that it will be useful, but WITHOUT ANY
%% WARRANTY; without even the implied warranty of MERCHANTABILITY or
%% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
%% License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with CJT. If not, see <http://www.gnu.org/licenses/>.
%%
%% For more information on the toolbox and contact to the authors visit
%% <https://github.com/geez0x1/CompliantJointToolbox>
function sys = getStateSpace(obj)
    [A, B, C, ~, ~, ~] = obj.getDynamicsMatrices();
    D = 0;
    sys = ss(A, B, C, D);
end
