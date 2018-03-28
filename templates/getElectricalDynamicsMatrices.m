%%GETELECTRICALDYNAMICSMATRICES Template for obtaining state space matrices for the linear electrical dynamics.
%%
%% [A, B, C, D] = jointObj.getElectricalDynamicsMatrices
%%
%% jointObj is the instance of the joint class object for which this
%% function has been called.
%%
%% Outputs::
%%
%%
%%   A: System matrix
%%   B: Input matrix
%%   C: Output matrix
%%   D: Direct Feedthrough matrix
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
function [A, B, C, D] = getElectricalDynamicsMatrices(obj)
    %s
end
