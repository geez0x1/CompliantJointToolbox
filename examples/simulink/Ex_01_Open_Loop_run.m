% EX_01 An example of using a joint mechanical subsystem with direct open-loop torque input.
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

% The text in the following line defines the display name for this Example
% #! 01: Open-Loop Mechanical Subsystem

dispText = ...
{'-----'
 'EX 01'
 '-----'
 'This example demonstrates the basic usage of a joint mechanical subsystem block. A torque is directly fed to'
 'joint input and several measurements are displayed using scopes.'
 ' ' 
};

displayFormattedText(dispText)

%% Configure example environment
% Standard configuration
Ex_00_example_config
mdlName = 'Ex_01_Open_Loop.mdl';

%% Run model
% Open block diagram
open(mdlName)

% Run it
simOut = sim(mdlName,'SrcWorkspace','current');

%% Plot results
% Reshape data
t = simOut.tout;
ref = simOut.yout(:,1);
res = simOut.yout(:,2);

% Plot some results
figure(1)
clf;
plot(t,[ref,res])
legend({'reference', 'response'})
xlabel('t [s]')
ylabel('torque [Nm]')
title('Torque chirp response')