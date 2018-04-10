% EX_04 An example of a simple PD controller
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
% #! A Test Simulink Example

dispText = ...
{'-----'
 'EX 04'
 '-----'
 'This example demonstrates the torque control of a compliant actuator with locked '
 'output using a simple PD controller.'
 ' ' 
};

displayFormattedText(dispText)

%% Configure example environment
% Standard configuration
Ex_00_example_config
mdlName = 'Ex_05_PD_Control.mdl';

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