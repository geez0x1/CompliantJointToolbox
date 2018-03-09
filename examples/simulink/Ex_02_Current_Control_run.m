% EX_02 An example of a PI plus feedforward current controller for a compliant actuator with locked output.
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
 'EX 02'
 '-----'
 'This example demonstrates the current control of a compliant actuator with locked '
 'output based on a PI controller with voltage feedforward action.'
 ' '
 'The model used in this example considers the electrical dynamics.' 
};

displayFormattedText(dispText)

%% Configure example environment

% Joint Builder
jb = jointBuilder;
jb.overwrite = 1;

jb.buildJoint('cjt_Orange_80_6000',...
       'output_fixed','coulomb_asym',...
       'electric_dyn',...
       'example_joint_locked');

addpath(jb.buildDir);

jObj = example_joint_locked;
jObj.Ts = 0.5e-4; % The current controller runs at 20 kHz

%% Configure Chirp Source
f_min = 0.01; % Initial frequency in Hz
f_max = 1000;  % Final frequency in Hz
t_max = 3; % Simulation stop time in s


mdlName = 'Ex_02_Current_Control.mdl';

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
hold on;
plot(t,ref,'r','linewidth',1.0,'DisplayName','reference');
plot(t,res,'k','linewidth',0.5,'DisplayName','response');
legend show
xlabel('t [s]')
ylabel('Current [A]')
title('Current chirp response')