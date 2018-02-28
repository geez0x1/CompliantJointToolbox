% CREATE DATASHEETS The example shows how to instantiate and configure a datasheetGenerator class object to create 
% characteristic plots and datasheets for joint model classes. The example is based on the parameter sets of the 
% *TREE Actuators* <www.treerobotics.eu>.
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also cjtExamples, jointBuilder, genericJoint.

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
% #! A First Joint


%% Instantiate a jointBuilder
jb = jointBuilder;

%% Build joint model classes
% A model of a TREE Robotics Pomegranate actuator with 9 kNm/rad sensor stiffness Coulomb friction as well as
% electric dynamics including winding inductance and rigid gearbox.
jb.buildJoint('cjt_Orange_80_6000', ... parameters
    'rigid_gearbox', ... linear dynamics
    {'coulomb',...  {nonlinear
    'viscous_asym'}, ... dynamics}
    [],... electro-dynamics
    'My_Orange');  % class name

% add build directory to search path
addpath(jb.buildDir) 

%% Instantiate an object of the new class
aJoint = My_Orange;


%% Instantiate a dataSheetGenerator
dsGen = dataSheetGenerator(aJoint);

figure(1)
clf;
dsGen.draw_torque_speed_curve;


figure(2)
clf;
dsGen.draw_thermal_characteristics;