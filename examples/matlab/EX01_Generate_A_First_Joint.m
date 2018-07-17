% GENERATE A FIRST JOINT The example shows how to instantiate and configure a jointBuilder class object to create a 
% first joint model class. The example is based on the parameter sets of the *TREE Actuators* <www.treerobotics.eu>.
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
% #! 01: A First Joint


%% Instantiate a jointBuilder
jb = jointBuilder;

%% Build joint model classes
% A model of a TREE Robotics Pomegranate actuator with 9 kNm/rad sensor stiffness Coulomb friction as well as
% electric dynamics including winding inductance and rigid gearbox.
jb.buildJoint('cjt_Pomegranate_160_9000', ... parameters
    'rigid_gearbox', ... linear dynamics
    {'coulomb',...  {nonlinear
    'viscous_asym'}, ... dynamics}
    'electric_dyn',... electro-dynamics
    'My_Pomegranate');  % class name

% add build directory to search path
addpath(jb.buildDir) 

%% Instantiate an object of the new class
aJoint = My_Pomegranate


%% Next, we clean up
% The joint builder can delete all the files it has genereated through a single command. That helps you keep your hard
% disk clean.
jb.purge();

