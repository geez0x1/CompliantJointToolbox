% BUILD MULTIPLE JOINTS Example script that builds a small number of simple and also more complex example joint models.
% The example is based on the parameter sets of the *TREE Actuators* <www.treerobotics.eu>.
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
% #! Build Multiple Joints

%% Instantiate a jointBuilder
jb = jointBuilder;

%% Build linear joint model classes
% First we build a set of linear models for different joint types from the TREE family <www.treerobotics.eu>. The models
% will feature the full linear model dynamics captured in the Compliant Joint Toolbox. This is a three mass system
% comprising the motor, gearbox and torque sensor inertia. Both, the gearbox and the sensor are elastic elements. The
% actuator output can move and is subject to load and motion disturbances.
%
% * The first joint class will implement linear model of a Lime type actuator with a gearbox ratio of 50:1 and a 
% stiffness of 21 kNm/rad.
% 
jb.buildJoint('cjt_Lime_50_21000',... parameters
    [],... full linear dynamics
    [],... no nonlinear dynamics
    [],... static electric subsystem
    'my_linear_lime'); % class name
% 
% * The second joint class will be an ordinary Orange type actuator with a gearbox ratio of 100:1 and a stiffness of 
% 6 kNm/rad.
%
jb.buildJoint('cjt_Orange_100_6000',... parameters
    [],... full linear dynamics
    [],... no nonlinear dynamics
    [],... static electric subsystem
    'my_ordinary_orange'); % class name

%% Build nonlinear and more complex joint model classes
% Now we build model classes with more complex joint dynamics. The actuators will be a Lemon and an Avocado type 
% actuator. 
% 
% * The third actuator model is a more advanced model of an Avocado type actuator, which features high rigidity with 
% asymmetric Coulomb and viscous friction as nonlinear dynamics terms.
%
jb.buildJoint('cjt_Avocado_100_21000',...parameters
    'rigid',...  linear dynamics
    {'coulomb_asym',...  nonlinear 
     'viscous_asym'},... dynamics
    [],...     static electric subsystem
    'my_advanced_avocado'); % class name
% 
% * The fourth actuator model is the Lemon type actuator, which also features asymmetric Coulomb and 
% viscous friction as nonlinear dynamics terms but now we incorporate also electric dynamics in the model and lock the 
% actuator output. Models for a locked actuator output simulate very large load impedances and are often used for torque
% controller evaluation.
%
jb.buildJoint('cjt_Lemon_50_6000', ... parameters
    'output_fixed_rigid_gearbox', ... linear dynamics
    {'coulomb_asym',...  nonlinear
    'viscous_asym'}, ... dynamics
    'electric_dyn',... electro-dynamics
    'my_locked_lemon');  % class name

%% Instantiate joint models
% After creating the joint model classes we can use them. First we need to add the build directory to the search path.
%
addpath(jb.buildDir) 

%
% Then we instantiate one actuator of each type
%
jointA = my_linear_lime;
jointB = my_ordinary_orange;
jointC = my_advanced_avocado;
jointD = my_locked_lemon;

%% Finish up
% If you wish to remove all the files that have been created by the builder, call the jointBuilder purge method:
%
%   |jb.purge|
%