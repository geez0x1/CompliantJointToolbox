% MYFIRSTJOINTS Example script that builds a small number of example joints
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
% #! My First Joints

%% Instantiate a jointBuilder
jb = jointBuilder;

%% Build joint model classes
disp('-------------------------------------------------');
disp('My First Joints: Building a few example joints.');
disp('-------------------------------------------------');
disp('1. my_linear_joint:');
disp('   A linear model of a WALKMAN leg actuator with 10 kNm/rad torque sensor stiffness.');
disp('   Command: jb.buildJoint(''WMBig10k'', [], [], [],''my_linear_joint'');');
jb.buildJoint('WMBig10k',... parameters
    [],... full linear dynamics
    [],... no nonlinear dynamics
    [],... static electric subsystem
    'my_linear_joint'); % class name
disp('-------------------------------------------------');
disp('2. My_joint_with_Coulomb_friction:');
disp('   The same actuator model with rigid gearbox and asymmetric Coulomb and viscous friction as nonlinear dynamics term.');
disp('   Command: jb.buildJoint(''WMBig10k'', ''rigid_gearbox'', {''coulomb_asym'', ''viscous_asym''}, [] ''my_nonlinear_joint'');');
jb.buildJoint('WMBig10k',...parameters
    'rigid_gearbox',...  linear dynamics
    {'coulomb_asym',...  nonlinear 
     'viscous_asym'},... dynamics
    [],...     static electric subsystem
    'my_nonlinear_joint'); % class name
disp('-------------------------------------------------');
disp('3. My_complex_joint:');
disp('   Next, the nonlinear rigid gearbox model uses a locked output and electric dynamics including winding inductance.');
disp('   Command: jb.buildJoint(''WMBig10k'', ''output_fixed_rigid_gearbox'', {''coulomb_asym'', ''viscous_asym''}, ''electric_dyn'', ''my_locked_joint'');');
disp('-------------------------------------------------');
% In addition, output is locked output 
%and a full linear electric dynamics
jb.buildJoint('WMBig10k', ... parameters
    'output_fixed_rigid_gearbox', ... linear dynamics
    {'coulomb_asym',...  nonlinear
    'viscous_asym'}, ... dynamics
    'electric_dyn',... electro-dynamics
    'my_locked_joint');  % class name

%% Instantiate joint models
% add build directory to search path
addpath(jb.buildDir) 
% declare joint objects
jointA = my_linear_joint;
jointB = my_nonlinear_joint;
jointC = my_locked_joint;

% Finish up
disp('Completed.');