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

% Create jointBuilder
jb = jointBuilder;
%jb.purge();

% Info
disp('-------------------------------------------------');
disp('My First Joints: Building a few example joints.');
disp('-------------------------------------------------');
disp('1. My_simple_joint:');
disp('   A simple output-fixed joint with only linear dynamics.');
disp('   Command: jb.buildJoint(''WMBig10k'', ''output_fixed'', '''', ''My_simple_joint'');');
disp('-------------------------------------------------');
disp('2. My_joint_with_Coulomb_friction:');
disp('   The same joint with Coulomb friction as nonlinear dynamics term.');
disp('   Command: jb.buildJoint(''WMBig10k'', ''output_fixed'', ''coulomb'', ''My_joint_with_Coulomb_friction'');');
disp('-------------------------------------------------');
disp('3. My_complex_joint:');
disp('   A joint model with multiple nonlinear dynamics terms: asymmetric Coulomb and asymmetric viscous friction.');
disp('   Command: jb.buildJoint(''WMBig10k'', ''output_fixed'', {''coulomb_asym'', ''viscous_asym''}, ''My_complex_joint'');');
disp('-------------------------------------------------');

% Build a few example joints
jb.buildJoint('WMBig10k', 'output_fixed', '', 'My_simple_joint');
jb.buildJoint('WMBig10k', 'output_fixed', 'coulomb', 'My_joint_with_Coulomb_friction');
jb.buildJoint('WMBig10k', 'output_fixed', {'coulomb_asym', 'viscous_asym'}, 'My_complex_joint');

% Finish up
disp('Completed.');