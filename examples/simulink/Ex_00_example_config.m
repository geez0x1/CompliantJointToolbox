% EX_00 An minimum example joint class generation to use the joint model blocks.
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
% #! 00: Example Joint Configuration

dispText = ...
{'-----'
 'EX 00'
 '-----'
 'This is an example of a minimal joint class generation for the usage of CJT Simulink blocks. This file is called in'
 'all other example files to configure the example joint.'
 ' ' 
};

% displayFormattedText(dispText)

%% Configure Joint Model
jb = jointBuilder;
jb.overwrite = 1;

jb.buildJoint(  'cjt_Orange_80_6000', ...
                'output_fixed_rigid_gearbox', ...
                'coulomb_asym', ...
                [], ...
                'example_joint'                 );

addpath(jb.buildDir);

jointObj = example_joint;

%% Configure Chirp Source
f_min = 0.01;   % Initial frequency in Hz
f_max = 100;    % Final frequency in Hz
t_max = f_max;  % Simulation stop time in s


