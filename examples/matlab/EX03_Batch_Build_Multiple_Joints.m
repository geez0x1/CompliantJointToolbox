% BATCH_BUILD_MULTIPLE_JOINTS Example for quick batch creations of joint 
% models in different parameter / dynamics configurations. 
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
% #! Batch Build Some Joints

%% Create jointBuilder
 % Instantiate jointBuilder
jb = jointBuilder;
 % Configure jointBuilder: overwrite existing files without asking.
jb.overwrite = 1;

%% Configure parameters 
% Parameter sets to use for the joint models
allParams = {'cjt_Pomegranate_160_9000';
    'cjt_Orange_100_6000';
    'cjt_Lemon_50_6000';
    'cjt_Avocado_100_21000';
    'cjt_Lime_50_21000'};
nPar = numel(allParams);

%% Run joint building

% For each parameter file..
for iPar = 1:nPar
    % * No friction
    jb.buildJoint(allParams{iPar}, 'full_dyn_no_friction');
    jb.buildJoint(allParams{iPar}, 'output_fixed_no_friction');
    jb.buildJoint(allParams{iPar}, 'rigid_gearbox_no_friction');
    jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox_no_friction');
    jb.buildJoint(allParams{iPar}, 'rigid_no_friction');
    
    % * No nonlinearities
    jb.buildJoint(allParams{iPar}, 'full_dyn');
    jb.buildJoint(allParams{iPar}, 'output_fixed');
    jb.buildJoint(allParams{iPar}, 'rigid_gearbox');
    jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox');
    
    % * Coulomb friction - symmetric
    jb.buildJoint(allParams{iPar}, 'full_dyn', 'coulomb');
    jb.buildJoint(allParams{iPar}, 'output_fixed', 'coulomb');
    jb.buildJoint(allParams{iPar}, 'rigid_gearbox', 'coulomb');
    jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox', 'coulomb');
    
    % * Coulomb friction - asymmetric
    jb.buildJoint(allParams{iPar}, 'full_dyn', 'coulomb_asym');
    jb.buildJoint(allParams{iPar}, 'output_fixed', 'coulomb_asym');
    jb.buildJoint(allParams{iPar}, 'rigid_gearbox', 'coulomb_asym');
    jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox', 'coulomb_asym');
    
    % * Multiple nonlinear terms
    jb.buildJoint(allParams{iPar}, 'full_dyn', {'coulomb', 'viscous_asym'});
    jb.buildJoint(allParams{iPar}, 'output_fixed', {'coulomb_asym', 'viscous_asym'});
    jb.buildJoint(allParams{iPar}, 'rigid_gearbox', {'coulomb', 'viscous_asym'});
    jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox', {'coulomb_asym', 'viscous_asym'});
end

%% Finish up
% Look how many nice joint models we have! You may wish to inspect and play with the generated models. Otherwise, 
% if you wish to remove all the files that have been created by the builder, call the jointBuilder purge method:
%
%   |jb.purge|
%