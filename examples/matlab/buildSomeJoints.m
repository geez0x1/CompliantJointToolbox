% BUILDSOMEJOINTS Example script that builds some example joints
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
% #! Build Some Joints

% clear all;
close all;

% Create jointBuilder
jb = jointBuilder;
jb.overwrite = 1;

%% Build various types of joint models

% Parameter names to use
allParams = {'CE_Big';
    'WMBig10k_ds';
    'WMBig10k';
    'WMBig2300_ds'};
nPar = numel(allParams);

disp('Creating various joint types..');

% For each parameter file..
for iPar = 1:nPar
    % No friction
    jb.buildJoint(allParams{iPar}, 'full_dyn_no_friction');
    jb.buildJoint(allParams{iPar}, 'output_fixed_no_friction');
    jb.buildJoint(allParams{iPar}, 'rigid_gearbox_no_friction');
    jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox_no_friction');
    jb.buildJoint(allParams{iPar}, 'rigid_no_friction');
    
    % No nonlinearities
    jb.buildJoint(allParams{iPar}, 'full_dyn');
    jb.buildJoint(allParams{iPar}, 'output_fixed');
    jb.buildJoint(allParams{iPar}, 'rigid_gearbox');
    jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox');
    
    % Coulomb friction - symmetric
    jb.buildJoint(allParams{iPar}, 'full_dyn', 'coulomb');
    jb.buildJoint(allParams{iPar}, 'output_fixed', 'coulomb');
    jb.buildJoint(allParams{iPar}, 'rigid_gearbox', 'coulomb');
    jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox', 'coulomb');
    
    % Coulomb friction - asymmetric
    jb.buildJoint(allParams{iPar}, 'full_dyn', 'coulomb_asym');
    jb.buildJoint(allParams{iPar}, 'output_fixed', 'coulomb_asym');
    jb.buildJoint(allParams{iPar}, 'rigid_gearbox', 'coulomb_asym');
    jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox', 'coulomb_asym');
    
    % Multiple nonlinear terms
    jb.buildJoint(allParams{iPar}, 'full_dyn', {'coulomb', 'viscous_asym'});
    jb.buildJoint(allParams{iPar}, 'output_fixed', {'coulomb_asym', 'viscous_asym'});
    jb.buildJoint(allParams{iPar}, 'rigid_gearbox', {'coulomb', 'viscous_asym'});
    jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox', {'coulomb_asym', 'viscous_asym'});
end

disp('Completed.');