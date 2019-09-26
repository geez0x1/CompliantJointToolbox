% SETCJTPATH Script that adds the Compliant Joint Toolbox (CJT) to your Matlab
% search path.
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

fun_path = mfilename('fullpath');
tbx_path = fileparts(fun_path);

% Add files in root
addpath(tbx_path);

% Create build directory if it doesn't exist
if ~exist([tbx_path, filesep, 'build'], 'dir')
    mkdir([tbx_path, filesep, 'build']);
end

% Add paths recursively
addpath(genpath([tbx_path, filesep, 'build']));
addpath(genpath([tbx_path, filesep, 'examples']));
addpath(genpath([tbx_path, filesep, 'lib']));
addpath(genpath([tbx_path, filesep, 'lib', filesep,'include']));
addpath(genpath([tbx_path, filesep, 'model']));
addpath(genpath([tbx_path, filesep, 'param']));
addpath(genpath([tbx_path, filesep, 'templates']));
addpath(genpath([tbx_path, filesep, 'tools']));

% Status
disp('Added Compliant Joint Toolbox (CJT) paths.');