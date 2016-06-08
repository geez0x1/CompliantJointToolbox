% setSEAPaths.m
% Include all the necessary files for using the SEA simulation library

% Add files in root
addpath(pwd);

% Create build directory if it doesn't exist
if ~exist([pwd, filesep, 'build'], 'dir')
    mkdir([pwd, filesep, 'build']);
end

% Add recursive paths
addpath(genpath([pwd, filesep, 'build']));
addpath(genpath([pwd, filesep, 'lib']));
addpath(genpath([pwd, filesep, 'model']));
addpath(genpath([pwd, filesep, 'param']));
addpath(genpath([pwd, filesep, 'templates']));
addpath(genpath([pwd, filesep, 'tools']));

% Status
disp('Added SEA simulation library paths.');