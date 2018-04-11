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


