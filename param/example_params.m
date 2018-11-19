%% Motor parameters
% Save these lines in example_params.m
% for later use in other examples.
%
% Motor rotor inertia [kg m^2]
params.('I_m') = 9.407000e-02;     
% Motor Damping [Nms/rad]
params.('d_m') = 2.188643e-03;
% Torque constant [Nm/A]
params.('k_t') = 4.100000e-02;

%% Gear parameters
% Gear transmission ratio [.]
params.('n')   = 100;

%% Sensor parameters
% Sensor inertia [kg m^2]
params.('I_l') = 1.137000e-04;
% Sensor stiffness [Nm/rad]
params.('k_b') = 21000;     
