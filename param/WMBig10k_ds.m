% WALK-MAN Big Motor with 10 kNm/rad torsion bar
% All inertiae/damping is reflected to link side using n^2
% These parameters are based on values obtained from datasheets (ds)
n = 100;
% Inertiae
params.('I_m')      = 2.63E-005 * n^2;      %% Motor rotor inertia [kg m^2]
params.('I_g')      = 5.48E-005 * n^2;      %% Gear inertia [kg m^2]
params.('I_b')      = 0.1;                  %% Torsion bar inertia [kg m^2]
% Stiffnesses
params.('k_g')      = 31e3;                 %% Gearbox stiffness [Nm/rad]
params.('k_b')      = 10000;                %% Torsion bar stiffness [Nm/rad]
% Linear viscous friction
params.('d_m')      = 4.178e-05 * n^2;     	%% Motor Damping [Nms/rad]
params.('d_g')      = 2e-4 * n^2;         	%% Gearbox damping [Nms/rad]
params.('d_b')      = 1e-5;                 %% Torsion bar damping [Nms/rad]
% Asymmetric viscous friction
params.('d_m_n')    = 4.178e-05 * n^2;      %% Motor Damping - negative direction [Nms/rad]
params.('d_g_n')    = 2e-4 * n^2;        	%% Gearbox Damping - negative direction [Nms/rad]
params.('d_b_n')    = 1e-5;                 %% Torsion bar damping - negative direction [Nms/rad]
% Linear internal viscous friction
params.('d_mg')     = 0.5;                  %% Gearbox internal damping [Nms/rad]
params.('d_gb')     = 0.5;                  %% Torsion bar internal damping [Nms/rad]
% Coulomb friction
params.('d_cm')     = 1e-5 * n^2;          	%% Motor Coulomb damping [Nm]
params.('d_cg')     = 1e-5 * n^2;          	%% Gearbox Coulomb damping [Nm]
params.('d_cb')     = 0.02;                 %% Torsion bar Coulomb damping [Nm]
% Asymmetric Coulomb friction
params.('d_cm_n')   = 5e-6 * n^2;         	%% Motor Coulomb damping - negative direction [Nm]
params.('d_cg_n')   = 5e-6 * n^2;          	%% Gearbox Coulomb damping - negative direction [Nm]
params.('d_cb_n')   = 0;                    %% Torsion bar Coulomb damping - negative direction [Nm]
% Misc
params.('n')        = n;                    %% Gear ratio []
params.('k_t')      = 0.0453;    			%% Torque constant [Nm/A]
params.('r')        = 0.099;      			%% Armature resistance [Ohm]
params.('x')        = 0.000168;     		%% Armature inductance [H]
params.('Ts')       = 1E-3;                 %% Sampling time [s]
params.('v_0')     	= 24;					%% Operating voltage [V]
params.('i_c')      = 40;					%% Max. continuous current [A]
params.('i_p')		= 40; 					%% Peak current [A]
params.('dq_c')     = 5.86;					%% Max. continuous speed (output) [rad/s]
params.('dq_p')		= 5.86; 				%% Max. peak speed (output) [rad/s]