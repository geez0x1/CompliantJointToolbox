% WALK-MAN Big Motor with 10000 Nm/rad torsion bar
n = 100;
% Inertiae
params.('I_m')      = 1.38 * 1.0059;       %% Motor rotor inertia [kg m^2]
params.('I_g')      = 1.38 * 0.67057;      %% Gear inertia [kg m^2]
params.('I_b')      = 0.0867;       %% Torsion bar inertia [kg m^2]
% Stiffnesses
params.('k_g')      = 75e3;        	%% Gearbox stiffness [Nm/rad]
params.('k_b')      = 13678;     	%% Torsion bar stiffness [Nm/rad]
% Linear viscous friction
params.('d_m')      = 3.59;     	%% Motor Damping [Nms/rad]
params.('d_g')      = 3.59;      	%% Gearbox damping [Nms/rad]
params.('d_b')      = 3.59;       	%% Torsion bar damping [Nms/rad]
% Asymmetric viscous friction
params.('d_m_n')    = 3.40;       	%% Motor Damping - negative direction [Nms/rad]
params.('d_g_n')    = 3.40;       	%% Gearbox Damping - negative direction [Nms/rad]
params.('d_b_n')    = 3.40;       	%% Torsion bar damping - negative direction [Nms/rad]
% Linear internal viscous friction
params.('d_mg')     = 296.0;       	%% Gearbox internal damping [Nms/rad]
params.('d_gb')     = 35.0;       	%% Torsion bar internal damping [Nms/rad]
% Coulomb friction
params.('d_cm')     = 11.347;     	%% Motor Coulomb damping [Nm]
params.('d_cg')     = 11.347;     	%% Gearbox Coulomb damping [Nm]
params.('d_cb')     = 0.0;        	%% Torsion bar Coulomb damping [Nm]
% Asymmetric Coulomb friction
params.('d_cm_n')   = 11.667;     	%% Motor Coulomb damping - negative direction [Nm]
params.('d_cg_n')   = 11.667;     	%% Gearbox Coulomb damping - negative direction [Nm]
params.('d_cb_n')   = 0.0;          %% Torsion bar Coulomb damping - negative direction [Nm]
% Misc
params.('n')        = n;          	%% Gear ratio []
params.('k_t')      = 0.0453;    	%% Torque constant [Nm/A]
params.('r')        = 0.099;      	%% Armature resistance [Ohm]
params.('x')        = 0.000168;     %% Armature inductance [H]
params.('Ts')       = 1E-3;       	%% Sampling time [s]
params.('v_0')     	= 24;			%% Operating voltage [V]
params.('i_c')      = 40;			%% Max. continuous current [A]
params.('i_p')		= 40; 			%% Peak current [A]
params.('dq_c')     = 5.86;			%% Max. continuous speed (output) [rad/s]
params.('dq_p')		= 5.86; 		%% Max. peak speed (output) [rad/s]