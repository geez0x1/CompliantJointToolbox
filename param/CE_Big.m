% Centauro Big Motor with 2300 Nm/rad torsion bar
% All inertiae/damping is reflected to link side using n^2
n = 100;
% Inertiae
params.('I_m')      = 3.04E-005 * n^2;    %% Motor rotor inertia [kg m^2]
params.('I_g')      = 1.33 * 6.6944E-005 * n^2;    %% Gear inertia [kg m^2] % TO BE UPDATED!
params.('I_b')      = 8.6719E-006 * n^2;    %% Torsion bar inertia [kg m^2] % TO BE UPDATED!
% Stiffnesses
params.('k_g')      = 31e3;             	%% Gearbox stiffness [Nm/rad]
params.('k_b')      = 31e3;                 %% Torsion bar stiffness [Nm/rad] % TO BE UPDATED! The strain based torque sensor is super stiff.
% Linear viscous friction
params.('d_m')      = 0.0064;               %% Motor Damping [Nms/rad]
params.('d_g')      = 6.0;                  %% Gearbox damping [Nms/rad]   % TO BE UPDATED!
params.('d_b')      = 7.0;                  %% Torsion bar damping [Nms/rad] % TO BE UPDATED!
% Asymmetric viscous friction
params.('d_m_n')    = 0.0064;               %% Motor Damping - negative direction [Nms/rad]
params.('d_g_n')    = 5.0;                  %% Gearbox Damping - negative direction [Nms/rad] % TO BE UPDATED!
params.('d_b_n')    = 5.0;                  %% Torsion bar damping - negative direction [Nms/rad] % TO BE UPDATED!
% Linear internal viscous friction
params.('d_mg')     = 144.0;                %% Gearbox internal damping [Nms/rad]     % TO BE UPDATED!
params.('d_gb')     = 2.0;                  %% Torsion bar internal damping [Nms/rad] % TO BE UPDATED!
% Coulomb friction
params.('d_cm')     = 3.2;                 	%% Motor Coulomb damping [Nm]  
params.('d_cg')     = 3.5;            	    %% Gearbox Coulomb damping [Nm] % TO BE UPDATED!
params.('d_cb')     = 0.0;                  %% Torsion bar Coulomb damping [Nm] % TO BE UPDATED!
% Asymmetric Coulomb friction
params.('d_cm_n')   = 0.0;                  %% Motor Coulomb damping - negative direction [Nm]
params.('d_cg_n')   = 0.0;                  %% Gearbox Coulomb damping - negative direction [Nm]
params.('d_cb_n')   = 0.0;                  %% Torsion bar Coulomb damping - negative direction [Nm]
% Misc
params.('n')        = n;                    %% Gear ratio []
params.('k_t')      = 0.0455;               %% Torque constant [Nm/A]
params.('r')        = 0.099;                %% Armature resistance [Ohm]
params.('x')        = 0;                    %% Armature inductance [H]
params.('Ts')       = 1E-3;                 %% Sampling time [s]
params.('v_0')     	= 24;					%% Operating voltage [V]
params.('i_c')      = 20;					%% Max. continuous current [A]
params.('i_p')		= 20; 					%% Peak current [A]
params.('dq_c')     = 5.86;					%% Max. continuous speed (output) [rad/s]
params.('dq_p')		= 5.86; 		     	%% Max. peak speed (output) [rad/s]