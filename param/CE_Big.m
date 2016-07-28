%CE_BIG Parameter script for the Centauro Big Motor with 2300 Nm/rad torsion bar
%
%
% Notes::
%  All inertiae/damping is reflected to link side using n^2
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also genericJoint, full_dyn, WMBig10k.

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

n = 100;
% Inertiae
params.('I_m')      = 3.04E-005 * n^2;              %% Motor rotor inertia [kg m^2]
params.('I_g')      = 1.33 * 6.6944E-005 * n^2;     %% Gear inertia [kg m^2] % TO BE UPDATED!
params.('I_b')      = 8.6719E-006 * n^2;            %% Torsion bar inertia [kg m^2] % TO BE UPDATED!
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
% Cogging
params.('cog_a1')   = 1e-2 * n;             %% Cosine amplitude [Nm]
params.('cog_a2')   = 0 * n;                %% Sine amplitude [Nm]
params.('cog_f')    = 6;                    %% Spatial frequency [periods/revolution]
% Misc
params.('n')        = n;                    %% Gear ratio []
params.('k_t')      = 0.0455;               %% Torque constant [Nm/A]
params.('r')        = 0.099;                %% Armature resistance [Ohm]
params.('x')        = 0;                    %% Armature inductance [H]
params.('Ts')       = 1E-3;                 %% Sampling time [s]
% Operating/max conditions
params.('v_0')     	= 24;					%% Operating voltage [V]
params.('i_c')      = 20;					%% Max. continuous current [A]
params.('i_p')		= 20; 					%% Peak current [A]
params.('dq_c')     = 5.86;					%% Max. continuous speed (output) [rad/s]
params.('dq_p')		= 5.86; 		     	%% Max. peak speed (output) [rad/s]