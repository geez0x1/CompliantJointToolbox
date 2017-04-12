%WMBIG2300_DS Parameter script for the WALK-MAN Big Motor with > 2000 Nm/rad torsion bar
%
% Notes::
%  All inertiae/damping is reflected to link side using n^2
%  Parameters are based on values obtained from datasheets (ds)
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also genericJoint, full_dyn, WMBig10k_ds.

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

% Inertiae
params.('m')         = 2;                    %% Total mass [kg]
params.('Im_m')      = 5.480E-5;             %% Motor rotor inertia [kg m^2]
params.('Im_g')      = 2.63E-5;              %% Gear inertia [kg m^2]
params.('I_l')       = 0.0867E-4;            %% Load inertia [kg m^2]
% Stiffnesses
params.('k_g')      = 31E3;                 %% Gearbox stiffness [Nm/rad]
params.('k_b')      = 2300;                 %% Torsion bar stiffness [Nm/rad]
% Linear viscous friction
params.('dm_m')     = 4.178E-05;            %% Motor Damping [Nms/rad]
params.('dm_g')     = 2E-4;                 %% Gearbox damping [Nms/rad]
params.('d_l')      = 1E-5;                 %% Load damping [Nms/rad]
% Asymmetric viscous friction
params.('dm_m_n')   = 4.178E-05;            %% Motor Damping - negative direction [Nms/rad]
params.('dm_g_n')   = 2E-4;                 %% Gearbox Damping - negative direction [Nms/rad]
params.('d_l_n')    = 1E-5;                 %% Load damping - negative direction [Nms/rad]
% Linear internal viscous friction
params.('dm_mg')    = 0.5E-4;               %% Gearbox internal damping [Nms/rad]
params.('d_gl')     = 0.5;                  %% Torsion bar internal damping [Nms/rad]
% Coulomb friction
params.('dm_cm')    = 1E-5;                 %% Motor Coulomb damping [Nm]
params.('dm_cg')    = 1E-5;                 %% Gearbox Coulomb damping [Nm]
params.('d_cl')     = 0.02;                 %% Load Coulomb damping [Nm]
% Asymmetric Coulomb friction
params.('dm_cm_n')  = 5E-6;                 %% Motor Coulomb damping - negative direction [Nm]
params.('dm_cg_n')  = 5E-6;                 %% Gearbox Coulomb damping - negative direction [Nm]
params.('d_cl_n')   = 0;                    %% Load Coulomb damping - negative direction [Nm]
% Stiction
params.('dm_s')     = 8.9E-4;               %% Break away torque [Nm]
params.('vm_s')     = 1;                    %% Stribeck velocity range [rad/s]
% Torque ripple sources
params.('rip_types')= [1, 2];               %% Torque ripple types (see torque_ripple())
params.('ripm_a1')  = [15E-3, 0.001];       %% Cosine amplitudes ([Nm] and/or [], see torque_ripple()) (second param to be updated!)
params.('ripm_a2')  = [0, 0];               %% Sine amplitudes [Nm] ([Nm] and/or [], see torque_ripple())
params.('ripm_f')   = [6, 2];               %% Spatial frequencies [periods/revolution]
% Misc
params.('n')        = 100;                  %% Transmission ratio []
params.('k_t')      = 0.0453;               %% Torque constant [Nm/A]
params.('r')        = 0.0885;               %% Armature resistance [Ohm]
params.('x')        = 0.000140;             %% Armature inductance [H]
params.('Ts')       = 1e-3;                 %% Sampling time [s]
% Operating/max conditions
params.('v_0')      = 24;                   %% Operating voltage [V]
params.('i_p')      = 40;                   %% Peak current [A]
params.('dq_p')     = 5.86;                 %% Max. peak speed (output) [rad/s]
