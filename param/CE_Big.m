%CE_BIG Parameter script for the Centauro Big Motor with strain gauge based torque sensor
%
%
% Notes::
%  All inertiae/damping is reflected to link side using n^2
%  The parameters in this script are partially obtained from experimental identification
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


% Inertiae
params.('m')        = 2;                    %% Total mass [kg]
params.('Im_m')     = 5.480E-5;             %% Motor rotor inertia [kg m^2] % TO BE UPDATED!
params.('Im_g')     = 2.63E-5;              %% Gear inertia [kg m^2] % TO BE UPDATED!
params.('I_l')      = 0.0867E-4;            %% Load inertia [kg m^2] % TO BE UPDATED!
% Stiffnesses
params.('k_g')      = 31E3;                 %% Gearbox stiffness [Nm/rad] % TO BE UPDATED!
params.('k_b')      = 31E3;                 %% Torsion bar stiffness [Nm/rad] % TO BE UPDATED! The strain based torque sensor is super stiff.
% Linear viscous friction
params.('dm_m')     = 0.0064E-4;            %% Motor Damping [Nms/rad] % TO BE UPDATED!
params.('dm_g')     = 6.0E-4;               %% Gearbox damping [Nms/rad] % TO BE UPDATED!
params.('d_l')      = 7.0;                  %% Load damping [Nms/rad] % TO BE UPDATED!
% Asymmetric viscous friction
params.('dm_m_n')   = 0.0064E-4;            %% Motor Damping - negative direction [Nms/rad] % TO BE UPDATED!
params.('dm_g_n')   = 5.0E-4;               %% Gearbox Damping - negative direction [Nms/rad] % TO BE UPDATED!
params.('d_l_n')    = 5.0;                  %% Load damping - negative direction [Nms/rad] % TO BE UPDATED!
% Linear internal viscous friction
params.('dm_mg')    = 144.0E-4;             %% Gearbox internal damping [Nms/rad] % TO BE UPDATED!
params.('d_gl')     = 2.0;                  %% Torsion bar internal damping [Nms/rad] % TO BE UPDATED!
% Coulomb friction
params.('dm_cm')    = 3.2E-4;               %% Motor Coulomb damping [Nm] % TO BE UPDATED!
params.('dm_cg')    = 3.5E-4;               %% Gearbox Coulomb damping [Nm] % TO BE UPDATED!
params.('d_cl')     = 0.0;                  %% Load Coulomb damping [Nm] % TO BE UPDATED!
% Asymmetric Coulomb friction
params.('dm_cm_n')  = 0.0;                  %% Motor Coulomb damping - negative direction [Nm] % TO BE UPDATED!
params.('dm_cg_n')  = 0.0;                  %% Gearbox Coulomb damping - negative direction [Nm] % TO BE UPDATED!
params.('d_cl_n')   = 0.0;                  %% Load Coulomb damping - negative direction [Nm] % TO BE UPDATED!
% Stiction
params.('dm_s')      = 8.9E-4;              %% Break away torque [Nm]
params.('vm_s')      = 1;                   %% Stribeck velocity range [rad/s]
% Torque ripple sources
params.('rip_types')= [1, 2];               %% Torque ripple types (see torque_ripple())
params.('ripm_a1')  = [15e-3, 0.001];       %% Cosine amplitudes ([Nm] and/or [], see torque_ripple()) (second param to be updated!)
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
