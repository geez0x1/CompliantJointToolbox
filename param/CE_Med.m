%CE_MED Parameter script for the Centauro Medium Motor with encoder based torque sensor
%
% Notes::
%  All inertiae/damping is reflected to link side using n^2
%  The parameters in this script are partially obtained from experimental identification
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also genericJoint, full_dyn.

% Copyright (C) 2017, by Joern Malzahn, Wesley Roozing
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

% Transmission ratio
n = 160;

% Inertiae
params.('I_m')      = 0.0 * n^2;            %% Motor rotor inertia [kg m^2] % TO BE UPDATED!
params.('I_g')      = 0.0 * n^2;            %% Gear inertia [kg m^2] % TO BE UPDATED!
params.('I_l')      = 0.0 * n^2;            %% Load inertia [kg m^2] % TO BE UPDATED!
% Stiffnesses
params.('k_g')      = 0.0;                  %% Gearbox stiffness [Nm/rad] % TO BE UPDATED!
params.('k_b')      = 0.0;                  %% Torsion bar stiffness [Nm/rad] % TO BE UPDATED! The strain based torque sensor is super stiff.
% Linear viscous friction
params.('d_m')      = 0.0;                  %% Motor Damping [Nms/rad] % TO BE UPDATED!
params.('d_g')      = 0.0;                  %% Gearbox damping [Nms/rad] % TO BE UPDATED!
params.('d_l')      = 0.0;                  %% Load damping [Nms/rad] % TO BE UPDATED!
% Asymmetric viscous friction - negative direction
params.('d_m_n')    = 0.0;                  %% Motor Damping - negative direction [Nms/rad] % TO BE UPDATED!
params.('d_g_n')    = 0.0;                  %% Gearbox Damping - negative direction [Nms/rad] % TO BE UPDATED!
params.('d_l_n')    = 0.0;                  %% Load damping - negative direction [Nms/rad] % TO BE UPDATED!
% Linear internal viscous friction
params.('d_mg')     = 0.0;                  %% Gearbox internal damping [Nms/rad] % TO BE UPDATED!
params.('d_gl')     = 0.0;                  %% Torsion bar internal damping [Nms/rad] % TO BE UPDATED!
% Coulomb friction
params.('d_cm')     = 0.0;                  %% Motor Coulomb damping [Nm] % TO BE UPDATED!
params.('d_cg')     = 0.0;                  %% Gearbox Coulomb damping [Nm] % TO BE UPDATED!
params.('d_cl')     = 0.0;                  %% Load Coulomb damping [Nm] % TO BE UPDATED!
% Asymmetric Coulomb friction - negative direction
params.('d_cm_n')   = 0.0;                  %% Motor Coulomb damping - negative direction [Nm] % TO BE UPDATED!
params.('d_cg_n')   = 0.0;                  %% Gearbox Coulomb damping - negative direction [Nm] % TO BE UPDATED!
params.('d_cl_n')   = 0.0;                  %% Load Coulomb damping - negative direction [Nm] % TO BE UPDATED!
% Stiction
params.('d_s')      = 0.0;                  %% Break away torque [Nm]
params.('v_s')      = 0.0;                  %% Stribeck velocity range [rad/s]
% Torque ripple sources
params.('rip_types')= [];                   %% Torque ripple types (see torque_ripple())
params.('rip_a1')   = [];                   %% Cosine amplitudes ([Nm] and/or [], see torque_ripple())
params.('rip_a2')   = [];                   %% Sine amplitudes [Nm] ([Nm] and/or [], see torque_ripple())
params.('rip_f')    = [];                   %% Spatial frequencies [periods/revolution]
% Misc
params.('n')        = n;                    %% Gear ratio []
params.('k_t')      = 0.0;                  %% Torque constant [Nm/A]
params.('r')        = 0.0;                  %% Armature resistance [Ohm]
params.('x')        = 0.0;                  %% Armature inductance [H]
params.('Ts')       = 1e-3;                 %% Sampling time [s]
% Operating/max conditions
params.('v_0')      = 48;                   %% Operating voltage [V]
params.('i_p')      = 0.0;                  %% Peak current [A]
params.('dq_p')     = 0.0;                  %% Max. peak speed (output) [rad/s]
