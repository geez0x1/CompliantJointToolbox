%B11 Parameter script for the WALK-MAN Big Motor B11 with 12823 Nm/rad torsion bar
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

params.('n')        = 100;                      %% Transmission ratio []
% Inertiae
params.('I_m')      = 0.55e-4 * params.('n')^2; %% Motor rotor inertia [kg m^2]
params.('I_g')      = 0.263e-4 * params.('n')^2;%% Gear inertia [kg m^2]
params.('I_l')      = 1.137e-4;                 %% Torsion bar inertia [kg m^2]
% Stiffnesses
params.('k_g')      = 31e3;                     %% Gearbox stiffness [Nm/rad]
params.('k_b')      = 12823;                    %% Torsion bar stiffness [Nm/rad]
% Linear viscous friction
params.('d_m')      = 0;                        %% Motor Damping [Nms/rad]
params.('d_g')      = 12.1997;                  %% Gearbox damping [Nms/rad]
params.('d_l')      = 0;                        %% Torsion bar damping [Nms/rad]
% Asymmetric viscous friction
params.('d_m_n')    = 0;                        %% Motor Damping - negative direction [Nms/rad]
params.('d_g_n')    = 12.0650;                  %% Gearbox Damping - negative direction [Nms/rad]
params.('d_l_n')    = 0;                        %% Torsion bar damping - negative direction [Nms/rad]
% Linear internal viscous friction
params.('d_mg')     = 252.5627;                 %% Gearbox internal damping [Nms/rad] (not identified)
params.('d_gl')     = 0;                        %% Torsion bar internal damping [Nms/rad] (not identified)
% Coulomb friction
params.('d_cm')     = 0;                        %% Motor Coulomb damping [Nm]
params.('d_cg')     = 3.3897;                   %% Gearbox Coulomb damping [Nm]
params.('d_cl')     = 0;                        %% Torsion bar Coulomb damping [Nm]
% Asymmetric Coulomb friction
params.('d_cm_n')   = 0;                        %% Motor Coulomb damping - negative direction [Nm]
params.('d_cg_n')   = 3.5562;                   %% Gearbox Coulomb damping - negative direction [Nm]
params.('d_cl_n')   = 0;                        %% Torsion bar Coulomb damping - negative direction [Nm]
% Stiction
params.('d_s')      = 1.8;                      %% Break away torque [Nm]
params.('v_s')      = 0.01;                     %% Stribeck velocity range [rad/s]
% Misc
params.('k_t')      = 0.0445;                   %% Torque constant [Nm/A]
params.('r')        = 0.0885;                   %% Armature resistance [Ohm]
params.('x')        = 0.000140;                 %% Armature inductance [H]
params.('Ts')       = 5e-4;                     %% Sampling time [s]
% Operating/max conditions
params.('v_0')      = 24;                       %% Operating voltage [V]
params.('i_p')      = 80;                       %% Peak current [A]
params.('dq_p')     = 10.66;                    %% Max. peak speed (output) [rad/s]
% Thermal parameters
params.('r_th1')    = 0.29;                     %% Thermal Resistance Windings to Housing [K/W]
params.('r_th2')    = 0.73;                     %% Thermal Resistance Housing to Air [K/W]
params.('T_thw')    = 0.56;                     %% Thermal Time Constant of the Windings [s]
params.('T_thm')    = 14456;                    %% Thermal Time Constant of the Motor [s]
params.('Tmp_WMax') = 155;                      %% Maximum Armature Temperature [�C]