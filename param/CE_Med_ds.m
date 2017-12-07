%CE_MED_DS Parameter script for the Centauro Medium Motor with encoder based torque sensor
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

% Motor type
% Color:            Orange &#x1F4D9
% Shaft inertia:    5.5375834	[kg mm^2]
% Motor inertia:    25.2        [kg mm^2]
% Gearbox inertia:  11.2        [kg mm^2]
% Motor and shaft are lumped together. Non-load inertia values are
% scaled to 110% (following experiments), to account for other
% (mounting) components.

% Components used
% Motor:            Kollmorgen TBMS6025-B (design voltage 24V, used at 48V)
% Gearbox:          CPL20-2A

% Friction calculations
% Motor:    Viscous friction estimated directly from datasheet:
%           d_rads = 5.36e-3 * (60/(1000*2*pi))
%           Coulomb friction: d_cm = 0.033 * n
% Gearbox:  Viscous friction estimated from 'no load running torque' at 3500rpm
%           input speed, 20 degrees Celcius, after removing compensation values
%           (table 29.1) and no load starting torque (table 29.2).
%           tau = (16-0.3-2.9)/100; omega = (3500/60)*(2*pi); d = n^2 * tau / omega;
%           'No load starting torque' used for Coulomb friction
%           values: d_cg = 2.9/100 * n

% Inertiae
params.('I_m')      = 3.3812e-05 * n^2;   	%% Motor rotor inertia [kg m^2] (motor+shaft)
params.('I_g')      = 1.2320e-05 * n^2; 	%% Gear inertia [kg m^2]
params.('I_l')      = 1.0;                  %% Load inertia [kg m^2] % TO BE UPDATED!
% Stiffnesses
params.('k_g')      = 25e3;                 %% Gearbox stiffness [Nm/rad] (using K_2 from datasheet, valid for 7..25 Nm)
params.('k_b')      = 6500;                 %% Torsion bar stiffness [Nm/rad]
% Linear viscous friction
params.('d_m')      = 5.1184e-05 * n^2;     %% Motor Damping [Nms/rad]
params.('d_g')      = 3.4923e-04 * n^2;     %% Gearbox damping [Nms/rad]
params.('d_l')      = 0.0;                  %% Load damping [Nms/rad] % TO BE UPDATED!
% Asymmetric viscous friction - negative direction
params.('d_m_n')    = 5.1184e-05 * n^2;     %% Motor Damping - negative direction [Nms/rad]
params.('d_g_n')    = 3.4923e-04 * n^2;     %% Gearbox Damping - negative direction [Nms/rad]
params.('d_l_n')    = 0.0;                  %% Load damping - negative direction [Nms/rad] % TO BE UPDATED!
% Linear internal viscous friction
params.('d_mg')     = 300.0;                %% Gearbox internal damping [Nms/rad] % TO BE UPDATED!
params.('d_gl')     = 0.0;                  %% Torsion bar internal damping [Nms/rad] % TO BE UPDATED!
% Coulomb friction
params.('d_cm')     = 0.033 * n;            %% Motor Coulomb damping [Nm]
params.('d_cg')     = 0.0290 * n;           %% Gearbox Coulomb damping [Nm]
params.('d_cl')     = 0.0;                  %% Load Coulomb damping [Nm] % TO BE UPDATED!
% Asymmetric Coulomb friction - negative direction
params.('d_cm_n')   = 0.033 * n;            %% Motor Coulomb damping - negative direction [Nm]
params.('d_cg_n')   = 0.0290 * n;           %% Gearbox Coulomb damping - negative direction [Nm]
params.('d_cl_n')   = 0.0;                  %% Load Coulomb damping - negative direction [Nm] % TO BE UPDATED!
% Stiction
params.('d_s')      = 0.0;                  %% Break away torque [Nm] % TO BE UPDATED!
params.('v_s')      = 0.0;                  %% Stribeck velocity range [rad/s] % TO BE UPDATED!
% Torque ripple sources
params.('rip_types')= [];                   %% Torque ripple types (see torque_ripple())
params.('rip_a1')   = [];                   %% Cosine amplitudes ([Nm] and/or [], see torque_ripple())
params.('rip_a2')   = [];                   %% Sine amplitudes [Nm] ([Nm] and/or [], see torque_ripple())
params.('rip_f')    = [];                   %% Spatial frequencies [periods/revolution]
% Misc
params.('n')        = n;                    %% Gear ratio []
params.('k_t')      = 0.078;                %% Torque constant [Nm/A]
params.('r')        = 0.318;             	%% Armature resistance [Ohm]
params.('x')        = 0.20e-3;              %% Armature inductance [H] % TO BE UPDATED!
params.('Ts')       = 1e-3;                 %% Sampling time [s]
% Operating/max conditions
params.('v_0')      = 48;                   %% Operating voltage [V]
params.('i_p')      = 0.0;                  %% Peak current [A] % TO BE UPDATED!
params.('dq_p')     = 0.0;                  %% Max. peak speed (output) [rad/s] % TO BE UPDATED!
