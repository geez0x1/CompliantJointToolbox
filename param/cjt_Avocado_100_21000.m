% Parameter script for a TREE Robotics Avocado actuator.
%
% The actuator uses a 100:1 transmission ratio and has a stiffness of 21000 Nm/rad.
% For more details and CAD files of the actuator visit www.treerobotics.eu.
%
% Notes::
%  All inertiae/damping is reflected to link side using n^2
%  The parameters in this script are partially obtained from experimental identification.
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also genericJoint, full_dyn.

% Copyright (C) since 2016, by Joern Malzahn, Wesley Roozing
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
params.('n')        =        100;  %% Gear transmission ratio [.]
% Dimensions
params.('diam')     =         82;  %% Actuator diameter [mm] 
params.('len')      =        127;  %% Actuator length [mm] 
% Inertia
params.('m')        =          1;  %% Actuator mass [kg] 
params.('I_m')      =     0.0941;  %% Motor rotor inertia [kg m^2]
params.('I_g')      =     0.0540;  %% Gear inertia [kg m^2]
params.('I_l')      =     0.0001;  %% Torsion bar inertia [kg m^2]
% Stiffnesses
params.('k_g')      =       8400;  %% Gearbox stiffness [Nm/rad]
params.('k_b')      =      21000;  %% Torsion bar stiffness [Nm/rad]
% Linear viscous friction
params.('d_m')      =     0.0022;  %% Motor Damping [Nms/rad]
params.('d_g')      =     3.4375;  %% Gearbox damping [Nms/rad]
params.('d_l')      =     0.0000;  %% Torsion bar damping [Nms/rad]
% Asymmetric viscous friction
params.('d_m_n')    =     0.0022;  %% Motor Damping - negative direction [Nms/rad]
params.('d_g_n')    =     3.4375;  %% Gearbox Damping - negative direction [Nms/rad]
params.('d_l_n')    =     0.0000;  %% Torsion bar damping - negative direction [Nms/rad]
% Linear internal viscous friction
params.('d_mg')     =   100.0000;  %% Gearbox internal damping [Nms/rad] (not identified)
params.('d_gl')     =     0.0000;  %% Torsion bar internal damping [Nms/rad] (not identified)
% Coulomb friction
params.('d_cm')     =     1.5000;  %% Motor Coulomb damping [Nm]
params.('d_cg')     =     3.3000;  %% Gearbox Coulomb damping [Nm]
params.('d_cl')     =     0.0000;  %% Torsion bar Coulomb damping [Nm]
% Asymmetric Coulomb friction
params.('d_cm_n')   =     1.5000;  %% Motor Coulomb damping - negative direction [Nm]
params.('d_cg_n')   =     3.3000;  %% Gearbox Coulomb damping - negative direction [Nm]
params.('d_cl_n')   =     0.0000;  %% Torsion bar Coulomb damping - negative direction [Nm]
% Electrical Properties
params.('k_t')      =     0.0410;  %% Torque constant [Nm/A]
params.('r')        =     0.6640;  %% Armature resistance [Ohm]
params.('x')        =     0.0003;  %% Armature inductance [H]
params.('p')        =     4.0000;  %% Number of pole pairs [.]
params.('Ts')       =     0.0010;  %% Sampling time [s]
params.('Ts_elec')  =     0.0001;  %% Sampling time for the electrical subsystem [s]
% Operating/max conditions
params.('v_0')      =         48;  %% Operating voltage [V]
params.('i_p')      =         20;  %% Peak current [A]
params.('dq_p')     =    11.7169;  %% Max. reachable speed (output) [rad/s]
% Thermal parameters
params.('r_th1')    =     3.8619;  %% Thermal Resistance Windings to Housing (theoretical value!)[K/W]
params.('r_th2')    =    40.1831;  %% Thermal Resistance Housing to Air (theoretical value!) [K/W]
params.('T_thw')    =    88.2581;  %% Thermal Time Constant of the Windings  (theoretical value!) [s]
params.('T_thm')    =  2611.7046;  %% Thermal Time Constant of the Motor  (theoretical value!) [s]
params.('Tmp_WMax') =        155;  %% Maximum Armature Temperature [�C]
