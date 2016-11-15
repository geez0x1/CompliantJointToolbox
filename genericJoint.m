% GENERICJOINT Abstract generic base class for joint models.
%
% This class is abstract meaning that instances of this class cannot be
% instantiated. Only subclasses of GENERICJOINT can be instantiated.
%
% Notes::
%
%
% Examples::
%
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also jointBuilder.

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

classdef genericJoint < handle
    
    properties (GetAccess = private)
        % Blacklist (non-variable property names)
        blacklist = {
            'name';
            'paramName';
            'modelName';
            'nonlinearModelName';
        };
    end
    
    properties (SetAccess = private)
        a_CU    = 0.0039;  % Resistance coefficient of copper [1/K]
        c_CU    = 380;     % Specific thermal capacitance of copper [J/kg/K]
        c_FE	= 460;     % Specific thermal capacitance of iron [J/kg/K]
    end
    
    properties (SetAccess = public)
        %
        % Mechanical Properties
        %
        % Dimensions
        diam 	= 100;      % Actuator diameter [mm]                                (default: 100)
        len  	= 150;      % Actuator length [mm]                                  (default: 150)
        % Inertiae
        m       = 2;        % Actuator mass [kg]                                    (default: 2)
        I_m  	= 0.5480;   % Link referred motor rotor inertia [kg m^2]            (default: 0.5480)
        I_g     = 0.2630;   % Link referred gear inertia [kg m^2]                   (default: 0.2630)
        I_l     = 0.0867    % Link referred load inertia [kg m^2]       	    	(default: 0.0867)
        % Stiffnesses
        k_g     = 31e3;     % Gearbox stiffness [Nm/rad]                            (default: 31e3)
        k_b     = 10e3;     % Torsion bar stiffness [Nm/rad]                        (default: 10e3)
        % Linear viscous friction
        d_m     = 14.786;   % Motor damping [Nms/rad]                               (default: 14.786)
        d_g     = 0;        % Gearbox damping [Nms/rad]                             (default: 0)
        d_l     = 0;        % Load damping [Nms/rad]                                (default: 0)
        % Asymmetric viscous friction
        d_m_n   = 14.786;   % Motor Damping - negative direction [Nms/rad]          (default: 14.786)
        d_g_n   = 0;        % Gearbox Damping - negative direction [Nms/rad]        (default: 0)
        d_l_n   = 0;        % Load damping - negative direction [Nms/rad]           (default: 0)
        % Linear internal viscous friction
        d_mg    = 300;      % Gearbox internal damping [Nms/rad]                    (default: 300)
        d_gl    = 35;       % Torsion bar internal damping [Nms/rad]                (default: 35)
        % Coulomb friction
        d_cm    = 0.1858;   % Motor Coulomb damping [Nm]                            (default: 0.1858)
        d_cg    = 0;        % Gearbox Coulomb damping [Nm]                          (default: 0)
        d_cl    = 0;        % Load Coulomb damping [Nm]                             (default: 0)
        % Asymmetric Coulomb friction
        d_cm_n  = 0.1858    % Motor Coulomb damping - negative direction [Nm]       (default: 0.1858)
        d_cg_n  = 0;        % Gearbox Coulomb damping - negative direction [Nm]     (default: 0)
        d_cl_n  = 0;        % Load Coulomb damping - negative direction [Nm]        (default: 0)
        % Stiction
        d_s     = 0;        % Break away torque [Nm]                                (default: 0)
        v_s     = 0;        % Stribeck velocity range [rad/s]                       (default: 0)
        % Cogging
        cog_a1	= 0;        % Cosine amplitude [Nm]                                 (default: 0)
        cog_a2	= 0;        % Sine amplitude [Nm]                                   (default: 0)
        cog_f	= 0;        % Spatial frequency [periods/revolution]                (default: 0)
        % Gear
        n       = 100;      % Transmission ratio [.]                                (default: 100)
        %
        % Electrical Properties
        %
        k_t  	= 0.0453;   % Torque constant [Nm/A]                                (default: 0.0453)
        r       = 0.0885;   % Armature resistance at normal temperature [Ohm]       (default: 0.0885)
        x       = 1.4e-4;   % Armature inductance [H]                               (default: 1.4e-4)
        %
        % Operating/max conditions
        %
        v_0     = 24;       % Operating [V]                                         (default: 24)
        i_p     = 80;       % Peak current (demagnetization point) [A]              (default: 80)
        dq_p    = 5.86;     % Max. peak speed (output) [rad/s]                      (default: 5.86)
        %
        % Thermal parameters
        %
        r_th1   = 0.29;     % Thermal Resistance Windings to Housing [K/W]          (default: 0.29)
        r_th2   = 3.45;     % Thermal Resistance Housing to Air [K/W]               (default: 3.45)
        T_thw   = 3.96;     % Thermal Time Constant of the Windings [s]             (default: 3.96)
        T_thm   = 1240;     % Thermal Time Constant of the Motor [s]                (default: 1240)
        Tmp_WMax = 120;  	% Maximum Armature Temperature [�C]                     (default: 120)
        Tmp_ANom = 25;      % Normal Ambient Temperature [�C]                       (default: 25)
        
        % Desciptive Properties
        name                % Joint name
        paramName           % Parameter name
        modelName           % Model name
        nonlinearModelName  % Nonlinear model name
        
        % Misc
        Ts      = 1e-3;     % Sampling time [s]
    end
    
    methods
        %__________________________________________________________________
        function this = genericJoint(params)
            % genericJoint Default constructor of the generic joint class.
            %
            %   gj = genericJoint(params)
            %
            % Inputs:
            %   params: parameter struct with fields corresponding to the
            %   genericJoint class parameters.
            %
            % Outputs:
            %   gj: the joint object gj.
            %
            % Notes::
            %   Since genericJoint is an abstract class, you can never manually
            %   call this constructor.
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also jointBuilder.
            
            % Gather information about the input parameter set
            parFields = fields(params);
            nFields = size(parFields,1);
            
            % Copy parameters into the class properties
            for iFields = 1:nFields
                % check if parameter is actually a property
                if isprop(this,parFields{iFields})
                    this.(parFields{iFields}) = params.(parFields{iFields});
                else
                    warning(['NOT A FIELD: ',parFields{iFields}, ' is not a field of genericJoint class.'])
                end
            end
            
        end
        
        
        %__________________________________________________________________
        
        function params = getParams(this)
            % GETPARAMS Return a struct with all joint parameters.
            %
            %   p = gj.getParams
            %
            % Inputs:
            %
            % Outputs:
            %   p: struct with fields identical to the parameters of the joint
            %      object gj.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also genericJoint, jointBuilder.
            
            params  = struct;
            p       = properties(this);
            
            % Put all properties
            for iP = 1:numel(p);
                params.(p{iP}) = this.(p{iP});
            end
            
        end
        
        %__________________________________________________________________
        %                     Secondary Properties
        %__________________________________________________________________
        %
        % We call those properties that are computed from the actual
        % physical acutator properties "Secondary Properties.
        
        %__________________________________________________________________
        function out = i_c(this)
            % I_C Continuous Current [A]
            %
            %   i_c = gj.i_c
            %
            % Inputs:
            %
            % Outputs:
            %   i_c: Value for the maximum permissible continuous current.
            %   The value is mainly determined by the Joule losses inside
            %   the motor windings at steady state current.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_r, p_rce, genericJoint, jointBuilder.
            
            % Shorthands
            aCU   = this.a_CU;          % Resistance coefficient for copper
            r_TA  = this.r;             % Winding resistance at normal ambient temperature
            T_A   = this.Tmp_ANom;      % Normal ambient temperature
            T_W   = this.Tmp_WMax;      % Maximum allowed winding temperature
            rth1  = this.r_th1;         % Thermal resistance Winding-Housing
            rth2  = this.r_th2;         % Thermal resistance Housing-Air
            
            dT = T_W-T_A;               % Allowed temperature rise
            
            % Compute maximum continuous current
            %
            % The temperature rise is mainly due to Joule losses PJ in the
            % windings. In steady state we obtain:
            % dT = (r_th1 + r_th2)*PJ,
            % With PJ = R(T_W)*I^2
            %
            % The winding resistance is temperature dependent:
            % R(T_W) = r_TA * (1 + aCU * dT)
            %
            % Solving for I^2 yields:
            iSquare = dT / (rth1 + rth2) / r_TA / (1 + aCU*dT);
            
            % The maximum continuous current thus computes to:
            out = sqrt(iSquare);
            
        end
        
        %__________________________________________________________________
        
        function out = t_r(this)
            % t_r Continuous torque [Nm]
            %
            %   t_r = gj.t_r
            %
            % Inputs:
            %
            % Outputs:
            %   t_r: Value for the rated operation torque obtained
            %   as the product of the maximum permissible continuous current
            %   i_c, the torque constant k_t as well as the transmission ratio n.
            %   It is the torque that can be continuously generated by the
            %   windings without overheating.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also i_c, t_p, genericJoint, jointBuilder.
            
            
            % Compute maximum continuous motor torque
            out = this.k_t * this.i_c * this.n;
        end
        
        %__________________________________________________________________
        
        function out = t_p(this)
            % T_P Peak torque [Nm]
            %
            %   t_p = gj.t_p
            %
            % Inputs:
            %
            % Outputs:
            %   t_p: Value for the peak torque obtained as the product of
            %   the maximum permissible peak current i_p, the torque constant
            %   k_t and the transmission ratio n.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also i_p, t_r, p_cm, genericJoint, jointBuilder.
            
            out = this.k_t * this.i_p * this.n;
        end
        
        %__________________________________________________________________
        
        function out = k_w(this)
            % K_W Speed constant [s^(-1)/V]
            %
            %   k_w = gj.k_w
            %
            % Inputs:
            %
            % Outputs:
            %   k_w: Value for the speed constant which equals the inverse
            %   of the torque constant k_t.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_r, genericJoint, jointBuilder.
            
            out = 1 / this.k_t;
        end

        %__________________________________________________________________
        
        function out = T_mech(this)
            % T_MECH Mechanical Time Constant [s]
            %
            %   T_mech = gj.T_mech
            %
            % Inputs:
            %
            % Outputs:
            %   T_mech: Mechanical time constant, that describes the time
            %   to arrive at 63 % of the steady state velocity when
            %   accelerating with constant armature voltage.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_r, p_rce, genericJoint, jointBuilder.
            
            out = (this.I_m + this.I_g + this.I_l) * this.r / this.n^2 / this.k_t^2;
        end
        
        %__________________________________________________________________
        
        function out = T_el(this)
            % T_EL Electrical Time Constant [s]
            %
            %   T_el = gj.T_el
            %
            % Inputs:
            %
            % Outputs:
            %   T_el: Electrical time constant, that describes the time
            %   to arrive at 63 % of the steady state current when
            %   applying a constant armature voltage.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_r, p_rce, genericJoint, jointBuilder.
            
            out = this.x / this.r;
        end
        
        %__________________________________________________________________
        function out = dq_0(this)
            % DQ_0 Zero-Torque Speed [rad/s]
            %
            %   dq_0 = gj.dq_0
            %
            % Inputs:
            %
            % Outputs:
            %   dq_0: Zero-Torque speed in rad/s;
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_r, p_ce, genericJoint, jointBuilder.
            
            out = this.k_w * this.v_0 / this.n;
        end
        
        %__________________________________________________________________
        function out = dq_NL(this)
            % DQ_NL No-Load Speed [rad/s]
            %
            %   dq_NL = gj.dq_NL
            %
            % Inputs:
            %
            % Outputs:
            %   dq_NL: Zero-Load speed in rad/s;
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_r, p_ce, genericJoint, jointBuilder.
            
            sumCoulomb = this.d_cm + this.d_cg + this.d_cl;
            sumViscous = this.d_m + this.d_g + this.d_l;
            
            out = ( this.dq_0 - this.dq_over_dm * sumCoulomb ) / (1 + this.dq_over_dm * sumViscous);
        end
        
        %__________________________________________________________________
        function out = t_NL(this)
            % t_NL No load torque [Nm]
            %
            %   t_NL = gj.t_NL
            %
            % Inputs:
            %
            % Outputs:
            %   iq_0: No load torque in Nm;
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_r, p_ce, genericJoint, jointBuilder.
             
            sumCoulomb = this.d_cm + this.d_cg + this.d_cl;
            sumViscous = this.d_m + this.d_g + this.d_l;
            
            
            out = sumCoulomb + sumViscous * this.dq_NL;
        end
        
        
        %__________________________________________________________________
        function out = i_NL(this)
            % i_NL No load current [A]
            %
            %   i_NL = gj.i_NL
            %
            % Inputs:
            %
            % Outputs:
            %   i_NL: No load current in A;
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_r, p_ce, genericJoint, jointBuilder.
             
            out = this.t_NL / this.n / this.k_t;
        end
        
        
        
        %__________________________________________________________________
        function out = dq_r(this)
            % DQ_r Rated speed [rad/s]
            %
            %   dq_r = gj.dq_r
            %
            % Inputs:
            %
            % Outputs:
            %   dq_r: Rated speed at which the motor turns, when the 
            %         maximum continuous current is applied at full supply 
            %         voltage rad/s.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_r, p_ce, genericJoint, jointBuilder.
            
            out = this.dq_0 - this.dq_over_dm * this.t_r;
        end
        
        
        %__________________________________________________________________
        
        function out = i_start(this)
            % I_START Starting current[A]
            %
            %   t_start = gj.i_start
            %
            % Inputs:
            %
            % Outputs:
            %   i_start: Current when nominal voltage is applied to the
            %   motor at rest. It is the current equivalent to the stall
            %   torque.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_r, t_stall, genericJoint, jointBuilder.
            
            out = this.v_0 / this.r;
        end
        
        %__________________________________________________________________
        
        function out = t_stall(this)
            % T_STALL Stall torque [Nm]
            %
            %   t_stall = gj.t_stall
            %
            % Inputs:
            %
            % Outputs:
            %   t_stall: Load torque in Nm at which the motor stops, if the
            %            full nominal voltage is applied, provided that the
            %            required current is delivered by the power source.
            %            
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also i_start, t_r, genericJoint, jointBuilder.
            
            out = this.i_start * this.k_t * this.n;
        end
        
        
        %__________________________________________________________________
        
        function out = dq_over_dm(this)
            % DQ_OVER_DM Motor slope [rad s^(-1) / Nm]
            %
            %   slope = gj.dq_over_dm
            %
            % Inputs:
            %
            % Outputs:
            %   slope: Slope defining the drop in motor speed per load torque
            %          increment.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_r, p_rce, genericJoint, jointBuilder.
            
            out =  this.dq_0 / this.t_stall;
        end
        
        %__________________________________________________________________
        
        function out = p_cm(this)
            % P_CM Continous power (mechanical) [W]
            %
            %   p_cm = gj.p_cm
            %
            % Inputs:
            %
            % Outputs:
            %   p_cm: Value for the rated mechanical continuous power 
            %   obtained as the product of rated speed and rated torque.
            %
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also dq_r, i_c, t_r, genericJoint, jointBuilder.
            
            out = this.dq_r * this.t_r;
        end
        
        %__________________________________________________________________
        
        function out = p_pm(this)
            % P_PM Peak power (mechanical) [W]
            %
            %   p_pm = gj.p_pm
            %
            % Inputs:
            %
            % Outputs:
            %   p_pm: Value for the mechanical peak power obtained as the 
            %   maximum surface are under the speed-torque curve. 
            %
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also dq_r, i_c, t_r, genericJoint, jointBuilder.
            
            % The maximum mechanical power occurs at half the stall torque.
            % Half the stall torque might however correspond to a winding
            % current beyond the peak current limit. Hence we have to
            % determine the minimum of the two:
            t_val = min(this.t_p, this.t_stall/2);
            
            % The maximum achievable speed at the torque defined by t_val
            % can be obtained from the speed-torque curve:
            v_tau = this.dq_0 - this.dq_over_dm * t_val;
            
            % The peak torque is the product of torque and speed.
            out = t_val * v_tau;
        end
        
        %__________________________________________________________________
        
        function out = p_ce(this)
            % P_CE Continous electrical power (electrical) [W]
            %
            %   p_ce = gj.p_ce
            %
            % Inputs:
            %
            % Outputs:
            %   p_ce: Value for the rated continuous power obtained as the
            %   product of operating voltage and maximum permissible continuous
            %   current.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also i_c, t_r, genericJoint, jointBuilder.
            
            out = this.v_0 * this.i_c;
        end
        
        %__________________________________________________________________
        
        function out = p_pe(this)
            % P_PE Peak electrical power (electrical) [W]
            %
            %   p_pe = gj.p_pe
            %
            % Inputs:
            %
            % Outputs:
            %   p_pe: Peak value for the electircal power obtained as the
            %   product of operating voltage and peak current. The 
            %   calculation thus assumes the motor at rest.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also i_c, t_r, genericJoint, jointBuilder.
            
            out = this.v_0 * this.i_p;
        end
        
        %__________________________________________________________________
        %                     Dynamic System Models
        %__________________________________________________________________
        %
        % The physical and secondary parameters are used to obtain danymic
        % system models for monitoring and control of the actuator.
        
        %__________________________________________________________________
        
        function sys = getStateSpace(obj)
            % GETSTATESPACE Get continuous time state-space representation of
            % the linear dynamics.
            %
            %   sys = gj.getStateSpace
            %
            % Inputs:
            %
            % Outputs:
            %   sys: continuous time state-space representation of the linear
            %   dynamics.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also getStateSpaceD, genericJoint, jointBuilder.
            [A, B, C, ~, ~, ~]  = obj.getDynamicsMatrices();
            D                   = 0;
            sys                 = ss(A, B, C, D);
        end
        
        %__________________________________________________________________
        function sysd = getStateSpaceD(this)
            % GETSTATESPACE Get discrete time state-space representation of
            % the linear dynamics.
            %
            %   sys = gj.getStateSpaceD
            %
            % Inputs:
            %
            % Outputs:
            %   sys: discrete time state-space representation of the linear
            %   dynamics.
            %
            % Notes::
            %   The function automatically uses the sampling time specified by
            %   the object parameters.
            %   Beyond that, the function discretizes the continous time
            %   description of the dynamics. Therefore implicitly calls
            %   getStateSpace and c2d afterwards.
            %
            % TODO::
            %  - Optionally allow different c2d methods
            %  - Optionally allow different sampling times
            %  - In particular, allow for both existing discrete time state
            %    space descriptions.
            %
            %   Since there exist two generally different forms of discrete
            %   time state space models, the functionality should be improved
            %   to allow decisions for any of them made by the user.
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn, (jorn.malzahn@iit.it)
            %  Wesley Roozing, (wesley.roozing@iit.it)
            %
            % See also getStateSpace, ss, c2d.
            
            sys     = this.getStateSpace();
            sysd    = c2d(sys, this.Ts, 'tustin');
        end
        
        %__________________________________________________________________
        function sys = getTF(this)
            % GETTF Get continuous time transfer function representation of the
            % linear dynamics.
            %
            %   sys = gj.getTF
            %
            % Inputs:
            %
            % Outputs:
            %   sys: continuous time transfer function representation of the linear
            %   dynamics.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn, jorn.malzahn@iit.it
            %  Wesley Roozing, wesley.roozing@iit.it
            %
            % See also getStateSpace, getTFd, genericJoint, jointBuilder.
            
            if ~this.isSym
                sys     = this.getStateSpace();
                sys     = tf(sys);
            else
                syms s
                [A B C D ] = this.getDynamicsMatrices;
                E = sym(eye(size(A,1)));
                sys = simplify( C*inv(s*E-A)*B );
                
            end
        end
        
        %__________________________________________________________________
        function sysd = getTFd(obj)
            % GETTFD Get discrete time transfer function representation of the
            % linear dynamics.
            %
            %   sys = gj.getTFd
            %
            % Inputs:
            %
            % Outputs:
            %   sys: discrete time transfer function representation of the linear
            %   dynamics.
            %
            % Notes::
            %   This function calls getStateSpaceD and converts the resulting
            %   state space model into a transfer function. For details about
            %   the behavior please look at getStateSpaceD.
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn, jorn.malzahn@iit.it
            %  Wesley Roozing, wesley.roozing@iit.it
            %
            % See also getStateSpaceD, getTFd, genericJoint, jointBuilder.
            sys     = obj.getStateSpaceD();
            sysd    = tf(sys);
        end
        
        %__________________________________________________________________
        function flag = isSym(this)
        
            % Get all properties
            props = properties(this);            
            
            % Get symbolic properties
            symProps    = setdiff(props, this.blacklist);
            nProps      = numel(symProps);
            
            flag = 0;
            
            for iProps = 1:nProps
                
                if isa(this.(symProps{iProps}),'sym')
                    flag = 1;   % Set the flag.
                    break;      % We found at least one symbolic property, we can break the loop hera.
                end
            end           
        
        end
        
        %__________________________________________________________________
        function makeSym(this,varargin)
            % MAKESYM Return a symbolic copy of the joint model.
            %
            %   gj.makeSym({doSparse})
            %
            % Inputs:
            %   doSparse:  If false, all parameters are turned into generic 
            %              symbolic  ariables. If true (default), 
            %              parameters that are exactly zero remain zero in 
            %              the symbolic version of the model.
            % Outputs:
            %   gj_sym: Converts of the original object into a symbolic model
            %   with all properties being symbolic variables.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also getStateSpaceD, getTFd, genericJoint, jointBuilder.
            
            doSparse = 1;
            if nargin == 2
                doSparse = varargin{1};
            end
            
            % Get all properties
            props = properties(this);
            
            
            % Get symbolic properties
            symProps    = setdiff(props, this.blacklist);
            nProps      = numel(symProps);
            
            % Set each symbolic property to symbolic
            for iProps = 1:nProps
                
                if doSparse && this.(symProps{iProps}) == 0
                    this.(symProps{iProps}) = 0;
                else
                    this.(symProps{iProps}) = sym(symProps{iProps},'real');
                end
            end

        end
        
        %__________________________________________________________________
        function makeNum(this)
            % MAKESYM Converts all properties into numeric variables by
            % reloading the default parameters from the class file.
            %
            %   gj.makeNum
            %
            %
            % Inputs:
            %
            % Outputs:
            %
            %
            % Notes::
            %   This method is effectively a wrapper to 'resetParams'. It
            %   exists as a complement to the method 'makeSym'.
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also makeSym, resetParams, genericJoint, jointBuilder.
            
            this.resetParams;
            
        end
        
        %__________________________________________________________________
        function resetParams(this, varargin)
            % RESETPARAMS Reload the original parameters from the class file.
            %
            %   gj.resetParams([params])
            %
            %
            % Inputs:
            %
            % Outputs:
            %
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also getStateSpaceD, getTFd, genericJoint, jointBuilder.
            
            
            % load default parameters from class file
            params = eval(this.name);
            
            % Gather information about the input parameter set
            parFields = fields(params);
            nFields = size(parFields,1);
            
            % Copy parameters into the class properties
            for iFields = 1:nFields
                % check if parameter is actually a property
                if isprop(this,parFields{iFields})
                    this.(parFields{iFields}) = params.(parFields{iFields});
                else
                    warning(['NOT A FIELD: ',parFields{iFields}, ' is not a field of genericJoint class.'])
                end
            end
        end
        
    end
    
    %__________________________________________________________________
    % Abstract methods - to be implemented by subclasses.
    methods(Abstract)
        %__________________________________________________________________
        % GETDYNAMICSMATRICES Get Dynamics Matrices for the Dynamics
        %
        %   [A, B, C, I, D, K] = gj.getDynamicsMatrices
        %
        % Inputs:
        %
        % Outputs:
        %   A: continuous time system matrix
        %   B: continuous time input matrix
        %   C: continuous time output matrix
        %   I: inertia matrix
        %   D: damping matrix
        %   K: stiffness matrix
        %
        % Notes::
        %   This is an abstract method. It has to be implemented by
        %   subclasses.
        %
        % Examples::
        %
        %
        % Author::
        %  Joern Malzahn
        %  Wesley Roozing
        %
        % See also getNonlinearDynamics, getStateSpace, getTFd, genericJoint.
        [A, B, C, I, D, K] = getDynamicsMatrices(obj)
        
        %__________________________________________________________________
        % GETDYNAMICSMATRICES Compute nonlinear dynamics
        %
        %   tau = gj.getNonlinearDynamics(x,dx)
        %
        % Inputs:
        %   x: position variable
        %  dx: temporal derivative of the position variable
        %
        % Outputs:
        %   tau: Generalized torque due to the nonlinear dynamics.
        %
        % Notes::
        %   This is an abstract method. It has to be implemented by
        %   subclasses.
        %
        % Examples::
        %
        %
        % Author::
        %  Joern Malzahn
        %  Wesley Roozing
        %
        % See also getDynamicsMatrices, getStateSpace, getTFd, genericJoint.
        tau = getNonlinearDynamics(obj, x, dx)
        
    end
  
end
