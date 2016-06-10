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
    
    properties
        verbose % verbose flag
        debug   % debug flag
    end
    
    properties (SetAccess = public)
        % Mechanical Properties
        % Inertiae
        I_m     % Motor rotor inertia [kg m^2] (link side)
        I_g     % Motor-side gear inertia [kg m^2] (link side)
        I_b     % Torsion bar inertia [kg m^2] (link side)
        % Stiffnesses
        k_g     % Gearbox stiffness [Nm/rad]
        k_b     % Torsion bar stiffness [Nm/rad]
        % Linear viscous friction
        d_m     % Motor Damping [Nms/rad]
        d_g     % Gearbox damping [Nms/rad]
        d_b     % Torsion bar damping [Nms/rad]
        % Asymmetric viscous friction
        d_m_n   % Motor Damping - negative direction [Nms/rad]
        d_g_n	% Gearbox Damping - negative direction [Nms/rad]
        d_b_n	% Torsion bar damping - negative direction [Nms/rad]
        % Linear internal viscous friction
        d_mg	% Gearbox internal damping [Nms/rad]
        d_gb	% Torsion bar internal damping [Nms/rad]
        % Coulomb friction
        d_cm    % Motor Coulomb damping [Nm]
        d_cg    % Gearbox Coulomb damping [Nm]
        d_cb    % Torsion bar Coulomb damping [Nm]
        % Asymmetric Coulomb friction
        d_cm_n  % Motor Coulomb damping - negative direction [Nm]
        d_cg_n  % Gearbox Coulomb damping - negative direction [Nm]
        d_cb_n  % Torsion bar Coulomb damping - negative direction [Nm]
        % Misc
        n       % Gear ratio []
        k_t     % Torque constant [Nm/A]
        r       % Armature resistance [Ohm]
        x       % Armature inductance [H]
        Ts      % Sampling time [s]
        
        v_0     % Operating [V]
        i_c     % Max. continuous current [A]
        i_p     % Peak current [A]
        dq_c    % Max. continuous speed (output)[rad/s]
        dq_p    % Max. peak speed (output) [rad/s]
            
        % Desciptive Properties
        name                % Joint name
        paramName           % Parameter name
        modelName           % Model name
        nonlinearModelName	% Nonlinear model name
        
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
            %  Joern Malzahn, jorn.malzahn@iit.it
            %  Wesley Roozing, wesley.roozing@iit.it
            %
            % See also jointBuilder, WMBigJoint.
            
            % Apply properties
            this.verbose	= params.verbose;
            this.debug      = params.debug;
            
            % Mechanical Properties
            % Inertiae
            this.I_m	= params.I_m;       % Motor rotor inertia [kg m^2] (link side)
            this.I_g    = params.I_g;       % Motor-side gear inertia [kg m^2] (link side)
            this.I_b    = params.I_b;       % Torsion bar inertia [kg m^2] (link side)
            % Stiffnesses
            this.k_g    = params.k_g;       % Gearbox stiffness [Nm/rad]
            this.k_b    = params.k_b;       % Torsion bar stiffness [Nm/rad]
            % Linear viscous friction
            this.d_m    = params.d_m;       % Motor Damping [Nms/rad]
            this.d_g    = params.d_g;       % Gearbox damping [Nms/rad]
            this.d_b    = params.d_b;       % Torsion bar damping [Nms/rad]
            % Asymmetric viscous friction
            this.d_m_n 	= params.d_m_n; 	% Motor Damping - negative direction [Nms/rad]
            this.d_g_n 	= params.d_g_n;   	% Gearbox Damping - negative direction [Nms/rad]
            this.d_b_n 	= params.d_b_n;   	% Torsion bar damping - negative direction [Nms/rad]
            % Linear internal viscous friction
            this.d_mg 	= params.d_mg;   	% Gearbox internal damping [Nms/rad]
            this.d_gb 	= params.d_gb;   	% Torsion bar internal damping [Nms/rad]
            % Coulomb friction
            this.d_cm   = params.d_cm;      % Motor Coulomb damping [Nm]
            this.d_cg   = params.d_cg;      % Gearbox Coulomb damping [Nm]
            this.d_cb   = params.d_cb;      % Torsion bar Coulomb damping [Nm]
            % Asymmetric Coulomb friction
            this.d_cm_n = params.d_cm_n;    % Motor Coulomb damping - negative direction [Nm]
            this.d_cg_n = params.d_cg_n;	% Gearbox Coulomb damping - negative direction [Nm]
            this.d_cb_n = params.d_cb_n; 	% Torsion bar Coulomb damping - negative direction [Nm]
            % Misc
            this.n      = params.n;         % Gear ratio []
            this.k_t	= params.k_t;       % Torque constant [Nm/A]
            this.r      = params.r;         % Armature resistance [Ohm]
            this.x      = params.x;         % Armature inductance [H]
            this.Ts     = params.Ts;        % Sampling time [s]
            
            this.v_0    = params.v_0;       % Operating voltage [V]
            this.i_c    = 65;%params.i_c;       % Max. continuous current [A]
            this.i_p    = 65;%params.i_p;       % Peak stall current [A]
            this.dq_c    = params.dq_c;     % Max. continuous speed (output)[rad/s]
            this.dq_p    = params.dq_p;     % Max. peak speed (output) [rad/s]
            
            % Desciptive Properties
            this.name               = params.name;                  % Joint descriptive name
            this.paramName          = params.paramName;             % Parameter name
            this.modelName          = params.modelName;             % Model name
            this.nonlinearModelName	= params.nonlinearModelName;    % Nonlinear model name
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
            %  Joern Malzahn, jorn.malzahn@iit.it
            %  Wesley Roozing, wesley.roozing@iit.it
            %
            % See also genericJoint, jointBuilder, WMBigJoint.
            
            params	= struct;
            p       = properties(this);
            
            % Put all properties
            for iP = 1:numel(p);
                params.(p{iP}) = this.(p{iP});
            end
            
        end
        
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
            %  Joern Malzahn, jorn.malzahn@iit.it
            %  Wesley Roozing, wesley.roozing@iit.it
            %
            % See also getStateSpaceD, genericJoint, jointBuilder.
            [A, B, C, ~, ~, ~] = obj.getDynamicsMatrices();
            D = 0;
            sys = ss(A, B, C, D);
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
        function sys = getTF(obj)
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
            sys     = obj.getStateSpace();
            sys     = tf(sys);
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
            sysd	= tf(sys);
        end
        
        
        %__________________________________________________________________
        function makeSym(obj)
            % MAKESYM Return a symbolic copy of the joint model.
            %
            %   gj_sym = gj.makeSym
            %
            % Inputs:
            %
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
            %  Joern Malzahn, jorn.malzahn@iit.it
            %  Wesley Roozing, wesley.roozing@iit.it
            %
            % See also getStateSpaceD, getTFd, genericJoint, jointBuilder.
            
            % Get all properties
            props = properties(obj);
            
            % Blacklist (non-variable property names)
            blacklist = {'verbose';
                'debug';
                'name';
                'paramName';
                'modelName';
                'nonlinearModelName';
                };
            
            % Get symbolic properties
            symProps	= setdiff(props, blacklist);
            nProps      = numel(symProps);
            
            % Set each symbolic property to symbolic
            for iProps = 1:nProps
                obj.(symProps{iProps}) = sym(symProps{iProps},'real');
            end
        end
        
        
        function out = t_c(obj)
            % T_C Continuous stall torque [Nm]
            %
            %   t_c = gj.t_c
            %
            % Inputs:
            %
            % Outputs:
            %   t_c: Value for the continuous operation stall torque obtained
            %   as the product of the maximum permissible continuous current
            %   i_c, the torque constant k_t as well as the transmission ratio n.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn, jorn.malzahn@iit.it
            %
            % See also t_p, p_rce, genericJoint, jointBuilder.
            
            out = obj.k_t*obj.i_c*obj.n;
            
            
        end
        
        function out = t_p(obj)
            % T_P Peak stall torque [Nm]
            %
            %   t_p = gj.t_p
            %
            % Inputs:
            %
            % Outputs:
            %   t_p: Value for the peak stall torque obtained as the product of
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
            %  Joern Malzahn, jorn.malzahn@iit.it
            %
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            out = obj.k_t*obj.i_p*obj.n;
        end
        
        function out = k_w(obj)
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
            %  Joern Malzahn, jorn.malzahn@iit.it
            %
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            out = 1/obj.k_t;
        end

        function out = dq_0(obj)
            % DQ_0 No load speed [rad/s]
            %
            %   dq_0 = gj.dq_0
            %
            % Inputs:
            %
            % Outputs:
            %   dq_0: No load speed in rad/s;
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn, jorn.malzahn@iit.it
            %
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            out = obj.k_w * obj.v_0/obj.n;
        end
       
        function out = dq_over_dm(obj)
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
            %  Joern Malzahn, jorn.malzahn@iit.it
            %
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            out =  obj.dq_0/obj.t_stall;
        end
        
        function out = dq_r(obj)
            % DQ_r rated speed [rad/s]
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
            %  Joern Malzahn, jorn.malzahn@iit.it
            %
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            out = obj.dq_0 - obj.dq_over_dm * obj.t_c;
        end
        
        
        function out = p_rce(obj)
            % P_RCE Rated continous power (electrical) [W]
            %
            %   p_rce = gj.p_rce
            %
            % Inputs:
            %
            % Outputs:
            %   p_rce: Value for the rated continuous power obtained as the
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
            %  Joern Malzahn, jorn.malzahn@iit.it
            %
            % See also i_c, i_p, genericJoint, jointBuilder.
            out = obj.v_0*obj.i_c;
            
        end
            
        function out = p_rcm(obj)
            % P_RCM Rated continous power (mechanical) [W]
            %
            %   p_rcm = gj.p_rcm
            %
            % Inputs:
            %
            % Outputs:
            %   p_rcm: Value for the rated continuous power obtained as the
            %   product of continuous speed and continuous torque.
            %   
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn, jorn.malzahn@iit.it
            %
            % See also i_c, i_p, genericJoint, jointBuilder.
            out = obj.dq_r*obj.t_c;
            
        end
        
        function out = p_peakm(obj)
            % P_PEAKM Peak continous power (mechanical) [W]
            %
            %   p_peakm = gj.p_peakm
            %
            % Inputs:
            %
            % Outputs:
            %   p_peakm: Value for the peak mechanicalpower obtained as the
            %   product of peak speed and peak torque.
            %   
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn, jorn.malzahn@iit.it
            %
            % See also i_c, i_p, genericJoint, jointBuilder.
            out = obj.dq_p*obj.t_p;
        end
        
        function out = t_stall(obj)
            % T_STALL Stall torque [Nm]
            %
            %   t_stall = gj.t_stall
            %
            % Inputs:
            %
            % Outputs:
            %   t_stall: Load torque in Nm at which the motor stops, if the 
            %            full  operating voltage is applied.
            %
            % Notes::
            %
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn, jorn.malzahn@iit.it
            %
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            out = obj.v_0/obj.r*obj.k_t*obj.n;
        end
        
                
        function h = draw_speed_torque_curve(obj)
            % DRAW_SPEED_TORQUE_CURVE Displays speed-torque-curve in a
            % figure.
            %
            %   p_rce = gj.p_rce
            %
            % Inputs:
            %
            % Outputs:
            %   p_rce: Value for the rated continuous power obtained as the
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
            %  Joern Malzahn, jorn.malzahn@iit.it
            %
            % See also i_c, i_p, genericJoint, jointBuilder.
            
            nVals = 100;
            xmax = 1.5*obj.t_c;
            ymax = 1.02*obj.dq_c;
            h = figure;
            hold on
            
                       
            %% torque speed line
            slope = obj.dq_over_dm;
            mVals = (0:1/100:1)*obj.t_stall;
            linCurve =  obj.dq_0 - slope*mVals;
            plot(mVals, linCurve,'k','DisplayName','Torque-Speed Line')
            
            %% Nominal operating point
            plot(obj.t_c,obj.dq_r,'bo', 'DisplayName', 'Nominal Operating Point')
            
            %% Friction
%             speedVals = obj.v_0 * obj.k_w - obj.dq_over_dm*mVals;
            speedVals = (0:1/nVals:1)*obj.dq_c;
            Mc = obj.d_cm + obj.d_cg + obj.d_cb; % Static part
            Mv = (obj.d_m + obj.d_g + obj.d_b)*speedVals;
            Mf = Mc + Mv;
            plot(Mf,speedVals,'r:','DisplayName','Friction Torque')
            

            
            %% Plot limits
            speedVals = obj.p_rcm./mVals;
            plot(mVals,speedVals,'r.', 'DisplayName','Rated Mechanical Power')
            speedVals = obj.p_peakm./mVals;
            plot(mVals,speedVals,'b.', 'DisplayName','Peak Mechanical Power')
            
            plot([0,xmax], obj.dq_c*[1,1],'k--','DisplayName','Maximum Continous Speed')
            plot(obj.t_c*[1,1], [0,ymax],'k:','DisplayName','Maximum Continous Torque')
            
            
            %% Annotations and Figure Style
            xlim([0,xmax]);
            ylim([0,ymax]);
            xlabel('torque [Nm]')
            ylabel('speed [rad/s]')
            box on
            
            legend show
                
            
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
        %  Joern Malzahn, jorn.malzahn@iit.it
        %  Wesley Roozing, wesley.roozing@iit.it
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
        %  Joern Malzahn, jorn.malzahn@iit.it
        %  Wesley Roozing, wesley.roozing@iit.it
        %
        % See also getDynamicsMatrices, getStateSpace, getTFd, genericJoint.
        tau = getNonlinearDynamics(obj, x, dx)
    end
    
    
end

