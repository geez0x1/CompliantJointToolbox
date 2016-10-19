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
        %
        % Mechanical Properties
        %
        % Inertiae
        m        = 2;      % Actuator mass [kg]                                    (default: 2)
        I_m      = 0.5480; % Link referred motor rotor inertia [kg m^2]            (default: 0.5480)
        I_g      = 0.2630; % Link referred gear inertia [kg m^2]                   (default: 0.2630)
        I_b      = 0.0867  % Link referred torsion bar inertia [kg m^2]       	  (default: 0.0867)
        % Stiffnesses
        k_g      = 31e3;   % Gearbox stiffness [Nm/rad]                            (default: 31e3)
        k_b      = 10e3;   % Torsion bar stiffness [Nm/rad]                        (default: 10e3)
        % Linear viscous friction
        d_m      = 14.786; % Motor damping [Nms/rad]                               (default: 14.786)
        d_g      = 0;      % Gearbox damping [Nms/rad]                             (default: 0)
        d_b      = 0;      % Torsion bar damping [Nms/rad]                         (default: 0)
        % Asymmetric viscous friction
        d_m_n    = 14.786; % Motor Damping - negative direction [Nms/rad]          (default: 14.786)
        d_g_n    = 0;      % Gearbox Damping - negative direction [Nms/rad]        (default: 0)
        d_b_n    = 0;      % Torsion bar damping - negative direction [Nms/rad]    (default: 0)
        % Linear internal viscous friction
        d_mg     = 300;    % Gearbox internal damping [Nms/rad]                    (default: 300)
        d_gb     = 35;     % Torsion bar internal damping [Nms/rad]                (default: 35)
        % Coulomb friction
        d_cm     = 0.1858; % Motor Coulomb damping [Nm]                            (default: 0.1858)
        d_cg     = 0;      % Gearbox Coulomb damping [Nm]                          (default: 0)
        d_cb     = 0;      % Torsion bar Coulomb damping [Nm]                      (default: 0)
        % Asymmetric Coulomb friction
        d_cm_n   = 0.1858  % Motor Coulomb damping - negative direction [Nm]       (default: 0.1858)
        d_cg_n   = 0;      % Gearbox Coulomb damping - negative direction [Nm]     (default: 0)
        d_cb_n   = 0;      % Torsion bar Coulomb damping - negative direction [Nm] (default: 0)
        % Stiction
        d_s      = 0;      % Break away torque [Nm]                                (default: 0)
        v_s      = 0;      % Stribeck velocity range [rad/s]                       (default: 0)
        % Cogging
        cog_a1	 = 0;      % Cosine amplitude [Nm]                                 (default: 0)
        cog_a2	 = 0;      % Sine amplitude [Nm]                                   (default: 0)
        cog_f	 = 0;      % Spatial frequency [periods/revolution]                (default: 0)
        % Gear
        n        = 100;    % Gear ratio [.]                                        (default: 100)
        %
        % Electrical Properties
        %
        k_t      = 0.0453;  % Torque constant [Nm/A]                               (default: 0.0453)
        r        = 0.0885;  % Armature resistance at normal temperature [Ohm]      (default: 0.0885)
        x        = 1.4e-4;  % Armature inductance [H]                              (default: 1.4e-4)
        %
        % Operating/max conditions
        %
        v_0      = 24;      % Operating [V]                                        (default: 24)
        i_c      = 40;      % Max. continuous current [A]                          (default: 40)
        i_p      = 40;      % Peak current [A]                                     (default: 40)
        dq_c     = 5.86;    % Max. continuous speed (output)[rad/s]                (default: 5.86)
        dq_p     = 5.86;    % Max. peak speed (output) [rad/s]                     (default: 5.86)
        %
        % Thermal parameters
        %
        r_th1    = 0.29;    % Thermal Resistance Windings to Housing [K/W]         (default: 0.29)
        r_th2    = 3.45;    % Thermal Resistance Housing to Air [K/W]              (default: 3.45)
        T_thw    = 3.96;    % Thermal Time Constant of the Windings [s]            (default: 3.96)
        T_thm    = 1240;    % Thermal Time Constant of the Motor [s]               (default: 1240)
        Tmp_WMax = 120;     % Maximum Armature Temperature [°C]                    (default: 120)
        Tmp_ANom = 25;      % Normal Ambient Temperature [°C]                      (default: 25)
       
        % Desciptive Properties
        name                % Joint name
        paramName           % Parameter name
        modelName           % Model name
        nonlinearModelName  % Nonlinear model name

        % Misc
        Ts      = 1e-3;   % Sampling time [s]
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
            sysd    = tf(sys);
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
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also getStateSpaceD, getTFd, genericJoint, jointBuilder.
            
            % Get all properties
            props = properties(obj);
            
            % Blacklist (non-variable property names)
            blacklist = {   'verbose';
                            'debug';
                            'name';
                            'paramName';
                            'modelName';
                            'nonlinearModelName';
                        };
            
            % Get symbolic properties
            symProps    = setdiff(props, blacklist);
            nProps      = numel(symProps);
            
            % Set each symbolic property to symbolic
            for iProps = 1:nProps
                obj.(symProps{iProps}) = sym(symProps{iProps},'real');
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
