classdef dataSheetGenerator
    %DATASHEETGENERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        paramStruct
    end
    
    methods
        function this = dataSheetGenerator(paramStruct)
            
            this.paramStruct = paramStruct;
            
        end
        
        function out = t_c(this)
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
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_p, p_rce, genericJoint, jointBuilder.
            
            pStruct = this.paramStruct; % Shorthand
            out = pStruct.k_t * pStruct.i_c * pStruct.n;
        end
        
        function out = t_p(this)
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
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            pStruct = this.paramStruct; % Shorthand
            out = pStruct.k_t * pStruct.i_p * pStruct.n;
        end
        
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
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            out = 1 / this.paramStruct.k_t;
        end

        function out = dq_0(this)
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
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            pStruct = this.paramStruct; % Shorthand
            out = this.k_w * pStruct.v_0 / pStruct.n;
        end
       
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
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            pStruct = this.paramStruct; % Shorthand
            out =  this.dq_0 / this.t_stall;
        end
        
        function out = dq_r(this)
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
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            pStruct = this.paramStruct; % Shorthand
            out = this.dq_0 - this.dq_over_dm * this.t_c;
        end
        
        
        function out = p_rce(this)
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
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also i_c, i_p, genericJoint, jointBuilder.
            
            pStruct = this.paramStruct; % Shorthand
            out = pStruct.v_0 * pStruct.i_c;
        end
            
        function out = p_rcm(this)
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
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also i_c, i_p, genericJoint, jointBuilder.
            
            pStruct = this.paramStruct; % Shorthand
            out = this.dq_r * this.t_c;
        end
        
        function out = p_peakm(this)
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
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also i_c, i_p, genericJoint, jointBuilder.
            
            pStruct = this.paramStruct; % Shorthand
            out = pStruct.dq_p * this.t_p;
        end
        
        function out = t_stall(this)
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
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            pStruct = this.paramStruct; % Shorthand
            out = pStruct.v_0 / pStruct.r * pStruct.k_t * pStruct.n;
        end
        
                
        function h = draw_speed_torque_curve(this)
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
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also i_c, i_p, genericJoint, jointBuilder.
            
            pStruct = this.paramStruct; % Shorthand
            
            nVals = 100;
            xmax = 1.5 * this.t_c;
            ymax = 1.02 * pStruct.dq_c;
            h = figure;
            hold on
                       
            % torque speed line
            slope = this.dq_over_dm;
            mVals = (0:1/100:1) * this.t_stall;
            linCurve =  this.dq_0 - slope * mVals;
            plot(mVals, linCurve, 'k', 'DisplayName', 'Torque-Speed Line')
            
            % Nominal operating point
            plot(this.t_c, this.dq_r, 'bo', 'DisplayName', 'Nominal Operating Point')
            
            % Friction
            %speedVals = obj.v_0 * obj.k_w - obj.dq_over_dm*mVals;
            speedVals = (0:1/nVals:1) * pStruct.dq_c;
            Mc = pStruct.d_cm + pStruct.d_cg + pStruct.d_cb; % Static part
            Mv = (pStruct.d_m + pStruct.d_g + pStruct.d_b) * speedVals;
            Mf = Mc + Mv;
            plot(Mf, speedVals, 'r:', 'DisplayName', 'Friction Torque')
            
            % Plot limits
            speedVals = this.p_rcm ./ mVals;
            plot(mVals,speedVals, 'r.', 'DisplayName', 'Rated Mechanical Power')
            speedVals = this.p_peakm ./ mVals;
            plot(mVals,speedVals, 'b.', 'DisplayName', 'Peak Mechanical Power')
            plot([0,xmax], pStruct.dq_c * [1,1], 'k--', 'DisplayName', 'Maximum Continous Speed')
            plot(this.t_c*[1,1], [0,ymax], 'k:', 'DisplayName', 'Maximum Continous Torque')
            
            % Annotations and Figure Style
            xlim([0,xmax]);
            ylim([0,ymax]);
            xlabel('torque [Nm]')
            ylabel('speed [rad/s]')
            box on
            
            legend show;
        end
    end
    
end

