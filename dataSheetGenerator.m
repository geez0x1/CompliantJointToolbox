classdef dataSheetGenerator
    %DATASHEETGENERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        paramStruct
        templateDir
        outputDir
        
        a_CU = 0.0039;  % Resistance coefficient of copper [1/K]
        c_CU = 380;     % Specific thermal capacitance of copper [J/kg/K]
        c_FE = 460;     % Specific thermal capacitance of iron [J/kg/K]
    end
    
    methods
        function this = dataSheetGenerator(paramStruct)
            
            this.paramStruct = paramStruct;
            
            tmpStr = which('dataSheetGenerator');
            dsgRoot = fileparts(tmpStr);
            this.templateDir = [dsgRoot, filesep, 'templates' ];
            
            this.outputDir = ['.', filesep];
            
        end
        
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
            % See also t_c, p_rce, genericJoint, jointBuilder.
            
            % Shorthands
            aCU   = this.a_CU;                      % Resistance coefficient for copper
            r_TA  = this.paramStruct.r;             % Winding resistance at normal ambient temperature
            T_A   = this.paramStruct.Tmp_ANom;      % Normal ambient temperature
            T_W   = this.paramStruct.Tmp_WMax;      % Maximum allowed winding temperature
            r_th1 = this.paramStruct.r_th1;         % Thermal resistance Winding-Housing
            r_th2 = this.paramStruct.r_th2;         % Thermal resistance Housing-Air
            
            dT = T_W-T_A;                           % Allowed temperature rise
            
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
            iSquare = dT / (r_th1 + r_th2) / r_TA / (1 + aCU*dT);
            
            % The maximum continuous current thus computes to:
            out = sqrt(iSquare);
            
        end
        
        function out = t_c(this)
            % T_C Continuous torque [Nm]
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
            
            % Shorthands
            k_t = this.paramStruct.k_t;     % Torque constant
            i_c = this.i_c;                 % Maximum continuous current
            n   = this.paramStruct.n;       % Transmission ratio
            
            % Compute maximum continuous motor torque
            out = k_t * i_c * n;
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
            
            % Shorthands
            pStruct = this.paramStruct; % Shorthand
            
            d_cm = this.paramStruct.d_cm;
            d_cg = this.paramStruct.d_cg;
            d_cb = this.paramStruct.d_cb;
            d_m  = this.paramStruct.d_m;
            d_g  = this.paramStruct.d_g;
            d_b  = this.paramStruct.d_b;
            
            % Plotting options
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
            plot(this.t_c, this.dq_r, 'ko', 'DisplayName', 'Nominal Operating Point')
            
            % Friction
            speedVals = (0:1/nVals:1) * pStruct.dq_c;
            Mc = d_cm + d_cg + d_cb;            % Static friction
            Mv = (d_m + d_g + d_b) * speedVals; % Velocity dependent friction
            Mf = Mc + Mv;

            fill([Mf, 0, 0], [speedVals speedVals(end) 0],0.8* [1 1 1],'LineStyle','none')
            alpha(0.25)
            plot(Mf, speedVals, 'k:', 'DisplayName', 'Friction Torque')
            
            % Plot limits
            speedVals = this.p_rcm ./ mVals;
            plot(mVals,speedVals, 'b.', 'DisplayName', 'Rated Mechanical Power')
            speedVals = this.p_peakm ./ mVals;
            plot(mVals,speedVals, 'r--', 'DisplayName', 'Peak Mechanical Power')
            plot([0,xmax], pStruct.dq_c * [1,1], 'r--', 'DisplayName', 'Maximum Continous Speed')
            plot(this.t_c*[1,1], [0,ymax], 'r--', 'DisplayName', 'Maximum Continous Torque')
            
            % Annotations and Figure Style
            xlim([-1,xmax]);
            ylim([0,ymax]);
            xlabel('torque [Nm]')
            ylabel('speed [rad/s]')
            box on
            
            %legend show;
        end
                

        
        function createDataSheet(this)
            texFName = 'cjtdsheet.tex';
            clsFName = 'cjtdsheet.cls';
            cfgFName = 'cjtdsheet.cfg';
            figFName = 'dummy';
            
            copyfile([this.templateDir, filesep, texFName],texFName);
            copyfile([this.templateDir, filesep, clsFName],clsFName);
            
            this.createDefFile(cfgFName);
            
            this.makeDataSheetPlots(figFName);
            
            this.compileTexFile(texFName)
            
            [~, fName] = fileparts(texFName);
            
            if ~exist(this.outputDir,'dir')
                mkdir(this.outputDir)
            end
            copyfile([fName,'.pdf'],[this.outputDir,filesep,'DATASHEET.pdf']);

            delete([fName,'.*'])
            delete([figFName,'.*'])
            
        end
        
        function makeDataSheetPlots(this,fName)
            h = this.draw_speed_torque_curve;
                        
            set(gcf,'Units','centimeters');
            set(gcf,'PaperUnits','centimeters');
            pos = get(gcf,'Position');
            pos(3) = 16;
            pos(4) = 9;
             
            set(gcf,'Position',pos)
            set(h,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
            
             print(gcf,fName,'-dpdf','-r600')
            
        end
        
        function createDefFile(this,cfgFName)

            pStruct = this.paramStruct; % Shorthand
            
            fid = fopen(cfgFName,'w+');
            
            fprintf(fid,'%s\n','% Mechanical Properties');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\gearratio{%d:1}\n', pStruct.n);
            fprintf(fid,'\\def \\stiffness{%5d}\n', pStruct.k_b);
            fprintf(fid,'\\def \\mass{%6.4f}\n',0);
            fprintf(fid,'\\def \\diameter{%6.4f}\n',0);
            fprintf(fid,'\\def \\actlength{%6.4f}\n',0);
            fprintf(fid,'\\def \\tmech{%4.2f}\n',0);
            fprintf(fid,'\\def \\inertiarotor{%6.4f}\n',pStruct.I_m);
            fprintf(fid,'\\def \\inertiagear{%6.4f}\n',pStruct.I_g);
            fprintf(fid,'\\def \\inertiaspring{%6.4f}\n',pStruct.I_b);
            fprintf(fid,'\\def \\viscousdamping{%6.4f}\n',pStruct.d_m+pStruct.d_g+pStruct.d_b);
            fprintf(fid,'\\def \\coulombdamping{%6.4f}\n',pStruct.d_cm+pStruct.d_cg+pStruct.d_cb);
            fprintf(fid,'\\def \\stribeckdamping{%6.4f}\n',pStruct.d_s);
            fprintf(fid,'\\def \\stribeckspeed{%6.4f}\n',pStruct.v_s);
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Electrical Properties');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\armaturresistance{%6.4f}\n',pStruct.r);
            fprintf(fid,'\\def \\armatureinductance{%6.4f}\n',pStruct.x);
            fprintf(fid,'\\def \\torqueconstant{%6.4f}\n',pStruct.k_t);
            fprintf(fid,'\\def \\speedconstant{%6.4f}\n',pStruct.k_t);
            fprintf(fid,'\\def \\speedtorquegradient{%6.4f}\n',this.dq_over_dm);
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Rated Operation');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\ratedvoltage{%4.2f}\n',pStruct.v_0);
            fprintf(fid,'\\def \\noloadspeed{%4.2f}\n',this.dq_0);
            fprintf(fid,'\\def \\noloadcurrent{%4.2f}\n',pStruct.v_0);
            fprintf(fid,'\\def \\ratedtorque{%4.2f}\n',pStruct.v_0);
            fprintf(fid,'\\def \\ratedcurrent{%4.2f}\n',pStruct.v_0);
            fprintf(fid,'\\def \\stalltorque{%4.2f}\n',pStruct.v_0);
            fprintf(fid,'\\def \\startingcurrent{%4.2f}\n',pStruct.v_0);
            fprintf(fid,'\\def \\maxefficiency{%4.2f}\n',pStruct.v_0);
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Specifications');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\maxspeed{%4.2f}\n',pStruct.dq_c);
            fprintf(fid,'\\def \\maxcurrent{%4.2f}\n',pStruct.i_p);
            fprintf(fid,'\\def \\maxtorque{%4.2f}\n',this.t_c);
            fprintf(fid,'\\def \\restherm1{%4.2f}\n',pStruct.r_th1);
            fprintf(fid,'\\def \\restherm2{%4.2f}\n',pStruct.r_th2);
            fprintf(fid,'\\def \\Tthw{%4.2f}\n',pStruct.T_thw);
            fprintf(fid,'\\def \\Tthm{%4.2f}\n',pStruct.T_thm);
            fprintf(fid,'\\def \\TmpWindMax{%4.2f}\n',pStruct.Tmp_WMax);
            fprintf(fid,'\\def \\TmpAmbMax{%4.2f}\n',pStruct.Tmp_AMax);
            fprintf(fid,'\\def \\TmpAmbMin{%4.2f}\n',pStruct.Tmp_AMin);
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Power Rating');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\contpower{%4.2f}\n',this.p_rcm);
            fprintf(fid,'\\def \\peakpower{%4.2f}\n', this.p_peakm);
            fprintf(fid,'\\def \\contcurrent{%4.2f}\n',this.i_c);
            fprintf(fid,'\\def \\conttorque{%4.2f}\n',this.t_c);
            
            fclose(fid);
        end
        
        function flag = compileTexFile(this,texFName)
            % Tested with MiKTeX/2.9
            cmd = ['lualatex.exe -synctex=-1 -interaction=nonstopmode ', texFName];
            flag = system(cmd);
                        
        end
        
    end
    
end

