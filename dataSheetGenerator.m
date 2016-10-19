classdef dataSheetGenerator
    %DATASHEETGENERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        jointModel
        templateDir
        outputDir
        
    end
    
    methods
        function this = dataSheetGenerator(jointModel)
            
            this.jointModel = jointModel;
            
            tmpStr = which('dataSheetGenerator');
            dsgRoot = fileparts(tmpStr);
            this.templateDir = [dsgRoot, filesep, 'templates' ];
            
            this.outputDir = ['.', filesep];
            
        end
       
        
      
%         function out = p_peakm(this)
%             % P_PEAKM Peak continous power (mechanical) [W]
%             %
%             %   p_peakm = gj.p_peakm
%             %
%             % Inputs:
%             %
%             % Outputs:
%             %   p_peakm: Value for the peak mechanicalpower obtained as the
%             %   product of peak speed and peak torque.
%             %   
%             %
%             % Notes::
%             %
%             %
%             % Examples::
%             %
%             %
%             % Author::
%             %  Joern Malzahn
%             %  Wesley Roozing
%             %
%             % See also i_c, i_p, genericJoint, jointBuilder.
%             
%             pStruct = this.jointModel; % Shorthand
%             out = pStruct.dq_p * this.t_p;
%         end
            
                
        function h = draw_speed_torque_curve(this)
            % DRAW_SPEED_TORQUE_CURVE Displays speed-torque-curve in a
            % figure.
            %
            %   fHandle = dsg.draw_speed_torque_curve
            %
            % Inputs:
            %
            % Outputs:
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
            % See also createDataSheet, genericJoint, jointBuilder.
            
            % Shorthands
            
            d_cm = this.jointModel.d_cm;
            d_cg = this.jointModel.d_cg;
            d_cb = this.jointModel.d_cb;
            d_m  = this.jointModel.d_m;
            d_g  = this.jointModel.d_g;
            d_b  = this.jointModel.d_b;
            slope = this.jointModel.dq_over_dm;
            dq_0 = this.jointModel.dq_0;
            t_stall = this.jointModel.t_stall;
            t_c = this.jointModel.t_c;
            dq_r = this.jointModel.dq_r;
            dq_c = this.jointModel.dq_c;
            p_cm = this.jointModel.p_cm;
            
            % Plotting options
            nVals = 100;
            xmax = 1.5 * t_c;
            ymax = 1.02 * dq_c;
            h = figure;
            hold on
                       
            % torque speed line
            
            mVals = (0:1/100:1) * t_stall;
            linCurve = dq_0 - slope * mVals;
            plot(mVals, linCurve, 'k', 'DisplayName', 'Torque-Speed Line')
            
            % Nominal operating point
            plot(t_c, dq_r, 'ko', 'DisplayName', 'Nominal Operating Point')
            
            % Friction
            speedVals = (0:1/nVals:1) * dq_c;
            Mc = d_cm + d_cg + d_cb;            % Static friction
            Mv = (d_m + d_g + d_b) * speedVals; % Velocity dependent friction
            Mf = Mc + Mv;

            fill([Mf, 0, 0], [speedVals speedVals(end) 0],0.8* [1 1 1],'LineStyle','none')
            alpha(0.25)
            plot(Mf, speedVals, 'k:', 'DisplayName', 'Friction Torque')
            
            % Plot limits
            speedVals = p_cm ./ mVals;
            plot(mVals,speedVals, 'b.', 'DisplayName', 'Rated Mechanical Power')
%             speedVals = this.p_peakm ./ mVals;
%             plot(mVals,speedVals, 'r--', 'DisplayName', 'Peak Mechanical Power')
            plot([0,xmax], dq_c * [1,1], 'r--', 'DisplayName', 'Maximum Continous Speed')
            plot(t_c*[1,1], [0,ymax], 'r--', 'DisplayName', 'Maximum Continous Torque')
            
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

            jM = this.jointModel; % Shorthand
            
            fid = fopen(cfgFName,'w+');
            
            fprintf(fid,'%s\n','% Mechanical Properties');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\gearratio{%d:1}\n', jM.n);
            fprintf(fid,'\\def \\stiffness{%5d}\n', jM.k_b);
            fprintf(fid,'\\def \\mass{%6.4f}\n',0);
            fprintf(fid,'\\def \\diameter{%6.4f}\n',0);
            fprintf(fid,'\\def \\actlength{%6.4f}\n',0);
            fprintf(fid,'\\def \\tmech{%4.2f}\n',0);
            fprintf(fid,'\\def \\inertiarotor{%6.4f}\n',jM.I_m);
            fprintf(fid,'\\def \\inertiagear{%6.4f}\n',jM.I_g);
            fprintf(fid,'\\def \\inertiaspring{%6.4f}\n',jM.I_b);
            fprintf(fid,'\\def \\viscousdamping{%6.4f}\n',jM.d_m+jM.d_g+jM.d_b);
            fprintf(fid,'\\def \\coulombdamping{%6.4f}\n',jM.d_cm+jM.d_cg+jM.d_cb);
            fprintf(fid,'\\def \\stribeckdamping{%6.4f}\n',jM.d_s);
            fprintf(fid,'\\def \\stribeckspeed{%6.4f}\n',jM.v_s);
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Electrical Properties');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\armaturresistance{%6.4f}\n',jM.r);
            fprintf(fid,'\\def \\armatureinductance{%6.4f}\n',jM.x);
            fprintf(fid,'\\def \\torqueconstant{%6.4f}\n',jM.k_t);
            fprintf(fid,'\\def \\speedconstant{%6.4f}\n',jM.k_t);
            fprintf(fid,'\\def \\speedtorquegradient{%6.4f}\n',jM.dq_over_dm);
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Rated Operation');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\ratedvoltage{%4.2f}\n',jM.v_0);
            fprintf(fid,'\\def \\noloadspeed{%4.2f}\n',jM.dq_0);
            fprintf(fid,'\\def \\noloadcurrent{%4.2f}\n',jM.v_0);
            fprintf(fid,'\\def \\ratedtorque{%4.2f}\n',jM.v_0);
            fprintf(fid,'\\def \\ratedcurrent{%4.2f}\n',jM.v_0);
            fprintf(fid,'\\def \\stalltorque{%4.2f}\n',jM.v_0);
            fprintf(fid,'\\def \\startingcurrent{%4.2f}\n',jM.v_0);
            fprintf(fid,'\\def \\maxefficiency{%4.2f}\n',jM.v_0);
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Specifications');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\maxspeed{%4.2f}\n',jM.dq_c);
            fprintf(fid,'\\def \\maxcurrent{%4.2f}\n',jM.i_p);
            fprintf(fid,'\\def \\maxtorque{%4.2f}\n',jM.t_c);
            fprintf(fid,'\\def \\restherm1{%4.2f}\n',jM.r_th1);
            fprintf(fid,'\\def \\restherm2{%4.2f}\n',jM.r_th2);
            fprintf(fid,'\\def \\Tthw{%4.2f}\n',jM.T_thw);
            fprintf(fid,'\\def \\Tthm{%4.2f}\n',jM.T_thm);
            fprintf(fid,'\\def \\TmpWindMax{%4.2f}\n',jM.Tmp_WMax);
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Power Rating');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\contpower{%4.2f}\n',jM.p_cm);
%            fprintf(fid,'\\def \\peakpower{%4.2f}\n', this.p_peakm);
            fprintf(fid,'\\def \\contcurrent{%4.2f}\n',jM.i_c);
            fprintf(fid,'\\def \\conttorque{%4.2f}\n',jM.t_c);
            
            fclose(fid);
        end
        
        function flag = compileTexFile(this,texFName)
            % Tested with MiKTeX/2.9
            cmd = ['lualatex.exe -synctex=-1 -interaction=nonstopmode ', texFName];
            flag = system(cmd);
                        
        end
        
    end
    
end

