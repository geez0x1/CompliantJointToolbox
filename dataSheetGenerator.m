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
            t_NL = this.jointModel.t_NL;
            dq_NL = this.jointModel.dq_NL;
            
            % Plotting options
            nVals = 200;
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
            
            % No-Load operating point
            plot(t_NL, dq_NL, 'ko', 'DisplayName', 'No-Load Operating Point')
            
           
            % Friction
            speedVals = (0:1/nVals:1) * dq_c;
            Mc = d_cm + d_cg + d_cb;            % Static friction
            Mv = (d_m + d_g + d_b) * speedVals; % Velocity dependent friction
            Mf = Mc + Mv;

            % Continuous operating range
            tOp =  [0, t_c,  t_c, t_NL ].';
            dqOp = [0,   0, dq_r, dq_NL].';
            fill(tOp,dqOp,0.8* [0 1 0],'LineStyle','none')
            
            fill([Mf, 0, 0], [speedVals speedVals(end) 0],0.8* [1 0 0],'LineStyle','none')
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
            xlim([0,xmax]);
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
            pos(3) = 18;
            pos(4) = 8;
             
            set(gcf,'Position',pos)
            set(h,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
            
             print(gcf,fName,'-dpdf','-r600')
             
        end
        
        function createDefFile(this,cfgFName)

            jM = this.jointModel; % Shorthand
            
            fid = fopen(cfgFName,'w+');
            
            fprintf(fid,'%s\n','% Mechanical Properties');
            fprintf(fid,'%s\n','%');
            %Dimensions
            fprintf(fid,'\\def \\diameter{%4.1f}\n'             , jM.diam                   );
            fprintf(fid,'\\def \\actlength{%4.1f}\n'            , jM.len                    );
            % Inertae
            fprintf(fid,'\\def \\mass{%6.4f}\n'                 , jM.m                      );
            fprintf(fid,'\\def \\inertiarotor{%6.4f}\n'         , jM.I_m                    );
            fprintf(fid,'\\def \\inertiagear{%6.4f}\n'          , jM.I_g                    );
            fprintf(fid,'\\def \\inertiaspring{%6.4f}\n'        , jM.I_b                    );
            fprintf(fid,'\\def \\Tmech{%6.4f}\n'                , jM.T_mech                 );
            % Stiffnesses
            fprintf(fid,'\\def \\springstiffness{%5d}\n'        , jM.k_b                    );
            fprintf(fid,'\\def \\gearstiffness{%5d}\n'          , jM.k_g                    );
            % Friction
            fprintf(fid,'\\def \\viscousdamping{%6.4f}\n'       , jM.d_m+jM.d_g+jM.d_b      );
            fprintf(fid,'\\def \\coulombdamping{%6.4f}\n'       , jM.d_cm+jM.d_cg+jM.d_cb   );
            fprintf(fid,'\\def \\stribeckdamping{%6.4f}\n'      , jM.d_s                    );
            fprintf(fid,'\\def \\stribeckspeed{%6.4f}\n'        , jM.v_s                    );
            % Electrical
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Electrical Properties');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\torqueconstant{%6.4f}\n'       , jM.k_t                    );
            fprintf(fid,'\\def \\speedconstant{%6.4f}\n'        , jM.k_w                    );
            fprintf(fid,'\\def \\armaturresistance{%6.4f}\n'    , jM.r                      );
            fprintf(fid,'\\def \\armatureinductance{%6.4f}\n'   , jM.x                      );
            fprintf(fid,'\\def \\Tel{%6.4f}\n'                  , jM.T_el                   );
            % Thermal
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Thermal Properties');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\resthermWH{%4.2f}\n'           , jM.r_th1                  );
            fprintf(fid,'\\def \\resthermHA{%4.2f}\n'           , jM.r_th2                  );
            fprintf(fid,'\\def \\Tthw{%4.2f}\n'                 , jM.T_thw                  );
            fprintf(fid,'\\def \\Tthm{%4.2f}\n'                 , jM.T_thm                  );
            fprintf(fid,'\\def \\TmpWindMax{%4.2f}\n'           , jM.Tmp_WMax               );
            fprintf(fid,'\\def \\TmpANom{%4.2f}\n'              , jM.Tmp_ANom               );


            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Rated Operation');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\ratedvoltage{%4.2f}\n'         , jM.v_0                    );
            fprintf(fid,'\\def \\ratedcurrent{%4.2f}\n'         , jM.i_c                    );
            fprintf(fid,'\\def \\ratedtorque{%4.2f}\n'          , jM.t_c                    );
            fprintf(fid,'\\def \\ratedspeed{%4.2f}\n'           , jM.dq_r                   );
            fprintf(fid,'\\def \\ratedpowere{%4.2f}\n'          , jM.p_ce                   );
            fprintf(fid,'\\def \\ratedpowerm{%4.2f}\n'          , jM.p_cm                   );
            fprintf(fid,'\\def \\noloadcurrent{%4.2f}\n'        , jM.i_NL                   );
            fprintf(fid,'\\def \\noloadtorque{%4.2f}\n'         , jM.t_NL                   );
            fprintf(fid,'\\def \\noloadspeed{%4.2f}\n'          , jM.dq_0                   );
            fprintf(fid,'\\def \\stalltorque{%4.2f}\n'          , jM.t_stall                );
            fprintf(fid,'\\def \\startingcurrent{%4.2f}\n'      , jM.i_start                );
            fprintf(fid,'\\def \\speedtorquegradient{%6.4f}\n'  , jM.dq_over_dm             );            

            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Peak Operation');
            fprintf(fid,'%s\n','%');
            fprintf(fid,'\\def \\maxcurrent{%4.2f}\n'           , jM.i_p                    );
            fprintf(fid,'\\def \\maxtorque{%4.2f}\n'            , jM.t_p                    );
%           fprintf(fid,'\\def \\maxspeed{%4.2f}\n'             , jM.dq_p                   );
            fprintf(fid,'\\def \\maxspeed{%4.2f}\n'             , -1                        );
            fprintf(fid,'\\def \\maxpowere{%4.2f}\n'            , jM.p_pe                   );
            fprintf(fid,'\\def \\maxpowerm{%4.2f}\n'            , jM.p_pm                   );
            fprintf(fid,'\\def \\contcurrent{%4.2f}\n'          , jM.i_c                    );
            fprintf(fid,'\\def \\conttorque{%4.2f}\n'           , jM.t_c                    );
            fprintf(fid,'\\def \\gearratio{%d:1}\n'             , jM.n                      );
            %            fprintf(fid,'\\def \\maxefficiency{%4.2f}\n',jM.v_0);
            fclose(fid);
        end
        
        function flag = compileTexFile(this,texFName)
            % Tested with MiKTeX/2.9
            cmd = ['lualatex.exe -synctex=-1 -interaction=nonstopmode ', texFName];
            flag = system(cmd);
                        
        end
        
    end
    
end

