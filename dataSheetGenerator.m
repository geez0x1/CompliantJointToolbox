classdef dataSheetGenerator
    %DATASHEETGENERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        jointModel
        templateDir
        outputDir
        
        texFName = 'cjtdsheet.tex';
        clsFName = 'cjtdsheet.cls';
        cfgFName = 'cjtdsheet.cfg';
        torqueSpeedFName = 'torqueSpeedCurve.pdf';
        
        outFName = 'DATASHEET.pdf';
    end
    
    methods
        function this = dataSheetGenerator(jointModel)
            
            this.jointModel = jointModel;
            
            tmpStr = which('dataSheetGenerator');
            dsgRoot = fileparts(tmpStr);
            this.templateDir = [dsgRoot, filesep, 'templates' ];
            
            this.outputDir = ['.', filesep];
            
        end
       
            
                
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
            t_r = this.jointModel.t_r;
            t_p = this.jointModel.t_p;
            dq_r = this.jointModel.dq_r;
            dq_p = this.jointModel.dq_p;
            p_cm = this.jointModel.p_cm;
            p_pm = this.jointModel.p_pm;
            t_NL = this.jointModel.t_NL;
            dq_NL = this.jointModel.dq_NL;
            
            % Plotting options
            nVals = 500;
            xmax = 1.2 * t_p;
            ymax = 1.02 * dq_p;
            h = figure;
            hold on
                       
            % torque speed line
            mVals = (0:1/nVals:1) * t_stall;
            linCurve = dq_0 - slope * mVals;
            plot(mVals, linCurve, 'k', 'DisplayName', 'Torque-Speed Line')
            
            % Rated operating point
            plot(t_r, dq_r, 'ko', 'DisplayName', 'Nominal Operating Point')
            
            % No-Load operating point
            plot(t_NL, dq_NL, 'ko', 'DisplayName', 'No-Load Operating Point')
            
           
            % Friction
            speedVals = (0:1/nVals:1) * dq_p;
            Mc = d_cm + d_cg + d_cb;            % Static friction
            Mv = (d_m + d_g + d_b) * speedVals; % Velocity dependent friction
            Mf = Mc + Mv;

            % Continuous operating range
            tOp =  [0, t_r,  t_r, 0 ].';
            dqOp = [0,   0, dq_r, dq_0].';
            fill(tOp,dqOp,0.9* [0.2 0.2 1],'LineStyle','none')
            
            fill([Mf, 0, 0], [speedVals speedVals(end) 0],0.8* [1 0 0],'LineStyle','none')
            alpha(0.25)
            plot(Mf, speedVals, 'k:', 'DisplayName', 'Friction Torque')
            
            % Peak Operation
            plot([t_p, t_p],[0, ymax],'r--', 'DisplayName', 'Peak Torque')
            speedVals = p_pm ./ mVals;
            plot(mVals,speedVals, 'r-', 'DisplayName', 'Peak Mechanical Power')
            plot([0,xmax], dq_p * [1,1], 'r--', 'DisplayName', 'Peak Speed')
            
            % Plot limits
            speedVals = p_cm ./ mVals;
            plot(mVals,speedVals, 'b-', 'DisplayName', 'Rated Mechanical Power')
%             speedVals = this.p_peakm ./ mVals;
%             plot(mVals,speedVals, 'r--', 'DisplayName', 'Peak Mechanical Power')
            plot([0,xmax], dq_r * [1,1], 'b--', 'DisplayName', 'Maximum Continous Speed')
            plot(t_r*[1,1], [0,ymax], 'b--', 'DisplayName', 'Maximum Continous Torque')
            
            % Annotations and Figure Style
            xlim([0,xmax]);
            ylim([0,ymax]);
            xlabel('torque [Nm]')
            ylabel('speed [rad/s]')
            box on
            
            %legend show;
        end
                

        
        function createDataSheet(this)
            
            copyfile([this.templateDir, filesep, this.texFName],this.texFName);
            copyfile([this.templateDir, filesep, this.clsFName],this.clsFName);
            
            this.createDefFile;
            
            this.makeDataSheetPlots;
            
            this.compileTexFile
            
            [~, fName] = fileparts(this.texFName);
            
            if ~exist(this.outputDir,'dir')
                mkdir(this.outputDir)
            end
            copyfile([fName,'.pdf'],[this.outputDir,filesep, this.outFName]);

            delete([fName,'.*'])
            delete([this.torqueSpeedFName,'.*'])
            
        end
        
        function makeDataSheetPlots(this)
            h = this.draw_speed_torque_curve;
                        
            set(gcf,'Units','centimeters');
            set(gcf,'PaperUnits','centimeters');
            pos = get(gcf,'Position');
            pos(3) = 18;
            pos(4) = 8;
             
            set(gcf,'Position',pos)
            set(h,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])          
            set(gca,'LooseInset',get(gca,'TightInset'))
            
            printpdf(gcf,this.torqueSpeedFName,'-r600')

             
        end
        
        function createDefFile(this)

            % The purpose of this function is to define a list of macros 
            % that can be called inside the Tex document to fill in the 
            % actuator parameters. The reason to do it in this way is to 
            % have single interface between the Matlab and the Tex world. 
            % 
            % It is good practice to define the actual symbols for
            % variables and properties as Tex commands or macros as well. 
            % If the a decision to change the symbol for a certain property 
            % quantity or quantity comes up, changes are required in just 
            % one place. 
            % As a convention, the macros representing the VALUE of the
            % property or quantity have 'val' as prefix. The macros 
            % describing the SYMBOL used in formulae have 'sym' as prefix.
            %
            % Both macro definitions go here in order to have them in just 
            % one place.
            
            jM = this.jointModel; % Shorthand
            
            fid = fopen(this.cfgFName,'w+');
            

            fprintf(fid,'%s\n','% Mechanical Properties');
            fprintf(fid,'%s\n','%');
            %Dimensions
            fprintf(fid,'\\def \\valDiameter{%4.1f}\n'             , jM.diam                   );
            fprintf(fid,'\\def \\symDiameter{%s}\n'                , 'D'                     );
            %
            fprintf(fid,'\\def \\valActlength{%4.1f}\n'            , jM.len                    );
            fprintf(fid,'\\def \\symActlength{%s}\n'            , 'l'                    );
            %
            % Inertae
            %
            fprintf(fid,'\\def \\valMass{%6.4f}\n'                 , jM.m                      );
            fprintf(fid,'\\def \\symMass{%s}\n'                 , 'm'                      );
            %
            fprintf(fid,'\\def \\valInertiarotor{%6.4f}\n'         , jM.I_m                    );
            fprintf(fid,'\\def \\symInertiarotor{%s}\n'         , 'I_m'                    );
            %
            fprintf(fid,'\\def \\valInertiagear{%6.4f}\n'          , jM.I_g                    );
            fprintf(fid,'\\def \\symInertiagear{%s}\n'          , 'I_g'                    );
            %
            fprintf(fid,'\\def \\valInertiaspring{%6.4f}\n'        , jM.I_b                    );
            fprintf(fid,'\\def \\symInertiaspring{%s}\n'        , 'I_m'                    );
            %
            fprintf(fid,'\\def \\valTmech{%6.4f}\n'                , jM.T_mech                 );
            fprintf(fid,'\\def \\symTmech{%s}\n'                , 'T_{mech}'                 );
            %
            % Stiffnesses
            %
            fprintf(fid,'\\def \\valSpringstiffness{%5d}\n'        , jM.k_b                    );
            fprintf(fid,'\\def \\symSpringstiffness{%s}\n'        , 'k_b'                    );
            fprintf(fid,'\\def \\valGearstiffness{%5d}\n'          , jM.k_g                    );
            fprintf(fid,'\\def \\symGearstiffness{%s}\n'          , 'k_g'                    );
            %
            % Friction
            %
            fprintf(fid,'\\def \\valViscousdamping{%6.4f}\n'       , jM.d_m+jM.d_g+jM.d_b      );
            fprintf(fid,'\\def \\symViscousdamping{%s}\n'       , 'd_v'      );
            %
            fprintf(fid,'\\def \\valCoulombdamping{%6.4f}\n'       , jM.d_cm+jM.d_cg+jM.d_cb   );
            fprintf(fid,'\\def \\symCoulombdamping{%s}\n'       , 'd_c'   );
            %
            fprintf(fid,'\\def \\valStribeckdamping{%6.4f}\n'      , jM.d_s                    );
            fprintf(fid,'\\def \\symStribeckdamping{%s}\n'      , 'd_s'                    );
            %
            fprintf(fid,'\\def \\valStribeckspeed{%6.4f}\n'        , jM.v_s                    );
            fprintf(fid,'\\def \\symStribeckspeed{%s}\n'        , '\dot{\theta}_s'                    );
            %
            % Electrical
            %
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Electrical Properties');
            fprintf(fid,'%s\n','%');
            %
            fprintf(fid,'\\def \\valTorqueconstant{%6.4f}\n'       , jM.k_t                    );
            fprintf(fid,'\\def \\symTorqueconstant{%s}\n'       , 'k_\tau'                    );
            %
            fprintf(fid,'\\def \\valGeneratorconstant{%6.4f}\n'        , jM.k_w                    );
            fprintf(fid,'\\def \\symGeneratorconstant{%s}\n'        , 'k_\omega'                    );
            %
            fprintf(fid,'\\def \\valArmatureresistance{%6.4f}\n'    , jM.r                      );
            fprintf(fid,'\\def \\symArmatureresistance{%s}\n'    , 'r_A'                      );
            %
            fprintf(fid,'\\def \\valArmatureinductance{%6.4f}\n'   , jM.x                      );
            fprintf(fid,'\\def \\symArmatureinductance{%s}\n'   , 'x_A'                      );
            %
            fprintf(fid,'\\def \\valTel{%6.4f}\n'                  , jM.T_el                   );
            fprintf(fid,'\\def \\symTel{%s}\n'                  , 'T_{el}'                   );
            %
            % Thermal
            %
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Thermal Properties');
            fprintf(fid,'%s\n','%');
            %
            fprintf(fid,'\\def \\valResthermWH{%4.2f}\n'           , jM.r_th1                  );
            fprintf(fid,'\\def \\symResthermWH{%s}\n'           , 'r_{th1}'                  );
            %
            fprintf(fid,'\\def \\valResthermHA{%4.2f}\n'           , jM.r_th2                  );
            fprintf(fid,'\\def \\symResthermHA{%s}\n'           , 'r_{th2}'                  );
            %
            fprintf(fid,'\\def \\valTthw{%4.2f}\n'                 , jM.T_thw                  );
            fprintf(fid,'\\def \\symTthw{%s}\n'                 , 'T_{th,\:w}'                  );
            %
            fprintf(fid,'\\def \\valTthm{%4.2f}\n'                 , jM.T_thm                  );
            fprintf(fid,'\\def \\symTthm{%s}\n'                 , 'T_{th,\:m}'                  );
            %
            fprintf(fid,'\\def \\valTmpWindMax{%4.2f}\n'           , jM.Tmp_WMax               );
            fprintf(fid,'\\def \\symTmpWindMax{%s}\n'           , '\nu_{w,\:Max}'               );
            %
            fprintf(fid,'\\def \\valTmpANom{%4.2f}\n'              , jM.Tmp_ANom               );
            fprintf(fid,'\\def \\symTmpANom{%s}\n'              , '\nu_{A,\:Nom}'               );
            %
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Rated Operation');
            fprintf(fid,'%s\n','%');
            %
            fprintf(fid,'\\def \\valRatedvoltage{%4.2f}\n'         , jM.v_0                    );
            fprintf(fid,'\\def \\symRatedvoltage{%s}\n'         , 'V_0'                    );
            %
            fprintf(fid,'\\def \\valRatedcurrent{%4.2f}\n'         , jM.i_c                    );
            fprintf(fid,'\\def \\symRatedcurrent{%s}\n'         , 'i_R'                    );
            %
            fprintf(fid,'\\def \\valRatedtorque{%4.2f}\n'          , jM.t_r                    );
            fprintf(fid,'\\def \\symRatedtorque{%s}\n'          , '\tau_R'                    );
            %
            fprintf(fid,'\\def \\valRatedspeed{%4.2f}\n'           , jM.dq_r                   );
            fprintf(fid,'\\def \\symRatedspeed{%s}\n'           , '\dot{\theta}_R'                   );
            %
            fprintf(fid,'\\def \\valRatedpowere{%4.2f}\n'          , jM.p_ce                   );
            fprintf(fid,'\\def \\symRatedpowere{%s}\n'          , 'P_{R,\:E}'                   );
            %
            fprintf(fid,'\\def \\valRatedpowerm{%4.2f}\n'          , jM.p_cm                   );
            fprintf(fid,'\\def \\symRatedpowerm{%s}\n'          , 'P_{R,\:M}'                   );
            %
            fprintf(fid,'\\def \\valNoloadcurrent{%4.2f}\n'        , jM.i_NL                   );
            fprintf(fid,'\\def \\symNoloadcurrent{%s}\n'        , 'i_{NL}'                   );
            %
            fprintf(fid,'\\def \\valNoloadtorque{%4.2f}\n'         , jM.t_NL                   );
            fprintf(fid,'\\def \\symNoloadtorque{%s}\n'         , '\tau_{NL}'                   );
            %
            fprintf(fid,'\\def \\valNoloadspeed{%4.2f}\n'          , jM.dq_NL                   );
            fprintf(fid,'\\def \\symNoloadspeed{%s}\n'          , '\dot{\theta}_{NL}'                  );
            %
            fprintf(fid,'\\def \\valStalltorque{%4.2f}\n'          , jM.t_stall                );
            fprintf(fid,'\\def \\symStalltorque{%s}\n'          , '\tau_{Stall}'                );
            %
            fprintf(fid,'\\def \\valStartingcurrent{%4.2f}\n'      , jM.i_start                );
            fprintf(fid,'\\def \\symStartingcurrent{%s}\n'      , 'i_{start}'                );
            %
            fprintf(fid,'\\def \\valSpeedtorquegradient{%6.4f}\n'  , jM.dq_over_dm             );            
            fprintf(fid,'\\def \\symSpeedtorquegradient{%s}\n'  , '\text{d}\dot{\theta}/\text{d}\tau'             );
            %
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Peak Operation');
            fprintf(fid,'%s\n','%');
            %
            fprintf(fid,'\\def \\valMaxcurrent{%4.2f}\n'           , jM.i_p                    );
            fprintf(fid,'\\def \\symMaxcurrent{%s}\n'           , 'i_p'                    );
            %
            fprintf(fid,'\\def \\valMaxtorque{%4.2f}\n'            , jM.t_p                    );
            fprintf(fid,'\\def \\symMaxtorque{%s}\n'            , '\tau_p'                    );
            %
            fprintf(fid,'\\def \\valMaxspeed{%4.2f}\n'             , jM.dq_p                   );
            fprintf(fid,'\\def \\symMaxspeed{%s}\n'             , '\dot{\theta}_p'                   );
            %
            fprintf(fid,'\\def \\valMaxpowere{%4.2f}\n'            , jM.p_pe                   );
            fprintf(fid,'\\def \\symMaxpowere{%s}\n'            , 'P_{p,\:E}'                   );
            %
            fprintf(fid,'\\def \\valMaxpowerm{%4.2f}\n'            , jM.p_pm                   );
            fprintf(fid,'\\def \\symMaxpowerm{%s}\n'            , 'P_{p,\:M}'                   );
            %
            fprintf(fid,'\\def \\valGearratio{%d:1}\n'                , jM.n                      );
            fprintf(fid,'\\def \\symGearratio{%s}\n'                , 'N'                      );
            %
            %            fprintf(fid,'\\def \\maxefficiency{%4.2f}\n',jM.v_0);
            %
            fclose(fid);
        end
        
        function flag = compileTexFile(this)
            % Tested with MiKTeX/2.9
            cmd = ['lualatex.exe -synctex=-1 -interaction=nonstopmode ', this.texFName];
            flag = system(cmd);
                        
        end
        
    end
    
end

