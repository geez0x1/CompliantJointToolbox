classdef dataSheetGenerator
%DATASHEETGENERATOR Module that creates datasheets for actuator models
%created with the CompliantJointToolbox. 
%
% The generated datasheets comprise a table listing the parameters of
% the actuator including a description for each parameter. The
% datasheet also includes figures that display the actuator
% characteristics such as the torque-speed curve.
%
% Example:
%
%     % Instantiate a jointBuilder class object
%     jb = jointBuilder;
% 
%     % Use the joint builder to create a joint model. 
%     motorName = 'dummyMotor';          % In this example we instantiate a model based on the dummyMotor parameter set.
%     dynName = 'full';                  % We chose the full dynamics.
%     jb.buildJoint(motorName, dynName); % And build the model...
%     addpath(jb.buildDir)               % ... and add the build directory to the search path.
% 
%     % We instantiate first the model ...
%     aJoint = eval([motorName,'_',dynName]); 
% 
%     % ... and then the data sheet generator.
%     dsg = dataSheetGenerator(aJoint);       
% 
%     % Finally, we generate the data sheet:
%     dsg.createDataSheet;                    
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

    properties (SetAccess = private)
        jointModel                                     % The joint object to generate the datasheet for.
        templateDir                                    % Directory that contains all template files for the datasheet.
        outputDir                                      % Directory to place the generated files.
        
        texFName = 'cjtdsheet.tex';                    % Datasheet source filename.
        clsFName = 'cjtdsheet.cls';                    % Name of the style file for the datasheet.
        cfgFName = 'cjtdsheet.cfg';                    % Name of the temporary file containing the macro configuration for the datasheet sources.
        torqueSpeedFName = 'torqueSpeedCurve.pdf';     % Name of the temporary file that contains the torque speed curve of the actuator.
        
    end
    
    properties
        
        outFName = 'datasheet.pdf';                    % Name stub for the generated datasheet file. The final name will 
                                                       % be composed of the "jointModel.name_outFName".
        plotResolution = 600;                          % Print resolution for the figures contained in the datasheet.
        nPlotVals = 10000;                               % Number of samples used to produce the graphs in the datasheet figures.
        
    end
    
    methods
        function this = dataSheetGenerator(jointModel)
            % DATASHEETGENERATOR Constructor for the dataSheetGenerator
            % class. 
            %
            %   dsg = dataSheetGenerator(jointModel)
            %
            % Inputs:
            %   jointModel: CompliantJointToolbox joint model class object.
            %
            % Outputs:
            %   dsg: dataSheetGenerator class object.
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
            % See also createDataSheet, genericJoint, jointBuilder. 
            
            this.jointModel = jointModel;
            
            tmpStr = which('dataSheetGenerator');
            dsgRoot = fileparts(tmpStr);
            this.templateDir = [dsgRoot, filesep, 'templates' ];
            
            this.outputDir = ['.', filesep];
            
        end
        
        function speed = torque_speed_curve(this,tau)
            % TORQUE_SPEED_CURVE Computes the maximum feasible motor speed
            % corresponding for a given torque level.
            %
            %   speed = dsg.torque_speed_curve(tau)
            %
            % Inputs:
            %   tau: Deliviered actuator torque.
            %
            % Outputs:
            %   speed: Speed value corresponding to the torque level.
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
            % See also createDataSheet, genericJoint, jointBuilder.
            
            % shorthands
            slope = this.jointModel.dq_over_dm;
            dq_0 = this.jointModel.dq_0;
            
            % output
            speed = dq_0 - slope * tau;
        end
        
        
        function [h, hAx, hLine1, hLine2] = draw_efficiency_curve(this)
            % DRAW_EFFICIENCY_CURVE Creates a plot with two y-axes. One 
            % axis displays the actuator efficiency. The other one displays
            % the delivered mechanical power. Both plots are shown with 
            % respect to the generated torque.
            %
            %   [hFigure hAx,hLine1,hLine2] = dsg.draw_efficiency_curve
            %
            % Inputs:
            %
            % Outputs:
            %   hFigure: Handle to the figure containing the plots
            %   hAx: Handle to the axes forming the plot (one for each y axis)
            %   hLine1: Handle to the lines of the power plot
            %   hLine2: Handle to the lines of the efficiency plot
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
            % See also createDataSheet, genericJoint, jointBuilder.      
            
            % Shorthands
            t_p = this.jointModel.t_p;
            v_0 = this.jointModel.v_0;            
            k_t = this.jointModel.k_t;
            r   = this.jointModel.r;
            t_NL = this.jointModel.t_NL;
            t_stall = this.jointModel.t_stall;
            N = this.jointModel.n;
            
            tauVals = (0:1/this.nPlotVals:1) * t_p;
            
            PL = v_0 / k_t/N * tauVals - r / k_t^2/N^2 * tauVals.^2;
            P_tot = v_0 / k_t/N * tauVals;
            speedVals = this.torque_speed_curve(tauVals);
            
            eta = 100*( speedVals .* (tauVals - t_NL) * N * k_t / v_0 ) ./ tauVals;
            eta_max = 100 * ( 1 - sqrt( t_NL / t_stall)  );
            
            h = figure;
            [hAx,hLine1,hLine2] = plotyy(tauVals(:), eta(:), ... Efficiency plot 
                tauVals(:),[PL(:) P_tot(:)]);                % Power plot
            
            % Manipulate style of the efficiency plot 
            set(hLine1,'Color','r')
            ylabel(hAx(1), 'Efficiency [%]','Color','r');
            set(hAx(1),'ylim',[0,100])
            
            % Manipulate style of the power plot 
            set(hLine2,'Color','b')
            set(hLine2(2),'LineStyle','--')
            ylabel(hAx(2),'Load Power [W]','Color','b');
            
            set(hAx,'xlim',[0,t_p]);
            xlabel('torque [Nm]')
            title(['Maximum Efficiency_ ', sprintf('%2.0f',eta_max),' %']);
            
        end

        function [h, hAx, hLine1, hLine2] = draw_thermal_characteristics(this)
            % DRAW_THERMAL_CHARACTERISTICS Creates a plot with two y-axes. 
            %
            %   [hFigure hAx,hLine1,hLine2] = dsg.draw_thermal_characteristics
            %
            % Inputs:
            %
            % Outputs:
            %   hFigure: Handle to the figure containing the plots
            %   hAx: Handle to the axes forming the plot (one for each y axis)
            %   hLine1: Handle to the lines of the 
            %   hLine2: Handle to the lines of the 
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
            % See also createDataSheet, genericJoint, jointBuilder.      

            % Shorthands
            r_th1    = this.jointModel.r_th1;      % Thermal Resistance Windings to Housing [K/W]
            r_th2    = this.jointModel.r_th2;      % Thermal Resistance Housing to Air [K/W]
            r_TA     = this.jointModel.r;          % Armature resistance at normal temperature [Ohm]
            N        = this.jointModel.n;          % Transmission ratio [.]    
            k_t      = this.jointModel.k_t;        % Torque constant [Nm/A]
            aCU      = this.jointModel.a_CU;       % Resistance coefficient of copper [1/K]
            T_thw    = this.jointModel.T_thw;      % Thermal Time Constant of the Windings [s] 
            Tmp_WMax = this.jointModel.Tmp_WMax;   % Maximum Armature Temperature [�C]  
            Tmp_ANom = this.jointModel.Tmp_ANom;  % Normal Ambient Temperature [�C]
            t_p      = this.jointModel.t_p;
            t_r = this.jointModel.t_r;
            
            TmeltCU = 1084; % Melting temperature of copper
            
            tauVals = (0:1/this.nPlotVals:1) * t_p; 

            % Steady State Temperature
            resCoeff = (r_th1 + r_th2) * r_TA / N^2 / k_t^2;
            Tmp_W = Tmp_ANom + (  ( resCoeff * tauVals.^2 ) ./ ( 1 - aCU * resCoeff *  min(tauVals,t_r).^2 )  );
            
            idx = find(Tmp_W > TmeltCU,1,'first');
            Tmp_W(idx:end) = TmeltCU;

            TmpLimW  = Tmp_WMax*ones(size(Tmp_W(:)));
            TmpLimCU = TmeltCU*ones(size(Tmp_W(:)));
            
            % Time to reach max winding temperature from normal conditions
            Tmax = Tmp_WMax;
            T0 = Tmp_ANom;
            dTend = Tmp_W - T0;
            tCrit = - T_thw * log( - (Tmax - T0  - dTend ) ./ dTend);
            
            idx = find(tauVals > t_r,1,'first');
            tCrit(1:idx) = inf;
            
            
            h = figure;
            [hAx,hLine1,hLine2] = plotyy(tauVals(:), [Tmp_W(:), TmpLimW(:),TmpLimCU(:) ], ... Steady State Temperature
                tauVals(:),[tCrit]);                % Time to critical temperature
%             
             % Manipulate style of the efficiency plot 
             set(hLine1,'Color','r')
             set(hLine1(2),'LineStyle','--')
             set(hLine1(3),'LineStyle','--')
             ylabel(hAx(1), 'Steady State Temp. [^\circC]','Color','r');
            set(hAx(1),'ylim',[0,1.05*max((Tmp_W))]);
%             
             % Manipulate style of the power plot 
             set(hLine2,'Color','b')

             ylabel(hAx(2),'Time to Crit Temperature [s]','Color','b');
             set(hAx(2),'ylim',[0,1.05*max(real(tCrit))]);
%             
             set(hAx,'xlim',[0,t_p]);
%             xlabel('torque [Nm]')
%             title(['Maximum Efficiency_ ', sprintf('%2.0f',eta_max),' %']);
            
        end
        
        function h = draw_torque_frequency_curve(this)
            
            h = figure;
            hold on
            
            t_stall = this.jointModel.t_stall;
            k = this.jointModel.k_b;
            slope = this.jointModel.dq_over_dm;
            dq_0 = this.jointModel.dq_0;
            t_p = this.jointModel.t_p;
            t_r = this.jointModel.t_r;
            t_NL = this.jointModel.t_NL;
            dq_NL = this.jointModel.dq_NL;
            
            I = this.jointModel.I_m + this.jointModel.I_g;
            w0 = sqrt(k/I);
            f0 = w0/2/pi;

            xmax = 1.2 * t_p;   
            
            torque = (0:1/this.nPlotVals:1) * t_stall;
            peakSpeeds = this.computeMaxPeakSpeed(torque);
            contSpeeds = this.computeMaxContSpeed(torque);
            
%             contW = ( contSpeeds * sqrt(I*k) ) ./ sqrt( torque.^2 + I^2 * contSpeeds.^2 );
            contW = w0 * sqrt( ( contSpeeds * k * I ) ./ ( torque + contSpeeds * k * I ) );
            contF = contW / 2 / pi;
%             
%             peakW = ( peakSpeeds * sqrt(I*k) ) ./ sqrt( torque.^2 + I^2 * peakSpeeds.^2 );
            peakW = w0 * sqrt( ( peakSpeeds * k * I ) ./ ( torque + peakSpeeds * k * I ) );
            peakF = peakW / 2 / pi;
            

            plot(torque,contF,'b--');
            plot(torque,peakF,'r--');

            
%%% -3dB torque
            contW = w0 * sqrt( ( contSpeeds * k * I ) ./ ( 0.5* torque + contSpeeds * k * I ) );
            contF = contW / 2 / pi;

            peakW = w0 * sqrt( ( peakSpeeds * k * I ) ./ ( 0.5* torque + peakSpeeds * k * I ) );
            peakF = peakW / 2 / pi;

            plot(torque,contF,'b');
            plot(torque,peakF,'r');
            
%%% output toruqe instead of input torque, yields the frequency behavior of
%%% the sensor!
%             contW = k*contSpeeds./torque;
%             contF = contW / 2 / pi;
% 
%             peakW = k*peakSpeeds./torque;
%             peakF = peakW / 2 / pi;
% 
%             plot(torque,contF,'bx');
%             plot(torque,peakF,'rx');
%             
%%%
            plot(torque, f0, 'k--' )
            
%             ymax = k*dq_NL /t_NL;
%             
            xlim([0,xmax]);
%             ylim([0,ymax]);
            xlabel('torque [Nm]')
            ylabel('frequency [Hz]')
            
        end
        
        function h = draw_torque_speed_curve(this)
            % draw_torque_speed_curve Displays speed-torque-curve in a
            % figure.
            %
            %   fHandle = dsg.draw_torque_speed_curve
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
            d_cl = this.jointModel.d_cl;
            d_m  = this.jointModel.d_m;
            d_g  = this.jointModel.d_g;
            d_l  = this.jointModel.d_l;
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
            xmax = 1.2 * t_p;
            ymax = 1.02 * dq_p;
            h = figure;
            hold on
                       
            % torque speed line
            mVals = (0:1/this.nPlotVals:1) * t_stall;
            linCurve = dq_0 - slope * mVals;
            plot(mVals, linCurve, 'k', 'DisplayName', 'Torque-Speed Line')
            
            % Rated operating point
            plot(t_r, dq_r, 'ko', 'DisplayName', 'Nominal Operating Point')
            
            % No-Load operating point
            plot(t_NL, dq_NL, 'ko', 'DisplayName', 'No-Load Operating Point')
            
           
            % Friction
            speedVals = (0:1/this.nPlotVals:1) * dq_p;
            Mc = d_cm + d_cg + d_cl;            % Static friction
            Mv = (d_m + d_g + d_l) * speedVals; % Velocity dependent friction
            Mf = Mc + Mv;

            % Continuous operating range
            tOp =  [0, t_r,  t_r,    0].';
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
            
            % Continuous operation 
            speedVals = p_cm ./ mVals;
            plot(mVals,speedVals, 'b-', 'DisplayName', 'Rated Mechanical Power')
%             speedVals = this.p_peakm ./ mVals;
%             plot(mVals,speedVals, 'r--', 'DisplayName', 'Peak Mechanical Power')
            plot([0,xmax], dq_r * [1,1], 'b--', 'DisplayName', 'Maximum Continous Speed')
            plot(t_r*[1,1], [0,ymax], 'b--', 'DisplayName', 'Maximum Continous Torque')
            
            
            maxSpeeds = this.computeMaxPeakSpeed(mVals);
            contSpeeds = this.computeMaxContSpeed(mVals);
            
            plot(mVals, maxSpeeds,'g')
            plot(mVals, contSpeeds,'m')
            
            % Annotations and Figure Style
            xlim([0,xmax]);
            ylim([0,ymax]);
            xlabel('torque [Nm]')
            ylabel('speed [rad/s]')
            box on
            
            %legend show;
        end
                
        
        function contSpeed = computeMaxContSpeed(this, torque)
            
            % shorthands
            t_r = this.jointModel.t_r;
            p_cm = this.jointModel.p_cm;
            
            % start off from the peak speed
            contSpeed = this.computeMaxPeakSpeed(torque);

            % rated torque
            peakSpeed(torque > t_r) = 0;
            
            % continuous power
            powerSpeed = p_cm ./ torque;
            
            contSpeed = min(contSpeed,powerSpeed);
            
            
        end
        
        function peakSpeed = computeMaxPeakSpeed(this, torque)

            
            % shorthands
            slope = this.jointModel.dq_over_dm;
            dq_0 = this.jointModel.dq_0;
            
            dq_p = this.jointModel.dq_p;
            d_cm = this.jointModel.d_cm;
            d_cg = this.jointModel.d_cg;
            d_cl = this.jointModel.d_cl;
            d_m  = this.jointModel.d_m;
            d_g  = this.jointModel.d_g;
            d_l  = this.jointModel.d_l;
            
            p_pm = this.jointModel.p_pm;
            t_p = this.jointModel.t_p;


            
            % first, look along the torque speed curve defined by the
            % electrical subsystem
            elSpeed = dq_0 - slope * torque;

            
            % Friction
            kC = d_cm + d_cg + d_cl;            % Static friction
            kV = (d_m + d_g + d_l); % Velocity dependent friction

            fricSpeed = (torque - kC)/kV;
            fricSpeed(torque < kC) = 0;

            peakSpeed = min(elSpeed, fricSpeed);
            
            % Power curve
            powerSpeed = p_pm ./ torque;
            
            peakSpeed = min(peakSpeed, powerSpeed);

            
        end

        
        function createDataSheet(this)
            % CREATEDATASHEET Main function to trigger the data sheet
            % generation.
            %
            %   dsg.createDataSheet
            %
            % As a result of this function call, a data sheet file will be
            % created and stored in dsg.outputDir. The filename will have 
            % the form  "jointModel.name_dsg.outFName".
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
            
            % Get local copies of all template files.
            copyfile([this.templateDir, filesep, this.texFName],this.texFName);
            copyfile([this.templateDir, filesep, this.clsFName],this.clsFName);
            
            % Collect all information about the actuator and prepare macros
            % for the use in the Tex source.
            this.createDefFile;

            % Produce all plots
            this.makeDataSheetPlots;
            
            % Compile the Tex source
            this.compileTexFile
            
            % Put the generated datasheet in the final destination.
            [~, fName] = fileparts(this.texFName);
            if ~exist(this.outputDir,'dir')
                mkdir(this.outputDir)
            end
            copyfile([fName,'.pdf'],[this.outputDir,filesep, this.assembleOutFileName ]);

            % Clean up.
            delete([fName,'.*'])
            delete([this.torqueSpeedFName])
            
        end
        
        function fName = assembleOutFileName(this)
            % ASSEMBLEOUTFILENAME Function to assembel the datasheet file
            % name.
            %
            %   fName dsg.assembleOutFileName
            %
            %
            % Inputs:
            %
            % Outputs:
            %   fName: Name for the datasheet file.
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
            
            fName = [this.jointModel.name, '_', this.outFName];
            
        end
        
        function makeDataSheetPlots(this)
            % MAKEDATASHEETPLOTS Calls the individual routines that plot
            % the actuator characteristics.
            %
            %   dsg.makeDataSheetPlots
            %
            % The function will create a number of figure files in the
            % output directory specified by dsg.outputDir.
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
            
            
            % Torque-Speed curve
            h = this.draw_torque_speed_curve;
                        
            set(gcf,'Units','centimeters');
            set(gcf,'PaperUnits','centimeters');
            pos = get(gcf,'Position');
            pos(3) = 18;
            pos(4) = 8;
             
            set(gcf,'Position',pos)
            set(h,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])          
            set(gca,'LooseInset',get(gca,'TightInset'))
            
            printpdf(gcf,this.torqueSpeedFName,['-r',num2str(this.plotResolution)])

            close(h);
             
        end
        
        function createDefFile(this)
            % CREATEDEFFILE Creates the temporary file containing the macro 
            % definitions for the datasheet sources.
            % 
            %   dsg.createDefFile
            %
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
            %
            
            jM = this.jointModel; % Shorthand
            
            fid = fopen(this.cfgFName,'w+');
            %
            fprintf(fid,'\\def \\valJointName{%s}\n'                 , jM.name                          );
            %
            fprintf(fid,'%s\n','% Mechanical Properties');
            fprintf(fid,'%s\n','%');
            %
            %Dimensions
            %
            fprintf(fid,'\\def \\valDiameter{%4.1f}\n'               , jM.diam                          );
            fprintf(fid,'\\def \\symDiameter{%s}\n'                  , 'D'                              );
            %
            fprintf(fid,'\\def \\valActlength{%4.1f}\n'              , jM.len                           );
            fprintf(fid,'\\def \\symActlength{%s}\n'                 , 'l'                              );
            %
            % Inertae
            %
            fprintf(fid,'\\def \\valMass{%6.4f}\n'                   , jM.m                             );
            fprintf(fid,'\\def \\symMass{%s}\n'                      , 'm'                              );
            %
            fprintf(fid,'\\def \\valInertiarotor{%6.4f}\n'           , jM.I_m                           );
            fprintf(fid,'\\def \\symInertiarotor{%s}\n'              , 'I_m'                            );
            %
            fprintf(fid,'\\def \\valInertiagear{%6.4f}\n'            , jM.I_g                           );
            fprintf(fid,'\\def \\symInertiagear{%s}\n'               , 'I_g'                            );
            %
            fprintf(fid,'\\def \\valInertiaspring{%6.4f}\n'          , jM.I_l                           );
            fprintf(fid,'\\def \\symInertiaspring{%s}\n'             , 'I_m'                            );
            % 
            fprintf(fid,'\\def \\valTmech{%6.4f}\n'                  , jM.T_mech                        );
            fprintf(fid,'\\def \\symTmech{%s}\n'                     , 'T_{mech}'                       );
            %
            % Stiffnesses
            %
            fprintf(fid,'\\def \\valSpringstiffness{%5d}\n'          , jM.k_b                           );
            fprintf(fid,'\\def \\symSpringstiffness{%s}\n'           , 'k_b'                            );
            fprintf(fid,'\\def \\valGearstiffness{%5d}\n'            , jM.k_g                           );
            fprintf(fid,'\\def \\symGearstiffness{%s}\n'             , 'k_g'                            );
            %
            % Friction
            %
            fprintf(fid,'\\def \\valViscousdamping{%6.4f}\n'         , jM.d_m+jM.d_g+jM.d_l             );
            fprintf(fid,'\\def \\symViscousdamping{%s}\n'            , 'd_v'                            );
            %
            fprintf(fid,'\\def \\valCoulombdamping{%6.4f}\n'         , jM.d_cm+jM.d_cg+jM.d_cl          );
            fprintf(fid,'\\def \\symCoulombdamping{%s}\n'            , 'd_c'                            );
            %
            fprintf(fid,'\\def \\valStribeckdamping{%6.4f}\n'        , jM.d_s                           );
            fprintf(fid,'\\def \\symStribeckdamping{%s}\n'           , 'd_s'                            );
            %
            fprintf(fid,'\\def \\valStribeckspeed{%6.4f}\n'          , jM.v_s                           );
            fprintf(fid,'\\def \\symStribeckspeed{%s}\n'             , '\dot{\theta}_s'                 );
            %
            % Electrical
            %
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Electrical Properties');
            fprintf(fid,'%s\n','%');
            %
            fprintf(fid,'\\def \\valTorqueconstant{%6.4f}\n'          , jM.k_t                          );
            fprintf(fid,'\\def \\symTorqueconstant{%s}\n'             , 'k_\tau'                        );
            %
            fprintf(fid,'\\def \\valGeneratorconstant{%6.4f}\n'       , jM.k_w                          );
            fprintf(fid,'\\def \\symGeneratorconstant{%s}\n'          , 'k_\omega'                      );
            %
            fprintf(fid,'\\def \\valArmatureresistance{%6.4f}\n'      , jM.r                            );
            fprintf(fid,'\\def \\symArmatureresistance{%s}\n'         , 'r_A'                           );
            %
            fprintf(fid,'\\def \\valArmatureinductance{%6.4f}\n'      , jM.x                            );
            fprintf(fid,'\\def \\symArmatureinductance{%s}\n'         , 'x_A'                           );
            %
            fprintf(fid,'\\def \\valTel{%6.4f}\n'                     , jM.T_el                         );
            fprintf(fid,'\\def \\symTel{%s}\n'                        , 'T_{el}'                        );
            %
            % Thermal
            %
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Thermal Properties');
            fprintf(fid,'%s\n','%');
            %
            fprintf(fid,'\\def \\valResthermWH{%4.2f}\n'              , jM.r_th1                        );
            fprintf(fid,'\\def \\symResthermWH{%s}\n'                 , 'r_{th1}'                       );
            %
            fprintf(fid,'\\def \\valResthermHA{%4.2f}\n'              , jM.r_th2                        );
            fprintf(fid,'\\def \\symResthermHA{%s}\n'                 , 'r_{th2}'                       );
            %
            fprintf(fid,'\\def \\valTthw{%4.2f}\n'                    , jM.T_thw                        );
            fprintf(fid,'\\def \\symTthw{%s}\n'                       , 'T_{th,\:w}'                    );
            %
            fprintf(fid,'\\def \\valTthm{%4.2f}\n'                    , jM.T_thm                        );
            fprintf(fid,'\\def \\symTthm{%s}\n'                       , 'T_{th,\:m}'                    );
            %
            fprintf(fid,'\\def \\valTmpWindMax{%4.2f}\n'              , jM.Tmp_WMax                     );
            fprintf(fid,'\\def \\symTmpWindMax{%s}\n'                 , '\nu_{w,\:Max}'                 );
            %
            fprintf(fid,'\\def \\valTmpANom{%4.2f}\n'                 , jM.Tmp_ANom                     );
            fprintf(fid,'\\def \\symTmpANom{%s}\n'                    , '\nu_{A,\:Nom}'                 );
            %
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Rated Operation');
            fprintf(fid,'%s\n','%');
            %
            fprintf(fid,'\\def \\valRatedvoltage{%4.2f}\n'            , jM.v_0                                          );
            fprintf(fid,'\\def \\symRatedvoltage{%s}\n'               , 'V_0'                                           );
            %
            fprintf(fid,'\\def \\valRatedcurrent{%4.2f}\n'            , jM.i_c                                          );
            fprintf(fid,'\\def \\symRatedcurrent{%s}\n'               , 'i_R'                                           );
            %
            fprintf(fid,'\\def \\valRatedtorque{%4.2f}\n'             , jM.t_r                                          );
            fprintf(fid,'\\def \\symRatedtorque{%s}\n'                , '\tau_R'                                        );
            %
            fprintf(fid,'\\def \\valRatedspeed{%4.2f}\n'              , jM.dq_r                                         );
            fprintf(fid,'\\def \\symRatedspeed{%s}\n'                 , '\dot{\theta}_R'                                );
            %
            fprintf(fid,'\\def \\valRatedpowere{%4.2f}\n'             , jM.p_ce                                         );
            fprintf(fid,'\\def \\symRatedpowere{%s}\n'                , 'P_{R,\:E}'                                     );
            %
            fprintf(fid,'\\def \\valRatedpowerm{%4.2f}\n'             , jM.p_cm                                         );
            fprintf(fid,'\\def \\symRatedpowerm{%s}\n'                , 'P_{R,\:M}'                                     );
            %
            fprintf(fid,'\\def \\valNoloadcurrent{%4.2f}\n'           , jM.i_NL                                         );
            fprintf(fid,'\\def \\symNoloadcurrent{%s}\n'              , 'i_{NL}'                                        );
            %
            fprintf(fid,'\\def \\valNoloadtorque{%4.2f}\n'            , jM.t_NL                                         );
            fprintf(fid,'\\def \\symNoloadtorque{%s}\n'               , '\tau_{NL}'                                     );
            %
            fprintf(fid,'\\def \\valNoloadspeed{%4.2f}\n'             , jM.dq_NL                                        );
            fprintf(fid,'\\def \\symNoloadspeed{%s}\n'                , '\dot{\theta}_{NL}'                             );
            %
            fprintf(fid,'\\def \\valStalltorque{%4.2f}\n'             , jM.t_stall                                      );
            fprintf(fid,'\\def \\symStalltorque{%s}\n'                , '\tau_{Stall}'                                  );
            %
            fprintf(fid,'\\def \\valStartingcurrent{%4.2f}\n'         , jM.i_start                                      );
            fprintf(fid,'\\def \\symStartingcurrent{%s}\n'            , 'i_{start}'                                     );
            %
            fprintf(fid,'\\def \\valSpeedtorquegradient{%6.4f}\n'     , jM.dq_over_dm                                   );            
            fprintf(fid,'\\def \\symSpeedtorquegradient{%s}\n'        , '\text{d}\dot{\theta}/\text{d}\tau'             );
            %
            fprintf(fid,'%s\n','%');
            fprintf(fid,'%s\n','% Peak Operation');
            fprintf(fid,'%s\n','%');
            %
            fprintf(fid,'\\def \\valMaxcurrent{%4.2f}\n'              , jM.i_p                    );
            fprintf(fid,'\\def \\symMaxcurrent{%s}\n'                 , 'i_p'                     );
            %
            fprintf(fid,'\\def \\valMaxtorque{%4.2f}\n'               , jM.t_p                    );
            fprintf(fid,'\\def \\symMaxtorque{%s}\n'                  , '\tau_p'                  );
            %
            fprintf(fid,'\\def \\valMaxspeed{%4.2f}\n'                , jM.dq_p                   );
            fprintf(fid,'\\def \\symMaxspeed{%s}\n'                   , '\dot{\theta}_p'          );
            %
            fprintf(fid,'\\def \\valMaxpowere{%4.2f}\n'               , jM.p_pe                   );
            fprintf(fid,'\\def \\symMaxpowere{%s}\n'                  , 'P_{p,\:E}'               );
            %
            fprintf(fid,'\\def \\valMaxpowerm{%4.2f}\n'               , jM.p_pm                   );
            fprintf(fid,'\\def \\symMaxpowerm{%s}\n'                  , 'P_{p,\:M}'               );
            %
            fprintf(fid,'\\def \\valGearratio{%d:1}\n'                , jM.n                      );
            fprintf(fid,'\\def \\symGearratio{%s}\n'                  , 'N'                       );
            %
            %            fprintf(fid,'\\def \\maxefficiency{%4.2f}\n',jM.v_0);
            %
            fclose(fid);
        end
        
        function [flag, cmdout] = compileTexFile(this)
            % COMPILETEXFILE Compiles the datasheet Tex source files using
            % LuaLatex.
            % 
            %   [status, cmdout] = dsg.compileTexFile
            %
            % The function calls LuaLatex to compile the Tex sources
            % forming the contents of the datasheet to be generated. This
            % requires that you have LuaLatex installed and registered in
            % your search path.
            % 
            % Make sure no PDF viewer blocks the overwriting of any existing 
            % previously generated datasheet file.
            %
            % Outputs:: 
            %  status: Status flag that is zero upon success and nonzero
            %  integer otherwise. 
            %
            %  cmdout: Output of the operating system command returned as
            %  as string.
            %
            % The outputs are the outputs passed by the Matlab built-in 
            % system command. See 'doc system' for details.
            %
            %
            % Notes::
            %   This method has been tested with MiKTeX version 2.9.
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also createDataSheet, genericJoint, jointBuilder.    
            %
            
            % Tested with MiKTeX/2.9
            cmd = ['lualatex.exe -synctex=-1 -interaction=nonstopmode ', this.texFName];
            
            % Make the system call invoking the LuaLatex compiler.
            [flag, cmdout] = system(cmd);
                        
        end
        
        
    end
    
end

