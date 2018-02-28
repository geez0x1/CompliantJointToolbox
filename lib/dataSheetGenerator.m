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
        nPlotVals = 10000;                             % Number of samples used to produce the graphs in the datasheet figures.
        plotNormalized = 0;                            % Plot graphs in normalized quantities (default: false)
        freqMax = 100;                                 % Maximum frequency for spectral plots in Hz (default: 100 Hz)
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
            this.templateDir = [dsgRoot(1:end-4), filesep, 'templates' ];
            
            this.outputDir = ['.', filesep];
            
        end
        
        function speed = torqueSpeedCurve(this,tau)
            % TORQUESPEEDCURVE Computes the maximum feasible motor speed
            % corresponding for a given torque level.
            %
            %   speed = dsg.TorqueSpeedCurve(tau)
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
        
        
        function [h, hAx, hLine1, hLine2] = drawEfficiencyCurve(this)
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
            
            % Linearly spaced vector of torques between 0 and peak torque
            tauVals = (0:1/this.nPlotVals:1) * t_p;
            
            % Power computation
            PL = v_0 / k_t/N * tauVals - r / k_t^2/N^2 * tauVals.^2;
            P_tot = v_0 / k_t/N * tauVals;
            speedVals = this.torqueSpeedCurve(tauVals);
            
            % Efficiency computation
            eta = 100*( speedVals .* (tauVals - t_NL) * N * k_t / v_0 ) ./ tauVals;
            eta_max = 100 * ( 1 - sqrt( t_NL / t_stall)  );
            
            % Initialize figure and plot
            h = gcf;
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

        function [h, hAx, hLine1, hLine2] = drawThermalCharacteristics(this)
            % DRAW_THERMAL_CHARACTERISTICS Creates a plot with two y-axes.
            % One axis displays the final temperature at the given
            % operating condition. The second axis displays the time to
            % reach the maximum allowed temperature.
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
            
            % Linearly spaced vector of torques between 0 and peak torque
            tauVals = (0:1/this.nPlotVals:1) * t_p; 

            % Steady State Temperature
            resCoeff = (r_th1 + r_th2) * r_TA / N^2 / k_t^2;
            Tmp_W = Tmp_ANom + ...
                (  ( resCoeff * tauVals.^2 ) ./ ( 1 - aCU * resCoeff *  min(tauVals,t_r).^2 )  );
            
            % Saturate temperature at melting point
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
            
            % Initialize figure and plot
            h = gcf;
            [hAx,hLine1,hLine2] = plotyy(tauVals(:), [Tmp_W(:), TmpLimW(:),TmpLimCU(:) ], ... Steady State Temperature
                tauVals(:),[tCrit]);                % Time to critical temperature
             
            % Manipulate line style
            set(hLine1,'Color','r','linewidth',1.5)
            set(hLine1(2),'LineStyle','--','linewidth',1)
            set(hLine1(3),'LineStyle',':','linewidth',1)
            ylabel(hAx(1), 'Steady State Temp. [^\circC]','Color','r');
            set(hAx(1),'ylim',[0,1.05*max((Tmp_W))]);

             % Manipulate plot style
             set(hLine2,'Color','b','linewidth',1.5)
             ylabel(hAx(2),'Time to Crit Temperature [s]','Color','b');
             set(hAx(2),'ylim',[0,1.05*max(real(tCrit))]);

             set(hAx,'xlim',[0,t_p]);  
             
             legend({'Steady state temperature [C]',...
                     'Critical winding temperature',...
                     'Copper melting temperature',...
                     'Time to critical winding temperature [s]'}, 'location','best')
                 
             xlabel('Torque [Nm]')
             title(['Steady-state temperature rise for ambient temperature ',num2str(Tmp_ANom),' [C]'])
             
             hold on;
             % Plot rated operation area
             Xr = [0 t_r t_r 0].';
             Yr = [0 0   Tmp_WMax Tmp_WMax].';
             rated_area = fill(Xr, Yr, 0.9* [0.2 0.2 1],'LineStyle','none');
             
             
             alpha(0.25)     % add some transparency
             plot(hAx(1),[t_r,t_r],[0,TmeltCU],'color','b','linestyle','--')
             hold off;
             
             % Plot peak operation area
             axes(hAx(2))
             hold on
             loc_subsamp = 10;
             Xp = [t_r tauVals(idx+1:loc_subsamp:end) fliplr(tauVals(idx+1:loc_subsamp:end)) t_r].';
             Yp = [0 tCrit(idx+1:loc_subsamp:end) zeros(size(fliplr(tauVals(idx+1:loc_subsamp:end)))) 0].';
             peak_area = fill(Xp, Yp, 0.8* [1 0 0],'LineStyle','none');
             alpha(0.25)     % add some transparency
             hold off
             
        end
                
        function h = drawTorqueFrequencyCurve(this, subtractFriction)
            
            % Check input parameters
            if ~exist('subtractFriction','var') % Should the friction torque be subtracted?
                subtractFriction = 1;
            end
                        
            % SHORTHANDS
            d_cm = this.jointModel.d_cm; % Coulomb frictions
            d_cg = this.jointModel.d_cg;
            d_cl = this.jointModel.d_cl;
            d_m  = this.jointModel.d_m;  % Viscous dampings
            d_g  = this.jointModel.d_g;
            d_l  = this.jointModel.d_l;
            d_gl = this.jointModel.d_gl;
            k = this.jointModel.k_b;     % Spring stiffness
            t_p = this.jointModel.t_p;   % Peak torque
            t_r = this.jointModel.t_r;   % Rated torque
            k_t = this.jointModel.k_t;   % Torque constant
            N = this.jointModel.n;       % gear ratio
            slope = this.jointModel.dq_over_dm; % torque speed slope
            dq_0 = this.jointModel.dq_0; % no torque speed
            v_0 = this.jointModel.v_0;   % supply voltage
            r = this.jointModel.r;       % electrical winding resistance
           
            % Compute characteristic parameters
            I = this.jointModel.I_m + this.jointModel.I_g; % inertia
            w0 = sqrt(k/I);                     % Resonance Frequency in rad/s
            f0 = w0/2/pi;                       % Resonance Frequency in Hz
            Mc = d_cm + d_cg + d_cl;            % Static friction
            dv = (d_m + d_g + d_l);             % Velocity dependent friction
            magDrop = db2mag(-3);               % -3dB, allowed magnitude drop for tracking
            
            % NORMALIZATION (by peak torque and natural frequency)
            if this.plotNormalized
                wNorm = 1/w0;
                fNorm = 2*pi*wNorm;
                tNorm = 1/t_p;
                xlabelStr = '$|\tau / \tau_p|$ [.]';
                ylabelStr = ' $f_c / f_0$ [.]';
            else
                fNorm = 1;
                tNorm = 1;
                xlabelStr = '$\tau$ [Nm]';
                ylabelStr = '$f_c$ [Hz]';
            end
            
            % Linearly spaced vectors of torques and normalized frequencies
            ymax =  this.freqMax;
            xmax = 1.1 * t_p;   
            torque = (1/this.nPlotVals:1/this.nPlotVals:1) * xmax;
            wn = (1/this.nPlotVals:1/this.nPlotVals:1) * 2*pi*ymax / w0;

            % LIMITATION DUE TO MAGNITUDE GAIN
            allTF = this.jointModel.getTF;
            % Torque transfer function, convert current input to torque
            torqueTF = allTF(7,1)/ k_t / N;  
            
            % Compute magnitude at given frequencies
            magGain = bode(torqueTF,wn*w0); 
            magGain = magGain(:);
            
            % Subtract friction toruqe if desired
            if subtractFriction
                dq_t_r = wn.*w0.*magGain.*t_r./k;
                dq_t_p = wn.*w0.*magGain.*t_p./k;
                fCorrC = Mc + dv.*dq_t_r;
                fCorrP = Mc + dv.*dq_t_p;
            else
                fCorrC = 0;
                fCorrP = 0;
            end
           
            % LIMITATION DUE TO TORQUE SPEED CURVE
            % This is the bandwidth limit due to the motor back-emf
            % generation and limited supply voltage
            contSpeeds = this.computeMaxContSpeed(torque);
            peakSpeeds = this.computeMaxPeakSpeed(torque);
            
            % Subtract friction toruqe if desired
            if subtractFriction
                fCorrC = Mc + dv*contSpeeds;
                fCorrP = Mc + dv*peakSpeeds;
            else
                fCorrC = 0;
                fCorrP = 0;
            end
            
            % Torque sensor transfer function
            springTF = tf([d_gl, k], [1 0]);
            springMag = bode(springTF,wn*w0);
            springMag = springMag(:);
            
            %% TOTAL FREQUENCY LIMIT
            fSpring = interp1(dq_0*springMag , wn*w0/2/pi, torque);
            fAmp    = interp1(magGain*t_p    , wn*w0/2/pi, torque);
             
            fTotal = min(fAmp,fSpring);
             
             
            % THERMAL TIME LIMITATION
            r_th2 = this.jointModel.r_th2;
            r     = this.jointModel.r;
            a_cu  =  0.0039;
            k_t   = this.jointModel.k_t;
            n     = this.jointModel.n;
            T_max = this.jointModel.Tmp_WMax;
            t_w   = this.jointModel.T_thw;
            t_r   = this.jointModel.t_r;
            
            % Torque linear space (some more lines for lower torques)
            tau = [0:9, 10:20:(t_p)];
            nLines = numel(tau);
            
            % Compute time to maximal temperature
            i_m2  = 0.25 * (tau / k_t / n).^2; % effective value -> multiply with 0.5            
            deltaTw =  ( r_th2 * r * i_m2 ) ./ (1 - a_cu * r_th2 * r *i_m2  );
            t_max = t_w * log( T_max ./ deltaTw );
                        
             %% PLOTTING
            figHandle = gcf;
            hold on

            % Contour plot
            %   Prepare data
            Z = repmat(magGain,1,nLines) .* tau;
            TMAX = repmat(magGain,1,nLines).*t_max;
            [F,TAU] = meshgrid(wn*w0/2/pi*fNorm, tau*tNorm/magDrop);
            %   Create plot
            [~,h] = contour(TAU.',F.', (TMAX),t_max(tau>t_r));
            %   Configure plot
            h.Fill = 'on';  % instead of lines use shaded colour surfaces
            set(h,'LineStyle','none')
            set(gcf,'renderer','painters') 
            alpha(0.25)     % add some transparency to make lines plotted on top clearly visible
            %   adjust colormap from cold (blue) to hot (red)
            nCVals = 64;    
            myMap = [(nCVals:-1:0).', 0*ones(nCVals+1,1),(0:nCVals).' ]/nCVals;
            colormap(myMap)
            %   add a color bar legend
            c = colorbar;      
            cpos = c.Position;
            % A hack to give the color bar the proper appearance. The color
            % bar does not feature transparency anymore.
            annotation('textbox',...
                c.Position,...
                'FitBoxToText','off',...
                'FaceAlpha',0.25,...
                'EdgeColor',[1 1 1],...
                'BackgroundColor',[1 1 1]);
            % Saturate the time to for rated operation
            tmp = c.TickLabels;
            tmp{end} = '\infty';
            c.TickLabels = tmp;
            c.Label.Interpreter = 'latex';
            c.Label.String = '$t_{max}$ [s]';

             
            % Plot total frequency limit (back-emf/voltage and sensor)
            plot(torque(torque<=t_p)*tNorm/magDrop,fTotal(torque<=t_p)*fNorm,'Color', [0.8 1 0.8],'LineWidth',2)
            
            % Plot torque transfer function magnitude
            plot((magGain*t_p - fCorrP)*tNorm/magDrop,    wn*w0/2/pi * fNorm,'k','LineWidth',0.5) % -3 dB
            
            % Plot voltage saturation limit
            plot(dq_0*springMag *tNorm/ magDrop , wn*w0/2/pi *fNorm,'k--','LineWidth',0.5)
            
            % plot resonance frequency
            plot(torque*tNorm, f0*ones(size(torque))*fNorm, '--','color',0.8*[1 1 1],'LineWidth',0.5 )

            %% Plot appearance and labels
            xlim([0,xmax]*tNorm);
            ylim([0,ymax]*fNorm);
            xlabel(xlabelStr,'Interpreter','latex')
            ylabel(ylabelStr,'Interpreter','latex')
            box on;
            set(gca,'TickLabelInterpreter','latex');
        end
        
        function h = drawTorqueSpeedCurve(this, subtractFriction, legendFontSize)
            % draw_TorqueSpeedCurve Displays speed-torque-curve in a
            % figure.
            %
            %   fHandle = dsg.draw_TorqueSpeedCurve(subtractFriction, legendFontSize)
            %
            % Inputs:
            %   subtractFriction: flag that controls if the torque axis
            %                     is the torque after friction subtraction
            %                     (true, default) or the torque including
            %                     fricion (false)
            %   legendFontSize:  Optional specifier for the friction font
            %                    size in pt. If unspecified, default font
            %                    size is used.
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
            
            if ~exist('legendFontSize','var')
                legendFontSize = get(0, 'DefaultTextFontSize');
            end
            
            if ~exist('subtractFriction','var')
                subtractFriction = 0;
            end
            
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
            
            h = gcf;
            hold on
                    
            % Compute Friction
            speedVals = (0:1/this.nPlotVals:1) * dq_p;
            Mc = d_cm + d_cg + d_cl;  % Static friction
            dv = (d_m + d_g + d_l);   % Velocity dependent friction
            Mv = dv * speedVals; 
            Mf = Mc + Mv;             % Total friction
            
            % Torque-Speed Gradient
            mVals = (0:1/this.nPlotVals:1) * t_stall;
            linCurve = dq_0 - slope * mVals;
            
            % Speed at peak torque
            dq_tp = dq_0 - slope * t_p;
            
            if subtractFriction
                fCorr =  Mc + dv * linCurve;
                fCorr_r = Mc + dv * dq_r;
                fCorr_NL = Mc + dv * dq_NL;
                fCorr_tp = Mc + dv * dq_tp;
            else
                fCorr    = 0;
                fCorr_r  = 0;
                fCorr_NL = 0;
                fCorr_tp = 0;
            end
            
            if this.plotNormalized
                tNorm = 1/t_p;
                dqNorm = 1/dq_NL;
                xlabelStr = '$\tau / \tau_p$ [.]';
                ylabelStr = '$\dot{q} / \dot{q}_{NL}$ [.]';
            else
                tNorm = 1;
                dqNorm = 1;
                xlabelStr = '$\tau$ [Nm]';
                ylabelStr = '$\dot{q}$ [rad/s]';
            end
            
            plot( (mVals - fCorr )*tNorm, linCurve*dqNorm, 'k', 'DisplayName', 'Torque-Speed Boundary','linewidth',1.5)
                        
            % OPERATING POINTS
            % Rated operating point
            plot((t_r - fCorr_r)*tNorm, dq_r*dqNorm, 'ko', 'DisplayName', 'Nominal Operating Point')
            
            % No-Load operating point
            plot((t_NL - fCorr_NL)*tNorm, dq_NL*dqNorm, 'ko', 'DisplayName', 'No-Load Operating Point')
            
            % Peak operating point
            plot((t_p - fCorr_tp)*tNorm, dq_tp*dqNorm, 'ko', 'DisplayName', 'Peak Operating Point')
            
            % COLOR FRICTION AND CONTINUOUS OPERATION REGIONS
            % Continuous operating range
            if subtractFriction
                fCorr_1 = Mc + dv * dq_r;
                fCorr_2 = Mc + dv * dq_0;
            else
                fCorr_1 = 0;
                fCorr_2 = 0;
            end
            tOp =  [0, t_r,  t_r-fCorr_1,    0-fCorr_2].';
            dqOp = [0,   0, dq_r, dq_0].';
            fill(tOp*tNorm, dqOp*dqNorm, 0.9* [0.2 0.2 1],'LineStyle','none')
            
            % Friction
            if ~subtractFriction
                fill([Mf, 0, 0]*tNorm, ...
                    [speedVals speedVals(end) 0]*dqNorm,...
                    0.8* [1 0 0],...
                    'LineStyle','none')
                plot(Mf*tNorm, speedVals*dqNorm, 'k:', 'DisplayName', 'Friction Torque','linewidth',1.5)
            end
            alpha(0.25)     % add some transparency
            
            % OPERATING BOUNDARIES
            % Peak Operation
            speedVals = p_pm ./ mVals;
            if subtractFriction
                fCorr_1 = Mc + dv * ymax;
                fCorr_2 = Mc + dv * speedVals;
            else
                fCorr_1 = 0;
                fCorr_2 = 0;
            end
            plot([t_p, t_p-fCorr_1]*tNorm,[0, ymax]*dqNorm,'r--', 'DisplayName', 'Peak Torque','linewidth',1.5)
            plot((mVals-fCorr_2)*tNorm,speedVals*dqNorm, 'r-', 'DisplayName', 'Peak Mechanical Power','linewidth',1)
            
            % Continuous operation 
            speedVals = p_cm ./ mVals;
            if subtractFriction
                fCorr_1 = Mc + dv * ymax;
                fCorr_2 = Mc + dv * speedVals;
            else
                fCorr_1 = 0;
                fCorr_2 = 0;
            end
            plot((t_r*[1,1] - fCorr_1*[0, 1])*tNorm, [0,ymax]*dqNorm, 'b--', 'DisplayName', 'Maximum Continous Torque','linewidth',1.5)
            plot((mVals-fCorr_2)*tNorm,speedVals*dqNorm, 'b-', 'DisplayName', 'Rated Mechanical Power','linewidth',1)
            
            
            
%             % ANNOTATIONS AND FIGURE STYLE
            xlim([0,xmax]*tNorm);
            ylim([0,ymax]*dqNorm);
            xlabel(xlabelStr,'interpreter','latex')
            ylabel(ylabelStr,'interpreter','latex')
            box on
            
            % CREATE COSTUMIZED LEGEND
            if legendFontSize ~= 0
                % The standard legend functionality provided by Matlab is
                % inconvenient for multiple reasons. It covers parts of the 
                % plot and positioning/resizing is painful...
                % So we have to create some custom solution. The custom 
                % solution uses an invisible dummy axes object.
                %
                % First, we place the main axes in normalized coordinates:
                pos = [0.11, 0.13 0.70, 0.85];
                set(gca,'position',pos)

                % Next, we place the dummy axes object:
                axSpace = 0.02;  % the space between the two axes
                dummyAx = axes;  
                newPos = [pos(1)+pos(3)+axSpace,pos(2),0.18-axSpace,pos(4)];
                set(dummyAx, 'position',newPos,'visible','off')

                % The legend entries are implemented with annotation objects.
                % Each entry is a pair of a line annotation and a textbox
                % annotation.
                fsize = 7;
                    
                % Torque Speed Gradient
                annotation(gcf,'line',[0.83 0.88],... % Create line
                    [0.93 0.93]);
                annotation(gcf,'textbox',...          % Create textbox
                    [0.88 0.91 0.07 0.06],...
                    'String',{'Torque Speed Boundary'},...
                    'FontSize', legendFontSize, ...
                    'FitBoxToText','off',...
                    'LineStyle','none');

                % Rated speed and torque limits
                posDec = 0.19;
                annotation(gcf,'line',[0.82 0.88],... % Create line
                    [0.93 0.93]-posDec,...
                    'LineStyle','--',...
                    'Color',[0 0 1]);
                annotation(gcf,'textbox',...          % Create textbox
                    [0.88 0.91-posDec 0.07 0.06],...
                    'String',{'Rated Torque'},...
                    'FontSize', legendFontSize, ...
                    'FitBoxToText','off',...
                    'LineStyle','none');

                % Peak speed and torque limits
                posDec = posDec + 0.14;
                annotation(gcf,'line',[0.82 0.88],...   % Create line
                    [0.93 0.93]-posDec,...
                    'LineStyle','--',...
                    'Color',[1 0 0]);
                annotation(gcf,'textbox',...            % Create textbox
                    [0.88 0.91-posDec 0.07 0.06],...
                    'String',{'Peak Torque'},...
                    'FontSize', legendFontSize, ...
                    'FitBoxToText','off',...
                    'LineStyle','none');

                % Rated power limits
                posDec = posDec + 0.14;
                annotation(gcf,'line',[0.83 0.88],...  % Create line
                    [0.93 0.93]-posDec,...
                    'LineStyle','-',...
                    'Color',[0 0 1]);
                annotation(gcf,'textbox',...           % Create textbox
                    [0.88 0.91-posDec 0.07 0.06],...
                    'String',{'Rated Power'},...
                    'FontSize', legendFontSize, ...
                    'FitBoxToText','off',...
                    'LineStyle','none');

                % Peak power limits
                posDec = posDec + 0.14;
                annotation(gcf,'line',[0.83 0.88],... % Create line
                    [0.93 0.93]-posDec,...
                    'LineStyle','-',...
                    'Color',[1 0 0]);
                annotation(gcf,'textbox',...          % Create textbox
                    [0.88 0.91-posDec 0.07 0.06],...
                    'String',{'Peak Power'},...
                    'FontSize', legendFontSize, ...
                    'FitBoxToText','off',...
                    'LineStyle','none');

                % Friction
                if ~subtractFriction
                    posDec = posDec + 0.14;
                    annotation(gcf,'line',[0.83 0.88],... % Create line
                        [0.93 0.93]-posDec,...
                        'LineStyle',':',...
                        'Color',[0 0 0]);
                    annotation(gcf,'textbox',...          % Create textbox
                        [0.88 0.91-posDec 0.07 0.06],...
                        'String',{'Friction'},...
                        'FontSize', legendFontSize, ...
                        'FitBoxToText','off',...
                        'LineStyle','none');
                end
                % End of the customized legend
            end
        end
        
                
        function contSpeed = computeMaxContSpeed(this, torque)
            
            
            % shorthands
            t_r = this.jointModel.t_r;
            p_cm = this.jointModel.p_cm;
            
            % start off from the peak speed
            contSpeed = this.computeMaxPeakSpeed(torque);

            % rated torque
            contSpeed(torque > t_r) = 0;
            
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
            
            peakSpeed(torque > t_p) = 0;

            
        end

        
        function destFName = createDataSheet(this)
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
            this.compileTexFile;
            
            % Put the generated datasheet in the final destination.
            [~, fName] = fileparts(this.texFName);
            if ~exist(this.outputDir,'dir')
                mkdir(this.outputDir)
            end
            destFName = [this.outputDir,filesep, this.assembleOutFileName ];
            copyfile([fName,'.pdf'],destFName);

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
            h = this.drawTorqueSpeedCurve;
                        
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
            
            % Tested with MiKTeX/2.9 under Windows 7 x64 and
            % TexLive 2015.20160320-1 under Linux (Ubuntu 16.04)
            if (isunix)
                cmd = ['pdflatex -synctex=-1 -interaction=nonstopmode ', this.texFName];
            else
                cmd = ['lualatex.exe -synctex=-1 -interaction=nonstopmode ', this.texFName];
            end
            
            % Make the system call invoking the LuaLatex compiler.
            [flag, cmdout] = system(cmd);
                        
        end
        
        
    end
    
end

