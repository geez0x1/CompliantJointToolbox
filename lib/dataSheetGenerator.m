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
%     dsg.generateDataSheet;                    
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

classdef dataSheetGenerator

    properties (SetAccess = private)
        jointModel                                     % The joint object to generate the datasheet for.
        templateDir                                    % Directory that contains all template files for the datasheet.
        
        texFName = 'cjtdsheet.tex';                    % Datasheet source filename.
        clsFName = 'cjtdsheet.cls';                    % Name of the style file for the datasheet.
        cfgFName = 'cjtdsheet.cfg';                    % Name of the temporary file containing the macro configuration for the datasheet sources.
        torqueSpeedFName = 'torqueSpeedCurve.pdf';     % Name of the temporary file that contains the torque speed curve of the actuator.
        efficiencyFName  = 'efficiencyCurve.pdf';      % Name of the temporary file that contains the efficiency plot of the actuator.
        thermalCharFName = 'thermalChar.pdf';          % Name of the temporary file that contains the curves with thermal operation ranges of the actuator.
        torFreqLoadFName = 'torFreqLoad.pdf';          % Name of the temporary file that contains the torque-frequency diagram for varying load inertia.
        torFreqLockFName = 'torFreqLock.pdf';          % Name of the temporary file that contains the torque-frequency diagram for varying locked output.
        
    end
    
    properties
        
        outFName = 'datasheet.pdf';                    % Name stub for the generated datasheet file. The final name will 
                                                       % be composed of the "jointModel.name_outFName".
        outputDir = ['.',filesep];                     % Directory to place the generated files. Defaults to the current directory.
        plotResolution = 600;                          % Print resolution for the figures contained in the datasheet.
        nPlotVals = 10000;                             % Number of samples used to produce the graphs in the datasheet figures.
        plotNormalized = 0;                            % Plot graphs in normalized quantities (default: false)
        freqMax = 100;                                 % Maximum frequency for spectral plots in Hz (default: 100 Hz)
        
        verbose = 0;                                   % Verbosity flag. 0 - (default) minimal console output, 1 - full console output. 
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
            % See also generateDataSheet, genericJoint, jointBuilder. 
            
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
            % See also generateDataSheet, genericJoint, jointBuilder.
            
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
            % See also generateDataSheet, genericJoint, jointBuilder.      
            
            % Shorthands
            t_p = this.jointModel.t_p;
            v_0 = this.jointModel.v_0;            
            k_t = this.jointModel.k_t;
            r   = this.jointModel.r;
            t_NL = this.jointModel.t_NL;
            i_NL = this.jointModel.i_NL;
            t_stall = this.jointModel.t_stall;
            N = this.jointModel.n;
            
            % Linearly spaced vector of torques between 0 and peak torque
            tauVals = (0:1/this.nPlotVals:1) * t_p;
            speedVals = this.torqueSpeedCurve(tauVals);
            
            % Power computation
            P_mech = speedVals .* (tauVals - i_NL*k_t*N);
            P_tot = v_0 / k_t/N * tauVals;

            % Efficiency computation
            eta = 100 * P_mech ./ P_tot;
            eta_max = 100 * ( 1 - sqrt( t_NL / t_stall)  )^2;
            
            % Initialize figure and plot
            h = gcf;
            [hAx,hLine1,hLine2] = plotyy(tauVals(:), eta(:), ... Efficiency plot 
                tauVals(:),[P_mech(:) P_tot(:)]);                % Power plot
            
            % Manipulate style of the efficiency plot 
            set(hLine1,'Color','r','linewidth',1.5)
            ylabel(hAx(1), 'Efficiency $\eta$ [\%]','Color','r','Interpreter','latex');
            set(hAx(1),'ylim',[0,100])
            set(hAx(1),'ytick',0:10:100)
            
            % Manipulate style of the power plot 
            set(hLine2,'Color','b','linewidth',1.5)
            set(hLine2(2),'LineStyle','--')
            ylabel(hAx(2),'Power $P$ [W]','Color','b','Interpreter','latex');
            
            set(hAx,'xlim',[0,t_p]);
            xlabel('Torque $\tau$ [Nm]','Interpreter','latex')
            title(['Maximum efficiency_ ', sprintf('%2.0f',eta_max),' %']);
            set(gca,'TickLabelInterpreter','latex');
            
            oldLim = get(hAx(2),'ylim');
            set(hAx(2),'ylim',[min(P_mech), oldLim(2)]);
            axes(hAx(1))
            hold on
            
            [~, idx_eta_max] = max(eta);
            stem(tauVals(idx_eta_max),eta_max,'or--');
            

            axes(hAx(2))
            hold on
            plot(tauVals,zeros(size(tauVals)),'--','color',0.5*[1 1 1])
            hAreaPTot = area(tauVals,P_tot(:),'FaceColor',0.8* [1 0 0],'FaceAlpha',0.25,'LineStyle','none');
            hAreaPMechWhite = area(tauVals,P_mech(:),'FaceColor', [1 1 1],'FaceAlpha',1,'LineStyle','none');
            hAreaPMech = area(tauVals,P_mech(:),'FaceColor',0.9* [0.2 0.2 1],'FaceAlpha',0.25,'LineStyle','none');
            plot(tauVals,P_mech,'b-','linewidth',1.5);
            hold off
            
           axes(hAx(1))
           hLegend = legend({'Efficiency [%]',...
               'Max. Efficiency',...
               'Power deliverd to load [W]',...
               'Total electrical power [W]'},...
               'location','southeast', 'box','off');
           
            
            
%             set(hLegend)
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
            % See also generateDataSheet, genericJoint, jointBuilder.      

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
            ylabel(hAx(1), 'Steady State Temp. [$^{\circ}$C]','Color','r','Interpreter','latex')
            set(hAx(1),'ylim',[0,1.05*max((Tmp_W))]);

             % Manipulate plot style
             set(hLine2,'Color','b','linewidth',1.5)
             ylabel(hAx(2),'Time to crit. temperature [s]','Color','b','Interpreter','latex')
             set(hAx(2),'ylim',[0,1.05*max(real(tCrit))]);

             set(hAx,'xlim',[0,t_p]);  
             
             legend({'Steady state temperature [$^{\circ}$C] ',...
                     'Critical winding temperature',...
                     'Copper melting temperature',...
                     'Time to critical winding temperature [s]'}, ...
                     'location','best',...
                     'Interpreter','latex')
                 
             xlabel('Torque $\tau$ [Nm]','Interpreter','latex')
             title(['Steady-state temperature rise from ',num2str(Tmp_ANom),' [${^\circ}$C]'],'Interpreter','latex')
             set(gca,'TickLabelInterpreter','latex');
             
             hold on;
             % Plot rated operation area
             Xr = [0 t_r t_r 0].';
             Yr = [0 0   Tmp_WMax Tmp_WMax].';
             rated_area = fill(Xr, Yr, 0.9* [0.2 0.2 1],'LineStyle','none');
             Xp = [t_r t_p t_p t_r].';
             Yp = [0 0  Tmp_WMax Tmp_WMax].';
             peak_area = fill(Xp, Yp, 0.8* [1 0 0],'LineStyle','none');
             
             alpha(0.25)     % add some transparency
             plot(hAx(1),[t_r,t_r],[0,TmeltCU],'color','b','linestyle','--')
             hold off;
             
             % Plot peak operation area
             axes(hAx(2))
             hold on
             alpha(0.25)     % add some transparency
             hold off
             
        end
                
        function h = drawTorqueFrequencyCurveLoad(this, fMax, plotNormalized)
            % DRAWTORQUEFREQUENCYCURVELOAD Theoretically feasible torque bandwidth under varying load inertia.
            %
            % Creates a contour plot that displays the -3dB cut-off frequency for the torque transfer function magnitude 
            % when the actuator generates peak torque. The contour lines are color coded on a logarithmic scale and 
            % correspond to different load inertiae in relation to the combined motor and gear inertia.
            %
            % A very light load will causes the actuator to reach the speed generation limit before a substantial torque
            % level can be established. A heavy load results in torques building up rapidly, with little actuator 
            % motion.
            %
            % The red curve corresponds to an infinite load inertia, that is equivalent to a locked output. The balck
            % solid curve corresponds to the load inertia being equal to the combined motor/gear inertia.
            % 
            % The plot includes the bandwidth limit due to the interplay of spring stiffness and back-EMF generation. 
            % The red area indicates the reachable operating conditions under the assumption of negligible spring 
            % damping. The zone expands with growing sping damping as indicated by black dotted curves. 
            %
            % One axis displays the final temperature at the given
            % operating condition. The second axis displays the time to
            % reach the maximum allowed temperature.
            %
            %   h = drawTorqueFrequencyCurveLoad(this, fMax, plotNormalized)
            %
            % Inputs:
            %   fMax           - The frequency in Hz up to which the magnitude is analyzed (default 100 Hz).
            %   plotNormalized - If true, the frequency axis is normalized to the frequency sqrt( k_ml/(I_m + I_g) ),
            %                    while the torque axis is normalized to the peak torque t_p (default 0). 
            %
            % Outputs:
            %   hContour: Handle to the contour plot
            %   
            %
            % Notes::
            %   For further reading, see: 
            %   J. Malzahn, N. Kashiri, W. Roozing, N. Tsagarakis and D. Caldwell, "What is the torque bandwidth of 
            %   this actuator?," 2017 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), 
            %   Vancouver, BC, 2017, pp. 4762-4768. doi: 10.1109/IROS.2017.8206351
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also generateDataSheet, drawTorqueFrequencyCurveLocked, drawTorqueSpeedCurve, genericJoint, jointBuilder.      
            
            % Check input parameters
            if ~exist('plotNormalized','var') % Normalize frequency and torque?
                plotNormalized = 0;
            end
            if ~exist('fMax','var') % Maximum frequency to plot (Hz)
                fMax = this.freqMax; 
            end
            
            % SHORTHANDS
            d_m  = this.jointModel.d_m;  % Viscous dampings
            d_g  = this.jointModel.d_g;
            d_l  = this.jointModel.d_l;
            d_gl = this.jointModel.d_gl;
            k_ml = this.jointModel.k_b;     % Spring stiffness
            t_p = this.jointModel.t_p;   % Peak torque
            dq_0 = this.jointModel.dq_0; % no torque speed
            I_m = this.jointModel.I_m;   % Motor and gear inertia
            I_g = this.jointModel.I_g;   % Motor and gear inertia
            nVals = this.nPlotVals;
            
            % Compute characteristic parameters
            I = I_m + I_g;
            w0 = sqrt(k_ml/I);                     % Resonance Frequency in rad/s
            dv = (d_m + d_g + d_l);             % Velocity dependent friction
            magDrop = db2mag(-3);%1/sqrt(2); % -3dB         % Allowed magnitude drop for tracking
            
            logScale = -3:0.15:1;
            I_l = I_m*10.^logScale;
            nInertias = numel(I_l);
            
            % NORMALIZATION
            if plotNormalized
                wNorm = 1/w0;
                fNorm = 2*pi*wNorm;
                tNorm = 1/t_p;
                xlabelStr = '$\tau_l / \tau_p$ [.]';
                ylabelStr = ' $f_c / f_0$ [.]';
            else
                fNorm = 1;
                tNorm = 1;
                xlabelStr = 'Torque $\tau$ [Nm]';
                ylabelStr = '3dB cut-off $f_c$ [Hz]';
            end
            
            xmax = 1.1 * t_p;
            wn = (1/nVals:1/nVals:1) * max(2*pi*fMax) / w0;
            
            w = wn*w0;
            
            % EVALUATE MAGNITUDES
            ZmagGain = zeros(nVals,nInertias);
            
            % N_l = (I_l(iI)*(1j*w).*(k_ml + d_gl*(1j*w)));
            % D_l = (I_m*I_l(iI)*(1j*w).^3 + (I_m*d_gl + I_l(iI)*d_gl + I_l(iI)*dv)*(1j*w).^2 + (I_m*k_ml + I_l(iI)*k_ml + d_gl*dv)*(1j*w) + k_ml*dv);
            % x:= w
            % y:= I_l
            N_l = @(x,y)(y.*(1j*x).*(k_ml + d_gl*(1j*x)));
            D_l = @(x,y)(I*y.*(1j*x).^3 + (I*d_gl + y.*d_gl + y.*dv)*(1j*x).^2 + (I*k_ml + y.*k_ml + d_gl*dv)*(1j*x) + k_ml*dv);
            tfMag = @(x,y) t_p*abs(N_l(x,y)./D_l(x,y));
            
            
            for iI = 1:nInertias
                ZmagGain(:,iI) = tfMag(w,I_l(iI));
            end
            
            
            figHandle = gcf;
            clf;
            hold on
                         
            % Plot the magnitude as a contour plot with colorbar legend
            [F,NORM_I] = meshgrid(w/2/pi*fNorm, log10(I_l/I));
            [~,h] = contour(tNorm/magDrop*ZmagGain.',F, (NORM_I), log10(I_l/I));
            
            myMap = flipud(winter(64));
            colormap(myMap)
            alpha(0.25)
            
            c = colorbar;%('YTick',log10(I_l/I));
            c.Label.Interpreter = 'latex';
            c.Label.String = '$\log_{10}(I_l / I)$ [.]';
            
            % A HACK TO GIVE THE COLOR BAR THE PROPER APPEARANCE
            % The color bar does not feature transparency anymore. So we put an empty textbox on top of it. The textbox
            % has the alpha property and we set it accordingly.
            hBox = annotation('textbox',...
                c.Position,...
                'FitBoxToText','off',...
                'FaceAlpha',0.25,...
                'EdgeColor',[1 1 1],...
                'BackgroundColor',[1 1 1]);

            % An issue arises, if we resize the figure. Then the textbox element would be in the wrong place. So we
            % register a callback to the figure. Whenever it resizes, it puts the textbox in the right spot.
                function sBarFcn(src,~)                  % Local function that implements the callback behavior.
                    txtBox = findall(src,'type','textboxshape');    % Locate the textbox handle
                    cBar = findall(gcf,'type','colorbar');          % Locate the colorbar handle
                    txtBox.Position = cBar.Position;                % Shift the textbox.
                end            
            set(gcf,'SizeChangedFcn',@sBarFcn);                     % Register handle to our callback function.
            % Hack completed...
            
            % PLOT A FEW SPECIAL MAGNITUDES
            hold on            
            
            % Load inertia equals motor+gear inertia
            magGain = tfMag(w,I);
            plot(magGain * tNorm/magDrop, w/2/pi*fNorm,'k','lineWidth',1.5)
            
            % Load inertia is quasi infinite
            magGain = tfMag(w,10^5*I);
            plot(magGain* tNorm/magDrop,w/2/pi*fNorm,'r','lineWidth',2)

            
            % LIMITATION DUE TO TORQUE SPEED CURVE
            % This is the bandwidth limit due to the motor back-emf generation           
            springTF = tf([0*d_gl, k_ml], [1 0]);
            springMag = bode(springTF,w);
            springMag = springMag(:);
            
            plot(dq_0*springMag *tNorm/ magDrop ,w/2/pi *fNorm,'k--','LineWidth',1.5)
            
            d =  2:2:20;
            
            for d_ = d
                magInt = abs( (1j*w*d_ + k_ml) ./ (1j*w) )  / magDrop;
                plot(dq_0*magInt * tNorm, w/2/pi * fNorm,':','color',0 * [ 1 1 1])
            end
            
            legendHandle = legend({'$I_l/I$', '$I_l/I = 1$', '$I_l/I \rightarrow \infty$', '$d_{gl} = 0$', '$d_{gl} > 0$'},'Interpreter','latex');
            set(legendHandle,'location','best')
            
            aHandle = area([0; dq_0*springMag] *tNorm/ magDrop , [fMax, w/2/pi]*fNorm );
            aHandle.LineStyle = 'none';
            aHandle.LineWidth = 1.5;
            aHandle.FaceColor = [0.8 0 0];
            aHandle.FaceAlpha = 0.15;
            aHandle.DisplayName = 'reachable ($d_{gl}=0$)';
            
            % PLOT APPEARANCE AND LABELS
            xlim([0,xmax]*tNorm);
            ylim([0,fMax ]*fNorm);
            xlabel(xlabelStr,'Interpreter','latex')
            ylabel(ylabelStr,'Interpreter','latex')
            title('Magnitude for harmonic peak torque operation with varying load inertia')
            box on;
            set(gca,'TickLabelInterpreter','latex');
        end
        
        function h = drawTorqueFrequencyCurveLocked(this, fMax, plotNormalized)
            % DRAWTORQUEFREQUENCYCURVELOAD Theoretically feasible torque bandwidth under locked actuator output.
            %
            % Creates a contour plot that displays the -3dB cut-off frequency for the torque transfer function magnitude 
            % and any torque amplitude from zero to peak torque. The color illustrates the thermally admissible 
            % operation time at a given frequency and torque amplitude, assuming the initial winding temperature being
            % normal temperature.
            %
            % Around the natural frequency f_0 (gray dashed line), theoretically any torque amplitude can be generated 
            % for an infinite amount of time.
            % Below and above the natural frequency, high torque amplitudes demand high motor currents that heat up the
            % motor windings. Beyond the rated torque (gray dotted lines), operating points can only be attained for
            % limited time. Operating points beyond the peak torque (black dotted line) are unreachable.
            % 
            % The plot includes the bandwidth limit due to the interplay of spring stiffness and back-EMF generation as 
            % black dashed line. The line rises towards higher frequencies with growing spring damping. 
            %
            % The line merges with the peak operation limit as indicated by the gray solid line. Operating points above 
            % this gray solid line are infeasible.
            %
            %
            %   h = drawTorqueFrequencyCurveLocked(this, fMax, plotNormalized)
            %
            % Inputs:
            %   fMax           - The frequency in Hz up to which the magnitude is analyzed (default 100 Hz).
            %   plotNormalized - If true, the frequency axis is normalized to the frequency sqrt( k_ml/(I_m + I_g) ),
            %                    while the torque axis is normalized to the peak torque t_p (default 0). 
            %
            % Outputs:
            %   hContour: Handle to the plot
            %   
            %
            % Notes::
            %   For further reading, see: 
            %   J. Malzahn, N. Kashiri, W. Roozing, N. Tsagarakis and D. Caldwell, "What is the torque bandwidth of 
            %   this actuator?," 2017 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), 
            %   Vancouver, BC, 2017, pp. 4762-4768. doi: 10.1109/IROS.2017.8206351
            %
            % Examples::
            %
            %
            % Author::
            %  Joern Malzahn
            %  Wesley Roozing
            %
            % See also generateDataSheet, drawTorqueFrequencyCurveLoad, drawTorqueSpeedCurve, genericJoint, jointBuilder. 
            
            % Check input parameters
            if ~exist('plotNormalized','var') % Normalize frequency and torque?
                plotNormalized = 0;
            end
            if ~exist('fMax','var') % Maximum frequency to plot (Hz)
                fMax = this.freqMax; 
            end
                        
            % SHORTHANDS
            d_m  = this.jointModel.d_m;       % Viscous dampings
            d_g  = this.jointModel.d_g;
            d_l  = this.jointModel.d_l;
            d_gl = this.jointModel.d_gl;
            I_m  = this.jointModel.I_m;
            k_ml = this.jointModel.k_b;       % Spring stiffness
            t_p = this.jointModel.t_p;        % Peak torque
            t_r = this.jointModel.t_r;        % Rated torque
            k_t = this.jointModel.k_t;        % Torque constant
            N = this.jointModel.n;            % Gear ratio
            dq_0 = this.jointModel.dq_0;      % No load speed
            r_th1 = this.jointModel.r_th1;    % Thermal resistance winding
            r_th2 = this.jointModel.r_th2;    % Thermal resistance housing
            r     = this.jointModel.r;        % Electrical winding resistance
            a_cu  =  0.0039;                  % Themperature coefficient copper
            T_max = this.jointModel.Tmp_WMax; % Max winding temperature
            t_w   = this.jointModel.T_thw;    % Thermal time constant widnings
            TAN   = this.jointModel.Tmp_ANom; % Normal ambient temperature.
            TmeltCU = 1084;                   % Melting temperature of copper
            
            nPoints = 100;
            
            % Compute characteristic parameters
            I = this.jointModel.I_m + this.jointModel.I_g; % Combined Inertia
            w0 = sqrt(k_ml/I);                             % Resonance Frequency in rad/s
            f0 = w0/2/pi;                       % Resonance Frequency in Hz
            dv = (d_m + d_g + d_l);             % Velocity dependent friction
            magDrop = db2mag(-3);               % -3dB, allowed magnitude drop for tracking
            
            % NORMALIZATION (by peak torque and natural frequency)
            if this.plotNormalized
                wNorm = 1/w0;
                fNorm = 2*pi*wNorm;
                tNorm = 1/t_p;
                xlabelStr = 'Torque $|\tau / \tau_p|$ [.]';
                ylabelStr = '3dB cut-off $f_c / f_0$ [.]';
            else
                fNorm = 1;
                tNorm = 1;
                xlabelStr = 'Torque $\tau$ [Nm]';
                ylabelStr = '3dB cut-off $f_c$ [Hz]';
            end
            
            % Linearly spaced vectors of torques and normalized frequencies
            ymax =  this.freqMax;
            xmax = 1.1 * t_p/magDrop;   
            tau = (0 : 1/nPoints : 1)*xmax;
            nLines = numel(tau);
            wn = (1/nPoints:1/nPoints:1) * 2*pi*ymax / w0;

            % BANDWIDTH LIMITATION DUE TO MAGNITUDE GAIN
            % N_l = (I_l(iI)*(1j*w).*(k_ml + d_gl*(1j*w)));
            % D_l = (I_m*I_l(iI)*(1j*w).^3 + (I_m*d_gl + I_l(iI)*d_gl + I_l(iI)*dv)*(1j*w).^2 + (I_m*k_ml + I_l(iI)*k_ml + d_gl*dv)*(1j*w) + k_ml*dv);
            % x:= w
            % y:= I_l
            N_l = @(x,y)(y.*(1j*x).*(k_ml + d_gl*(1j*x)));
            D_l = @(x,y)(I*y.*(1j*x).^3 + (I*d_gl + y.*d_gl + y.*dv)*(1j*x).^2 + (I*k_ml + y.*k_ml + d_gl*dv)*(1j*x) + k_ml*dv);
            tfMag = @(x,y) abs(N_l(x,y)./D_l(x,y));
            magGain = tfMag(wn*w0,1e9).'; % 1e9 ~ infinite impedance reflects locked output
            
           
            % LIMITATION DUE TO TORQUE SPEED CURVE
            % This is the bandwidth limit due to the motor back-emf
            % generation and limited supply voltage
            % Torque sensor transfer function
            springTF = tf([d_gl, k_ml], [1 0]);
            springMag = bode(springTF,wn*w0);
            springMag = springMag(:);
            
            % TOTAL FREQUENCY LIMIT          
            fSpring = dq_0*springMag;
            fAmp = magGain*t_p;
            fIdx = find(fSpring < fAmp,1,'last');
            fTot = [fSpring(1:fIdx); fAmp(fIdx+1:end)];
            
            
            % THERMAL TIME LIMITATION
            % Steady State Temperature
            resCoeff = (r_th1 + r_th2) * r / N^2 / k_t^2; % thermal resistance coefficient
            Tmp_W = TAN + ...  
                (  ( resCoeff * tau.^2 ) ./ ( 1 - a_cu * resCoeff *  min(tau,t_r).^2 )  );
            
            % Saturate temperature at melting point
            idx = find(Tmp_W > TmeltCU,1,'first');
            Tmp_W(idx:end) = TmeltCU;

            % Time to reach max winding temperature from normal conditions
            Tmax = this.jointModel.Tmp_WMax;
            T0 = this.jointModel.Tmp_ANom;
            dTend = Tmp_W - T0;
            tCrit = - this.jointModel.T_thw * log( - (Tmax - T0  - dTend ) ./ dTend);
            
            % Within the rated operation, the time to reach the maximum allowed temperature is infinite.
            idx = find(tau > t_r,1,'first');
            tCrit(1:idx) = inf;
            
              
             %% PLOTTING
            figHandle = gcf;
            clf;
            hold on

            % Contour plot
            %   Prepare grid data
            Z = repmat(magGain,1,nLines) .* tau;
            TMAX = repmat(magGain,1,nLines).*tCrit;
            [F,TAU] = meshgrid(wn*w0/2/pi*fNorm, tau*tNorm  /magDrop);

            
            %   Create plot
             h = surf(Z,F.', (TAU.'),min(TMAX,0.5*t_w),'DisplayName','t_{max} [s] to max. temperature');

            %   Configure plot
            set(h,'LineStyle','none') % no grid lines
            view(0,-90)               % since we are using surf with additionally plotted lines, we should look from the bottom
            set(gca,'YDir','reverse') % reverse frequency axis to compensate for the view
            set(gcf,'renderer','opengl')
            
            %   adjust colormap from cold (blue) to hot (red)
            nCVals = 64;    
            myMap = [(nCVals:-1:0).', 0*ones(nCVals+1,1), (0:nCVals).']/nCVals;
            colormap((myMap))
            %   add a color bar legend
            c = colorbar;      
            % Saturate the time to for rated operation
            tmp = c.TickLabels;
            tmp{end} = '\infty';
            c.TickLabels = tmp;
            c.Label.Interpreter = 'latex';
            c.Label.String = '$t_{max}$ [s]';
                         
            
             
            % Plot total frequency limit (back-emf/voltage and sensor)
            plot(fTot*tNorm/magDrop, wn*w0/2/pi * fNorm,'Color', 0.75* [1 1 1 ],'LineWidth',5,'DisplayName','feasibility limit')
            
            % Plot torque transfer function magnitude
            plot((magGain*t_p)*tNorm/magDrop,    wn*w0/2/pi * fNorm,':','Color',0*[1 1 1],'LineWidth',1.5,'DisplayName', 'peak amplitude \tau_p') % -3 dB
            plot((magGain*t_r)*tNorm/magDrop,    wn*w0/2/pi * fNorm,':','Color',0.75*[1 1 1],'LineWidth',1.5,'DisplayName', 'rated amplitude \tau_r') % -3 dB
            
            % Plot voltage saturation limit
            plot(dq_0*springMag *tNorm/ magDrop , wn*w0/2/pi *fNorm,'k--','LineWidth',1.5,'DisplayName','back-EMF limit')
            
            % plot resonance frequency
            plot(tau*tNorm, f0*ones(size(tau))*fNorm, '--','color',0.8*[1 1 1],'LineWidth',1.0,'DisplayName','f_0 [Hz]')

            %% Plot appearance and labels
            xlim([0,xmax]*tNorm);
            ylim([ wn(1)*w0/2/pi ,ymax]*fNorm);
            xlabel(xlabelStr,'Interpreter','latex')
            ylabel(ylabelStr,'Interpreter','latex')
            box on;
            set(gca,'TickLabelInterpreter','latex');
            legend('location','best')
            title('Bandwidth for locked output with thermally permissible operation time.')
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
            % See also generateDataSheet, drawTorqueFrequencyCurveLoad, genericJoint, jointBuilder.
            
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
                xlabelStr = 'Torque $\tau / \tau_p$ [.]';
                ylabelStr = 'Speed $\dot{q} / \dot{q}_{NL}$ [.]';
            else
                tNorm = 1;
                dqNorm = 1;
                xlabelStr = 'Torque $\tau$ [Nm]';
                ylabelStr = 'Speed $\dot{q}$ [rad/s]';
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
            title('Torque-speed diagram.')
            
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

        
        function destFName = generateDataSheet(this)
            % GENERATEDATASHEET Main function to trigger the data sheet
            % generation.
            %
            %   dsg.generateDataSheet
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
            % See also genericJoint, jointBuilder.            
            
            % Get local copies of all template files.
            copyfile([this.templateDir, filesep, this.texFName],this.texFName);
            copyfile([this.templateDir, filesep, this.clsFName],this.clsFName);
            
            % Collect all information about the actuator and prepare macros
            % for the use in the Tex source.
            this.createDefFile;

            % Produce all plots
            this.makeDataSheetPlots;
            
            % Make sure, the output directory exists
            if ~exist(this.outputDir,'dir')
                mkdir(this.outputDir)
            end
            
            % Compile the Tex source. 
            [flag, cmdout] = this.compileTexFile;

            % Ouptut last console output, if desired.
            if this.verbose
                disp(cmdout);
            end
            
            % Put the generated datasheet in the final destination.
            [~, fName] = fileparts(this.texFName);

            destFName = [this.outputDir,filesep, this.assembleOutFileName ];
            srcFName =  [this.outputDir,filesep, fName, '.pdf'];
            copyfile(srcFName,destFName);

            % Clean up.
            delete([fName,'.*'])
            delete([this.torqueSpeedFName])
            delete([this.efficiencyFName])
            delete([this.thermalCharFName])
            delete([this.torFreqLoadFName])
            delete([this.torFreqLockFName])
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
            % See also generateDataSheet, genericJoint, jointBuilder.   
            
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
            % See also generateDataSheet, genericJoint, jointBuilder.    
            
            close all;
            
            % Torque-Speed curve
            h = this.drawTorqueSpeedCurve;
                      
            % In the datasheet we are going to have a figure caption. Remove the axes title to save space.
            ax = get(gcf,'children');
            set(ax,'Title',[])
            
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
            
            % Efficiency curve
            h = this.drawEfficiencyCurve;
            
            set(gcf,'Units','centimeters');
            set(gcf,'PaperUnits','centimeters');
            pos = get(gcf,'Position');
            pos(3) = 18;
            pos(4) = 10.5;
             
            set(gcf,'Position',pos)
            set(gcf,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])          
            set(gca,'LooseInset',get(gca,'TightInset'))
            
            printpdf(gcf,this.efficiencyFName,['-r',num2str(this.plotResolution)])

            close(gcf);
             
            % Thermal characteristics
            h = this.drawThermalCharacteristics;

            % In the datasheet we are going to have a figure caption. Remove the axes title to save space.
            title('') 
            
            set(gcf,'Units','centimeters');
            set(gcf,'PaperUnits','centimeters');
            pos = get(gcf,'Position');
            pos(3) = 18;
            pos(4) = 8;
             
            set(gcf,'Position',pos)
            set(gcf,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])          
            set(gca,'LooseInset',get(gca,'TightInset'))
            
            printpdf(gcf,this.thermalCharFName,['-r',num2str(this.plotResolution)])

            close(gcf);
            
            % Torque Frequency Characteristics varying load
            h = this.drawTorqueFrequencyCurveLoad;

            % In the datasheet we are going to have a figure caption. Remove the axes title to save space.
            title('') 
            
            set(gcf,'Units','centimeters');
            set(gcf,'PaperUnits','centimeters');
            pos = get(gcf,'Position');
            pos(3) = 18;
            pos(4) = 9;
             
            set(gcf,'Position',pos)
            set(gcf,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])          
            set(gca,'LooseInset',get(gca,'TightInset'))
            
            hCbar = findobj(gcf,'Type','colorbar');
            hCbar.AxisLocation = 'out';
            
            printpdf(gcf,this.torFreqLoadFName,['-r',num2str(this.plotResolution)])

            close(gcf);
            
            % Torque Frequency Characteristics varying load
            h = this.drawTorqueFrequencyCurveLocked;

            % In the datasheet we are going to have a figure caption. Remove the axes title to save space.
            title('') 
            
            set(gcf,'Units','centimeters');
            set(gcf,'PaperUnits','centimeters');
            pos = get(gcf,'Position');
            pos(3) = 18;
            pos(4) = 9;
             
            set(gcf,'Position',pos)
            set(gcf,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])          
            set(gca,'LooseInset',get(gca,'TightInset'))
            
            hCbar = findobj(gcf,'Type','colorbar');
            hCbar.FontSize = get(gca,'FontSize');
            hCbar.AxisLocation = 'in';
                        
            printpdf(gcf,this.torFreqLockFName,['-r',num2str(this.plotResolution)])

            close(gcf);
            
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
            % See also generateDataSheet, genericJoint, jointBuilder.    
            %
            
            jM = this.jointModel; % Shorthand
            
            fid = fopen(this.cfgFName,'w+');
            %
            fprintf(fid,'\\def \\valJointName{%s}\n'                 , strrep(jM.name,'_',' ')          );
            %
            % File names of plots to include
            %
            fprintf(fid,'\\def \\valTorqueSpeedFName{%s}\n'           , this.torqueSpeedFName           );
            fprintf(fid,'\\def \\valEfficiencyFName{%s}\n'            , this.efficiencyFName            );
            fprintf(fid,'\\def \\valThermalCharFName{%s}\n'           , this.thermalCharFName           );
            fprintf(fid,'\\def \\valTorFreqLoadFName{%s}\n'           , this.torFreqLoadFName           );
            fprintf(fid,'\\def \\valTorFreqLockFName{%s}\n'           , this.torFreqLockFName           );       
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
            % See also system, generateDataSheet, genericJoint, jointBuilder.    
            %
            
            % Tested with MiKTeX/2.9 under Windows 7 x64 and
            % TexLive 2015.20160320-1 under Linux (Ubuntu 16.04)
            if (isunix)
%                 cmd = ['pdflatex -synctex=-1 -interaction=nonstopmode ', this.texFName];
                  cmd = sprintf('pdflatex -synctex=-1 -interaction=nonstopmode -output-directory %s %s',...
                      this.outputDir, this.texFName);
            else
%                 cmd = ['lualatex.exe -synctex=-1 -interaction=nonstopmode ', this.texFName];
                cmd = sprintf('lualatex.exe -synctex=-1 -interaction=nonstopmode -output-directory=%s %s',...
                      this.outputDir, this.texFName);
            end
            
            % Make the system call invoking the \Latex compiler.
            % We compile three times to get all references right.
            cmdout = [];
            flag = 0;
            for it = 1:3 
            [flag, cmdout] = system(cmd);
                if flag ~= 0
                    disp('Datasheet compilation failed. See compiler output for details: ')
                    disp(cmdout)
                    
                    if (isunix)
                        disp('-------');
                        disp('Warning: If the above error references libstdc++, this issue arises due to MATLAB including libraries which are too old for your system. Try removing the following symlinks:');
                        disp('- /MATLAB/sys/os/glnxa64/libstdc++.so.6');
                        disp('- /MATLAB/bin/glnxa64/libtiff.so.5');
                        disp('See also: https://mathworks.com/matlabcentral/answers/329796-issue-with-libstdc-so-6');
                        disp('-------');
                    end
                end
            end
            
            % Ouptut last console output, if desired.
            if this.verbose
                disp(cmdout);
            end
                                    
        end
        
        
    end
    
end

