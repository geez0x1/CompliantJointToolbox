% GETMODEL_FROMDATA Estimate linear model from experimental data.
%
%   [ P ] = getModel_fromData(t, u, y, id_Np, id_Nz [, roi])
%
%  This function estimates a linear model (transfer function) from experi-
%  mental input/output data, with a specified number of zeroes and poles.
%  The timestep is determined from the median difference between the
%  elements of t.
%
% Inputs::
%   t: Time data vector
%   u: Input data vector
%   y: Output data vector
%   id_Np: Number of poles in model identification
%   id_Nz: Number of zeroes in model identification
%   roi: Frequency range of interest in [Hz], sets x limits in produced plots (default [0.1,100])
%
% Outputs::
%   P: Estimated transfer function
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
% See also getLinearDOB_fromData.

% Copyright (C) 2017, by Joern Malzahn, Wesley Roozing
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

function [ P ] = getModel_fromData(t, u, y, id_Np, id_Nz, roi)
    %% Default parameters
    if (~exist('id_Np', 'var') || isequal(id_Np,[]))
        id_Np   = 4;        % Model number of poles []
    end
    if (~exist('id_Nz', 'var') || isequal(id_Nz,[]))
        id_Nz   = 1;        % Model number of zeros []
    end
    if (~exist('roi', 'var') || isequal(roi,[]))
        roi = [0.1, 100];	% Region of interest [[Hz], [Hz]]
    end
    
    % Bode options
    bodeOpt = bodeoptions;
    bodeOpt.FreqUnits = 'Hz';


    %% Get variables

    % Resample data to obtain uniform sampling for tfest()
    Ts      = median(diff(t));      % Sampling time [s]
    t_RS    = t(1):Ts:t(end);       % Resampled time
    u       = interp1(t, u, t_RS)'; % Resampled input
    y       = interp1(t, y, t_RS)'; % Resampled output
    t       = t_RS';                % Replace time vector

    % Plot bode plot of original data
    [f, mag_db, phase] = bode_tuy(t, u, y);

    
    %% Identification

    % Generate iddata object of data
    d = iddata(y, u, [], 'SamplingInstants', t);

    % Identify transfer function P
    Options             = tfestOptions;
    Options.Display     = 'on';
    Options.InitMethod  = 'all';
    P                   = tfest(d, id_Np, id_Nz, Options);

    % Get magnitude and phase of Pc over f
    [mag_P, phase_P]    = bode(P, 2*pi*f);
    mag_db_P            = mag2db(mag_P(:));
    phase_P             = phase_P(:);


    %% Plot original data and P approximation
    figure();
    
    % Magnitude
    subplot(2,1,1);
    semilogx(f,mag_db); hold on;
    semilogx(f,mag_db_P, 'r');
    grid on
    xlim(roi);
    ylabel('Magnitude [dB]');
    legend('Data', 'Model');

    % Phase
    subplot(2,1,2);
    semilogx(f,phase); hold on;
    semilogx(f,phase_P, 'r');
    grid on;
    xlim(roi);
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');

    % Meta
    subplot(2,1,1);
    title(['Data vs approximation: Np = ' num2str(id_Np) ', Nz = ' num2str(id_Nz) ' (fit: ' num2str(P.Report.Fit.FitPercent, 3) '%)']);
    
end
