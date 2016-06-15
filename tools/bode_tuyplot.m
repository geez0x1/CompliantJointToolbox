% BODE_TUYPLOT Compute and plot magnitude and phase in frequency domain from equidistantly 
% sampled I/O signals.
%
%   [f, mag_db, phase] = bode_tuyplot(t, u, y [,lineseries properties])
%
% Inputs::
%   t: time vector
%   u: input data vector
%   y: output data vector
%
% Outputs::
%   f: frequency vector
%   mag_db: output magnitude vector in [db]
%   phase: output phase in [deg]
%
% Notes::
%   varargin holds additional plotting arguments passed to the plot command
%
% Examples::
%
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also bode_tuy.


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


function [f, mag_db, phase] = bode_tuyplot(t, u, y, varargin)
    
    % Get FFT
    [f, mag_db, phase] = bode_tuy(t, u, y);
    
    % Magnitude
    subplot(2,1,1); hold on;
    plot(f, mag_db,varargin{:});
    grid on
    ylabel('Magnitude [dB]');
    set(gca,'XScale','log');

    % Phase
    subplot(2,1,2); hold on;
    plot(f, phase,varargin{:});
    grid on;
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');
    set(gca,'XScale','log');
    
end