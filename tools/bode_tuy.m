% BODE_TUY Compute magnitude and phase in frequency domain from equidistantly 
% sampled I/O signals.
%
%   [f, mag_db, phase] = bode_tuy(t, u, y)
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
%
% Notes::
%
%
% Examples::
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also bode_tuyplot.

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


function [f, mag_db, phase] = bode_tuy(t, u, y)
    
    % Get timestep
    dt      = mean(diff(t));
    
    % Resample data based on timestep
    t_RS    = min(t):dt:max(t);
    u       = interp1(t, u, t_RS);
    y       = interp1(t, y, t_RS);
    
    % FFT
    Fs      = 1/dt;
    L       = length(t);
    NFFT    = 2^nextpow2(L);
    Y       = fft(y,NFFT)/L;
    f       = Fs/2*linspace(0,1,NFFT/2+1);
    U       = fft(u,NFFT)/L;
    H       = Y./U;
    
    % Calculate amplitude in decibels
    mag     = abs(H(1:NFFT/2+1));
    mag_db  = mag2db(mag);
    
    % Calculate phase in degrees and unwrap until the phase starts within
    % -90..+90 degrees.
    phase = (180/pi) * unwrap(angle(H(1:NFFT/2+1)));
    if (phase(1) > 90)
        while (phase(1) > 90)
            phase = phase - 180;    
        end
    elseif (phase(1) < -90)
        while (phase(1) < -90)
            phase = phase + 180;
        end
    end
    
end