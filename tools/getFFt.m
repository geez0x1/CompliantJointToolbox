function [f, amp] = getFFt(t, signal)
% GETFRUEQUENCIES Perform Fast Fourier Transform on a signal.
% =========================================================================
%
% [f, amp] = getfft(t, sig)
%
%  Description:
%    Time and signal vector must have equal length. Time steps are assumed
%    to be equidistant.
%
%  Inputs::
%    t:      Time vector
%    sig:    Signal vector
%
%  Outputs::
%    f1: vector of single-sided frequencies in Hz (0..fs/2)
%    f2: vector of two-sided frequencies in Hz (-fs/2..fs/2);
%
%
%  Authors:
%    Joern Malzahn
%    Wesley Roozing
%
%  See also fft, ifft.
%
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
%
    
    % compute sampling time, constant samping intervals are assumed
    Ts = mean(diff(t));

    % fft size
    l = length(signal);
    b=2^nextpow2(l);

    % two-sided spectrum
    F = fft(signal,b)/l;

    % create frequency vector up to Nyquist frequency
    f = 0.5/Ts*linspace(0,1,b/2);

    % single-sided amplitude spectrum
    amp = 2*abs(F(1:b/2));

end

