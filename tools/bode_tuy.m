%% [f, mag_db, phase] = bode_tuy(t, u, y)
% Compute magnitude and phase in frequency domain from equidistantly 
% sampled I/O signals.

function [f, mag_db, phase] = bode_tuy(t, u, y)
    
    % Get timestep
    dt      = mean(diff(t));
    
    % Resample data based on timestep
    t_RS	= min(t):dt:max(t);
    u       = interp1(t, u, t_RS);
    y       = interp1(t, y, t_RS);
    
    % FFT
    Fs      = 1/dt;
    L       = length(t);
    NFFT   	= 2^nextpow2(L);
    Y       = fft(y,NFFT)/L;
    f       = Fs/2*linspace(0,1,NFFT/2+1);
    U       = fft(u,NFFT)/L;
    H       = Y./U;
    
    % Calculate amplitude in decibels
    mag     = abs(H(1:NFFT/2+1));
    mag_db	= mag2db(mag);
    
    % Calculate phase in degrees and unwrap until the phase starts within
    % -90..+90 degrees.
    phase = unwrap((180/pi) * angle(H(1:NFFT/2+1)));
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