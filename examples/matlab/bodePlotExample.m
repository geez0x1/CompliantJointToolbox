% BODEPLOTEXAMPLE Exemplifies how to use the bode plot functionality
% shipped with CJT.
%
% The Compliant Joint Toolbox ships with tools to compute and display
% Bode plots from either simulated or experimental time domain data.
% This example provides a quick introduction on the provided
% functionalities on how to use them.
%
% Author::
%  Joern Malzahn
%  Wesley Roozing
%
% See also buildSomeJoints.

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

% The text in the following line defines the display name for this Example
% #! Bode Plot Example
clc;

disp('The Compliant Joint Toolbox ships with tools to compute and display')
disp('Bode plots from either simulated or experimental time domain data.')
disp('This example provides a quick introduction on the provided')
disp('functionalities on how to use them.')

pause
disp('')

%% Create input/output data
% Plant model
disp('The example uses a second order continuous time transfer function as a plant model.')
disp('The plant has a resonance frequency of 11 Hz and a damping ration of 0.5');
echo on

w0   = 2*pi*11;                      % resonance frequency [rad/s]
zeta = 0.5;                          % damping ratio
sys  = tf([w0^2], [1 2*zeta w0^2] ); % continuous time transfer function

echo off; disp(' ')
pause

% Data generation
disp('As test signal, an input chirp with a duration of t1 = 100 s is used. The chirp starts at f0 = 0.01 Hz and');
disp('linearly rises to f1 = 100 Hz. The test signal is sampled and applied at 1/dt = 2 kHz:');

echo on

dt  = 5e-4;     % sampling time
f0  = 0.01;     % initial frequency
f1  = 100;      % final frequency 
t   = 0:dt:100; % time vector

% chirp generation
u = chirp(t,f0,t(end),f1); % the chirp signal

% plant simulation
y = lsim(sys, u, t);       % the system response

echo off; disp(' ')

%% Visualize input/output data and bode diagram
% Plot input/output data
disp('First, the data is visualized to check if it is what we wanted.')
echo on;
figure(1)
clf;
hold on
plot(t(:),y(:),'k')
plot(t(:), u(:), 'r')
legend({'output', 'input'})
xlabel('Time [s]')
ylabel('Amplitude [.]')
echo off;
echo off; disp(' ')

% Bode diagram
disp('Configure the Bode diagram to appear in figure 2 with the x axis limits')
disp('set according to the range of our chirp signal and frequency units given in Hz.')
echo on;
figure(2)
clf;
hold on;
bodeOpt           = bodeoptions;     % configure the bode diagram appearance
bodeOpt.XLim      = [f0 f1];         % set the plot limits to minimum and maximum frequency, default would be [1, 10]
bodeOpt.FreqUnits = 'Hz';            % change frequency units, default frequency unit is rad/s
echo off; disp(' ')

disp('Now we display the actual Bode diagram for this data:');
pause()

echo on;
bode_tuyplot(t, ... time vector
             u, ... plant input
             y, ... plant output
             [],... resample, we don't do that here
             [],... filtering, we dont't do that here
             bodeOpt ); % configure plot appearance
echo off; disp(' ')
         
         
%% Bode diagram for imperfect data
disp('The data samples can be subject to jitter. That means that the sampling')
disp('intervals are not perfectly 0.5 ms. ')

disp('We simulate some jitter with 30% standard deviation in the sampling intervals')
pause()

echo on;
tj = abs(t + dt*0.3*randn(size(t)));  % jittered time instants
echo off;

disp('The jitter deteriorates the bode diagram quality visibly.')
pause()

echo on;
bode_tuyplot(tj, ...     time vector
             u, ...      plant input
             y, ...      plant output
             [],...      resample, we don't do that here, because we wish to see the effect of the jitter.
             [],...      filtering, we dont't do that here
             bodeOpt,... configure plot appearance
             'k');  ...  color this data
echo off;

disp('However, the the bode_tuyplot method can interpolate and resample the')
disp('data to improve this situation. ')
echo on;
bode_tuyplot(tj, ...     time vector
             u, ...      plant input
             y, ...      plant output
             1,...       do resample now, because we wish to see the improvements due to the resampling in the presence of the jitter 
             [],...      filtering, we dont't do that here
             bodeOpt,... configure plot appearance
             'g');  ...  color this data
echo off;

disp('The data samples are often subject to noise. The bode_tuyplot features a filtering functionality, to improve the bode diagram appearance.')
disp('Lets add 5 % zero mean Gaussian noise to the signal:')
echo on;
yn = y .* (1 + 0.05*randn(size(y)) );  % noisy signal
echo off;

disp('This noise has a clear effect on the bode diagram.')
echo on;
bode_tuyplot(t, ...     time vector
             u, ...      plant input
             yn, ...      plant output
             1,...       do resample now, to reduce the effect of jitter.
             [],...      filtering, we dont't do that here, we wish to see the impact of noise on the bode diagram.
             bodeOpt,... configure plot appearance
             'r');  ...  color this data
echo off;

disp('The built-in filtering function can smoothen the diagram to provide a better idea of how a diagram based on noise free data wouldlook like.')
echo on;
bode_tuyplot(t, ...     time vector
             u, ...      plant input
             yn, ...      plant output
             1,...      do resample now, to reduce the effect of jitter, filtering implies the use of the resampling anyways.
             1,...      now do filtering to improve the diagram
             bodeOpt,... configure plot appearance
             'b');  ...  color this data
echo off;