%CJTBODETEST Test implementation of the bode plotting functions.
%
% About this main function:
% The main function collects all of the local test functions
% into a test array. Since it is the main function, the function name
% corresponds to the name of your .m file and follows the naming convention
% of starting or ending in the word 'test', which is case-insensitive.
%
% To run tests from the command prompt, use the runtests command with your
% MATLAB test file as input. For example:
%
%   results = runtests('cjtBodeTest.m')
%
% Alternatively, you can run tests using the run function.
%
%   results = run(cjtBodeTest)
%
% To analyze the test results, examine the output structure from runtests
% or run. For each test, the result contains the name of the test function,
% whether it passed, failed, or didn't complete, and the time it took to
% run the test.
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

function tests = cjtBodeTest
    tests = functiontests(localfunctions);
end

% Setup and teardown code, also referred to as test fixture functions,
% set up the pretest state of the system and return it to the original
% state after running the test. There are two types of these functions:
% FILE FIXTURE functions that run once per test file, and FRESH FIXTURE
% functions that run before and after each local test function. These
% functions are not required to generate tests. In general, it is
% preferable to use fresh fixtures over file fixtures to increase unit test
% encapsulation.
%
% A function test case object, testCase, must be the only input to file
% fixture and fresh fixture functions. The Unit Test Framework
% automatically generates this object.
%
% The TestCase object is a means to pass information between setup
% functions, test functions, and teardown functions. Its TestData property
% is, by default, a struct, which allows easy addition of fields and data.
% Typical uses for this test data include paths and graphics handles.


% ----------------- FILE FIXTURE -----------------------
function setupOnce(testCase)  % do not change function name
% set a new path, for example

% cleanup console
close all;
clc;

% Create some data
dt = 5e-4;
f0 = 0.01;
f1 = 100;
t = 0:dt:100;
u = chirp(t,f0,t(end),f1);
w0 = 2*pi*11;
zeta = 0.5;
sys = tf([w0^2], [1 2*zeta w0^2] );
y = lsim(sys, u, t);

testCase.('TestData').t = t;
testCase.('TestData').u = u;
testCase.('TestData').y = y;
testCase.('TestData').f0 = f0;
testCase.('TestData').f1 = f1;
testCase.('TestData').dt = dt;

end

function teardownOnce(testCase)  % do not change function name
% change back to original path, for example



end

% % ----------------- FRESH FIXTURE -----------------------
% function setup(testCase)  % do not change function name
% % open a figure, for example
% end
% % 
% function teardown(testCase)  % do not change function name
% % close figure, for example
% end
% % -----------------------------------------------

% Individual test functions are included as local functions in the same
% MATLAB file as the main (test-generating) function. These test function
% names must begin or end with the case-insensitive word, 'test'. Each of
%the local test functions must accept a single input, which is a function
% test case object, testCase. The Unit Test Framework automatically
% generates this object.
%
% A test function is also called a "Qualification". There exist different
% conceptual types of qualifications.
function testBasicPlot(testCase)
% Test specific code

% shorthands
t = testCase.('TestData').t;
u = testCase.('TestData').u;
y = testCase.('TestData').y;
f0 = testCase.('TestData').f0;
f1 = testCase.('TestData').f1;

% DISPLAY SIGNAL
figure(1)
clf;
hold on
plot(t(:),y(:),'k')
plot(t(:), u(:), 'r')
legend({'output', 'input'})
xlabel('Time [s]')
ylabel('Amplitude [.]')

% BODEPLOT
figure(2)
clf;
bodeOpt             = bodeoptions;
opt                 = bode2options;
bodeOpt.XLim        = [f0 f1];          % set the plot limits to minimum and maximum frequency, default would be [1, 10]
bodeOpt.FreqUnits   = 'Hz';             % change frequency units, default frequency unit is rad/s

% actual bodeplot
bode_tuyplot(t, u, y, 0, 0, bodeOpt, opt);

verifyTrue(testCase,true) % If we arrive here, everything is fine.
    
end

function testResamplePlot(testCase)
% Test specific code

% shorthands
t = testCase.('TestData').t;
u = testCase.('TestData').u;
y = testCase.('TestData').y;
f0 = testCase.('TestData').f0;
f1 = testCase.('TestData').f1;
dt = testCase.('TestData').dt;

% DISPLAY SIGNAL
figure(1)
clf;
hold on
plot(t(:),y(:),'k')
plot(t(:), u(:), 'r')
legend({'output', 'input'})
xlabel('Time [s]')
ylabel('Amplitude [.]')

% Now we simulate some jitter in the sampling intervals
tj = abs(t + dt*0.3*randn(size(t)));  % jittered time instants

% BODEPLOT
figure(2)
clf;
hold on
bodeOpt             = bodeoptions;
opt                 = bode2options;
bodeOpt.XLim        = [f0 f1];      	% set the plot limits to minimum and maximum frequency, default would be [1, 10]
bodeOpt.FreqUnits   = 'Hz';             % change frequency units, default frequency unit is rad/s

% actual bodeplot
bode_tuyplot(t, u, y, [], [], bodeOpt, opt, 'r');  % result with the perfect signal
bode_tuyplot(tj, u, y, [], [], bodeOpt, opt, 'k'); % result with the jittered signal
bode_tuyplot(tj, u, y, 1, 0, bodeOpt, opt, 'g');   % result with the resampled jittered signal
bode_tuyplot(tj, u, y, 1, 1, bodeOpt, opt, 'b');   % reesult when the resampled jittered signal is filtered

legend({'ideal','jittered','resampled', 'filtered'},'location','SouthWest')

verifyTrue(testCase,true) % If we arrive here, everything is fine.
    
end