%DATASHEETGENERATORTEST Test implementation of the datasheet generator.
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
%   results = runtests('dataSheetGeneratorTest.m')
%
% Alternatively, you can run tests using the run function.
%
%   results = run(dataSheetGeneratorTest)
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
% See also dataSheetGenerator.

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

function tests = dataSheetGeneratorTest
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

jb = jointBuilder;
jb.overwrite = 1;

motorName = 'dummyMotor';
dynName = 'full_dyn';
jb.buildJoint(motorName, dynName);
addpath(jb.buildDir)

testCase.('TestData').aJoint = eval([motorName,'_',dynName]);
testCase.('TestData').aJoint.name = 'Dummy Motor';

testCase.('TestData').jb = jb;

end

function teardownOnce(testCase)  % do not change function name
% change back to original path, for example

testCase.('TestData').jb.purge;

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
function testTorqueSpeedCurve(testCase)
% Test specific code

% shorthands
% allParams = testCase.('TestData').allParams;
% nPar = numel(allParams);
% 
% for iPar = 1:nPar
%    eval(allParams{iPar});
%    
%    dsg = dataSheetGenerator(params);
%    
%    dsg.createDataSheet;
%    
% end

dsg = dataSheetGenerator(testCase.('TestData').aJoint);

dsg.generateDataSheet;

delete( dsg.assembleOutFileName );
    
end