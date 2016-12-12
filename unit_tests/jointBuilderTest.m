%JOINTBUILDERTEST Test implementation of the joint builder.
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
%   results = runtests('jointBuilderTest.m')
%
% Alternatively, you can run tests using the run function.
%
%   results = run(jointBuilderTest)
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

function tests = jointBuilderTest
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

% instantiate a joint builder
testCase.('TestData').JB = jointBuilder;
testCase.('TestData').JB.overwrite = 1; % Overwrite existing model files.

testCase.('TestData').allParams = {'WMBig10k_ds';
    'WMBig10k';
    'WMBig2300_ds'};


end

function teardownOnce(testCase)  % do not change function name
% change back to original path, for example

testCase.('TestData').JB.purge; % remove all files created during the tests

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
function testNoFrictionModels(testCase)
% Test specific code

% shorthands
jb = testCase.('TestData').JB;
allParams = testCase.('TestData').allParams;
nPar = numel(allParams);

% tests
    for iPar = 1:nPar
        jb.buildJoint(allParams{iPar},'full_dyn_no_friction');
        jb.buildJoint(allParams{iPar},'output_fixed_no_friction');
        jb.buildJoint(allParams{iPar},'rigid_gearbox_no_friction');
        jb.buildJoint(allParams{iPar},'output_fixed_rigid_gearbox_no_friction');
        jb.buildJoint(allParams{iPar},'rigid_no_friction');
    end

    verifyTrue(testCase,true) % If we arrive here, everything is fine.
    
end

function testCoulombModels(testCase)
% Test to build models with coulomb friction.

% shorthands
jb = testCase.('TestData').JB;
allParams = testCase.('TestData').allParams;
nPar = numel(allParams);

% tests
    for iPar = 1:nPar
        jb.buildJoint(allParams{iPar}, 'full_dyn', 'coulomb');
        jb.buildJoint(allParams{iPar}, 'output_fixed', 'coulomb');
        jb.buildJoint(allParams{iPar}, 'rigid_gearbox', 'coulomb');
        jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox', 'coulomb');

        % Coulomb friction - asymmetric
        jb.buildJoint(allParams{iPar}, 'full_dyn', 'coulomb_asym');
        jb.buildJoint(allParams{iPar}, 'output_fixed', 'coulomb_asym');
        jb.buildJoint(allParams{iPar}, 'rigid_gearbox', 'coulomb_asym');
        jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox', 'coulomb_asym');
    end
    verifyTrue(testCase,true) % If we arrive here, everything is fine.
end

function testMultipleNonlinearTerms(testCase)
% Test to build models with combinations of nonlinear terms

% shorthands
jb = testCase.('TestData').JB;
allParams = testCase.('TestData').allParams;
nPar = numel(allParams);

% tests
    for iPar = 1:nPar
    
        jb.buildJoint(allParams{iPar}, 'full_dyn', {'coulomb', 'viscous_asym'});
        jb.buildJoint(allParams{iPar}, 'output_fixed', {'coulomb_asym', 'viscous_asym'});
        jb.buildJoint(allParams{iPar}, 'rigid_gearbox', {'coulomb', 'viscous_asym'});
        jb.buildJoint(allParams{iPar}, 'output_fixed_rigid_gearbox', {'coulomb_asym', 'viscous_asym'});
        
    end
    verifyTrue(testCase,true) % If we arrive here, everything is fine.
end

function testInputParameters(testCase)
    % Test to build models with correct and incorrect input arguments.

    % shorthands
    jb = testCase.('TestData').JB;
    allParams = testCase.('TestData').allParams;


    nTest = 6;              % We do a number nTest of tests here.
    flags = zeros(1,nTest);  % This vector

    % TESTS
    % this should work
    jb.buildJoint(allParams{1}, 'full_dyn', [],'electric_dyn');
    flags(1) = 1;

    % this should work
    jb.buildJoint(allParams{1}, 'full_dyn', [],'electric_dyn_zero_inductance');
    flags(2) = 1;

    % this should work
    jb.buildJoint(allParams{1}, 'full_dyn', [],[],'elTEST');
    flags(3) = 1;

    try
        % this should NOT work, electrical dynamics specified as nonlinear model
        jb.buildJoint(allParams{1}, 'full_dyn', 'electric_dyn_zero_inductance');
    catch
        flags(4) = 1;
    end

    try
        % this should NOT work, wrong electrical dynamics model string
        jb.buildJoint(allParams{1}, 'full_dyn', [], 'dyn_electric');
    catch
        flags(5) = 1;
    end

    try
        % this should NOT work, electrical dynamics model specified as linear mechanical model
        jb.buildJoint(allParams{1}, 'dyn_electric' );
    catch
        flags(6) = 1;
    end

    verifyTrue(testCase,all(flags)) % If we arrive here with all flags == 1, everything is fine.
end
    
function testElectricalSubsystem(testCase)
    
    % shorthands
    jb = testCase.('TestData').JB;
    allParams = testCase.('TestData').allParams;

    nTest = 4;              % We do a number nTest of tests here.
    flags = zeros(1,nTest);  % This vector should have just ones at the end of this test.
    
    % Build a joint
    jb.buildJoint(allParams{1}, 'full_dyn', [],'electric_dyn','elTestJoint');
    % Update path
    rehash
    % Instantiate joint
    jObj = elTestJoint;
    
    % Get electrical dynamics model
    [A B C D] = jObj.getElectricalDynamicsMatrices;

    % Check model against preevaluated values
    flags(1) = A == -6.321428571428571e+02;
    flags(2) = all( B ==    1.0e+03 * [  7.142857142857143   -0.004530000000000 ] );
    flags(3) = C ==    1;
    flags(4) = all( D ==   [0, 0] );
    
    verifyTrue(testCase,all(flags)) % If we arrive here with all flags == 1, everything is fine.
end