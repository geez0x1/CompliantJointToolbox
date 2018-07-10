%CJTCHECKSUMTEST Test implementation of the cjtChecksumComputer
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
%   results = runtests('cjtChecksumTest.m')
%
% Alternatively, you can run tests using the run function.
%
%   results = run(cjtChecksumTest)
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

function tests = cjtChecksumTest
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

testCase.('TestData').LoremIpsum = ...
    strcat('Lorem ipsum dolor sit amet, consetetur sadipscing elitr,',...
           'sed diam nonumy eirmod tempor invidunt ut labore et dolore',...
           'magna aliquyam erat, sed diam voluptua. At vero eos et ',...
           'accusam et justo duo dolores et ea rebum. Stet clita kasd',...
           'gubergren, no sea takimata sanctus est Lorem ipsum dolor sit',...
           'amet. Lorem ipsum dolor sit amet, consetetur sadipscing elitr,',...
           'sed diam nonumy eirmod tempor invidunt ut labore et dolore',...
           'magna aliquyam erat, sed diam voluptua. At vero eos et accusam',...
           'et justo duo dolores et ea rebum. Stet clita kasd gubergren,',...
           'no sea takimata sanctus est Lorem ipsum dolor sit amet.');

testCase.('TestData').Magic10 = magic(10);
       
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
function testLoremIpsumArray(testCase)
% This test verifies that the checksums of character arrays are compared correctly


    res_1 = cjtCompareChecksum(testCase.('TestData').LoremIpsum,...
        testCase.('TestData').LoremIpsum);                            % self-check
    res_2 = cjtCompareChecksum(testCase.('TestData').LoremIpsum, ...
                                        'CEE2BD6F0C8368B4783A9B3B7F905951'); 
    res_3 = cjtCompareChecksum(testCase.('TestData').LoremIpsum, ...
                                        'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'); 
    
    verifyTrue(testCase, res_1);    % Self-check
    verifyTrue(testCase, res_2);    % Matching checksum
    verifyTrue(testCase, ~res_3);   % Check of cjtCompareChecksum (must not always return true)
    
end

function testLoremIpsumFile(testCase)
% This test verifies that the checksums of text files are compared correctly

    LoremIpsum = testCase.('TestData').LoremIpsum;
    fName = 'tmp_cjt_checksum_test.txt';
    save(fName,...
        'LoremIpsum',...
        '-ascii');

    res_1 = cjtCompareChecksum(fName, fName); % self-check
    res_2 = cjtCompareChecksum(fName, '804F52A6C2563816396963B5C8B88F1E'); 
    res_3 = cjtCompareChecksum(fName, '804F52A6C2563816396963B5C8B88FE1'); 

    delete(fName);

    verifyTrue(testCase,res_1 && res_2 && ~res_3); % Checksums must always match
    
end

function testMagic10Array(testCase)
% This test verifies that the identical checksums of numeric arrays are compared correctly
    
    mysum = cjtComputeChecksum(testCase.('TestData').Magic10);          % compute checksum
    res = cjtCompareChecksum(testCase.('TestData').Magic10, mysum);     % self-check
    
    verifyTrue(testCase,res);                                           % validate self-check
    verifyEqual(testCase,mysum, '2FEC13058FA99395A58FB4C0FE3F4BE9');    % validate reference checksum
    
end

function testMagic10File(testCase)
% This test verifies that the checksums of mat files are compared correctly

    magic10 = testCase.('TestData').Magic10;
    fName = 'tmp_cjt_checksum_test.mat';
    save(fName,'magic10');

    res_1 = cjtCompareChecksum(fName, fName); % self-check
    % For binary files, the checksum also includes the file creation date.
    % That makes it hard to define a reference checksum.

    delete(fName);

    verifyTrue(testCase,res_1); % Checksums must always match
    
end
