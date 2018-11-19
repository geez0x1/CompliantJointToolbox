clc;
setCJTPaths;

%% Run Tests
results.JointBuilder        = runtests('cjtJointBuilderTest');
results.GenericJoint        = runtests('cjtGenericJointTest');
results.DataSheetGenerator  = runtests('cjtDataSheetGeneratorTest');
results.bodeTest            = runtests('cjtBodeTest');
results.bode2Test           = runtests('cjtBode2Test');
results.SimulinkLibrary     = runtests('cjtSimulinkLibraryTest');
results.ChecksumTest        = runtests('cjtChecksumTest');


%% Display Results
display('')
display('#####################################')
display('TEST SUMMARY')
display('#####################################')

nPassed = 0;
nFailed = 0;
nIncomplete = 0;

allFields = fields(results);
for iField = 1:numel(allFields)
    display('----------------')
    display(allFields{iField})
    display('----------------')
    curResult = results.(allFields{iField});
    
    nTest = numel(curResult);
    for iTest = 1:nTest
        nPassed = nPassed + curResult(iTest).Passed;
        nFailed = nFailed + curResult(iTest).Failed;
        nIncomplete = nIncomplete + curResult(iTest).Incomplete;
    end
    curResult.table
end
disp(['Passed:     ', num2str(nPassed)])
disp(['Failed:     ', num2str(nFailed)])
disp(['Incomplete: ', num2str(nIncomplete)])
