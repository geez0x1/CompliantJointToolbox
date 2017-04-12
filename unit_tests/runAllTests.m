clc;
setCJTPaths;

%% Run Tests
results_JointBuilder        = runtests('jointBuilderTest');
results_GenericJoint        = runtests('genericJointTest');
results_DataSheetGenerator  = runtests('dataSheetGeneratorTest');
results_bodeTest            = runtests('bodeTest');


%% Display Results
display('#####################################')
display('TEST SUMMARY')
display('#####################################')

summaryTable = ...
    [ results_JointBuilder.table;
      results_GenericJoint.table;
      results_DataSheetGenerator.table;
      results_bodeTest.table;
    ]
